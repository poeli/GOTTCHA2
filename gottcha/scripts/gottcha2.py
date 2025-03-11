#!/usr/bin/env python3

__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__credits__   = ["Po-E Li", "Anna Chernikov", "Jason Gans", "Tracey Freites", "Patrick Chain"]
__version__   = "2.1.8.11"
__date__      = "2018/10/07"
__copyright__ = """
Copyright (2019). Traid National Security, LLC. This material was produced
under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Laboratory
(LANL), which is operated by Los Alamos National Security, LLC for the U.S.
Department of Energy. The U.S. Government has rights to use, reproduce, and
distribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE
USE OF THIS SOFTWARE.  If software is modified to produce derivative works, such
modified software should be clearly marked, so as not to confuse it with the
version available from LANL.

Additionally, this program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free
Software Foundation; either version 3 of the License, or (at your option) any
later version. Accordingly, this program is distributed in the hope that it will
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
License for more details.
"""

import argparse as ap, textwrap as tw
import sys, os, time, subprocess
import pandas as pd
import gc
from re import search,findall
from multiprocessing import Pool
from itertools import chain
import math
import logging

try:
    # Try relative import first (for package usage)
    from . import taxonomy as gt
except ImportError:
    # Fall back to direct import (for script usage)
    import taxonomy as gt

def parse_params(ver, args):
    """
    Parse and validate command line arguments for GOTTCHA2.
    
    This function sets up the argument parser, defines all possible command-line
    options, parses the provided arguments, and performs validation to ensure
    the configuration is valid and complete.
    
    Parameters:
        ver (str): Version string to display in help messages
        args (list): Command line arguments to parse
        
    Returns:
        argparse.Namespace: Object containing all validated arguments
        
    Raises:
        SystemExit: If validation fails or --version is specified
    """
    p = ap.ArgumentParser( prog='gottcha2.py', description="""Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an
            annotation-independent and signature-based metagenomic taxonomic profiling tool
            that has significantly smaller FDR than other profiling tools. This program
            is a wrapper to map input reads to pre-computed signature databases using minimap2
            and/or to profile mapped reads in SAM format. (VERSION: %s)""" % ver)

    eg = p.add_mutually_exclusive_group( required=True )

    eg.add_argument( '-i','--input', metavar='[FASTQ]', nargs='+', type=ap.FileType('r'),
                    help="Input one or multiple FASTQ/FASTA file(s). Use space to separate multiple input files.")

    eg.add_argument( '-s','--sam', metavar='[SAMFILE]', nargs=1, type=ap.FileType('r'),
                    help="Specify the input SAM file. Use '-' for standard input.")

    p.add_argument( '-d','--database', metavar='[GOTTCHA2_db]', type=str, default=None,
                    help="The path and prefix of the GOTTCHA2 database.")

    p.add_argument( '-l','--dbLevel', metavar='[LEVEL]', type=str, default='',
                    choices=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'],
                    help="""Specify the taxonomic level of the input database. You can choose one rank from "superkingdom", "phylum", "class", "order", "family", "genus", "species" and "strain". The value will be auto-detected if the input database ended with levels (e.g. GOTTCHA_db.species).""")

    p.add_argument( '-ti','--taxInfo', metavar='[PATH]', type=str, default='',
                    help="""Specify the path to the taxonomy information directory or file. The program will attempt to locate a matching .tax.tsv file for the specified database. If it cannot find one, it will use the ‘taxonomy_db’ directory located in the same directory as the executable by default.""")

    p.add_argument( '-np','--nanopore', action="store_true",
                    help="Adjust options for Nanopore reads. It will overwrite the other options to '-xm map-ont -mr 1 -mf 0 -mz 0'.")

    p.add_argument( '-pm','--mismatch', metavar='<INT>', type=int, default=10,
                    help="Mismatch penalty for the aligner. [default: 10]")

    p.add_argument('-e', '--extract', metavar='[TAXID[,TAXID...] or @FILE]', type=str, default=None,
                help="""Extract reads mapping to specific TAXIDs. Either comma-separated list 
                or @filename with one taxid per line. Examples: 
                -e "1234,5678" or -e @taxids.txt [default: None]""")

    p.add_argument( '-fm','--format', metavar='[STR]', type=str, default='tsv',
                    choices=['tsv','csv','biom'],
                    help='Format of the results; available options include tsv, csv or biom. [default: tsv]')

    p.add_argument( '-r','--relAbu', metavar='[FIELD]', type=str, default='ROLLUP_DOC',
                    choices=['ROLLUP_DOC','READ_COUNT','TOTAL_BP_MAPPED'],
                    help='The field will be used to calculate relative abundance. You can specify one of the following fields: "LINEAR_LENGTH", "TOTAL_BP_MAPPED", "READ_COUNT" and "LINEAR_DOC". [default: ROLLUP_DOC]')

    p.add_argument( '-t','--threads', metavar='<INT>', type=int, default=1,
                    help="Number of threads [default: 1]")

    p.add_argument( '-o','--outdir', metavar='[DIR]', type=str, default='.',
                    help="Output directory [default: .]")

    p.add_argument( '-p','--prefix', metavar='<STR>', type=str, required=False,
                    help="Prefix of the output file [default: <INPUT_FILE_PREFIX>]")

    p.add_argument( '-xm','--presetx', metavar='<STR>', type=str, required=False, default='sr',
                    choices=['sr','map-pb','map-ont'],
                    help="The preset option (-x) for minimap2. Default value 'sr' for short reads. [default: sr]")

    p.add_argument( '-mc','--minCov', metavar='<FLOAT>', type=float, default=0.005,
                    help="Minimum linear coverage to be considered valid in abundance calculation. [default: 0.005]")

    p.add_argument( '-mr','--minReads', metavar='<INT>', type=int, default=3,
                    help="Minimum number of reads to be considered valid in abundance calculation. [default: 3]")

    p.add_argument( '-ml','--minLen', metavar='<INT>', type=int, default=60,
                    help="Minimum unique length to be considered valid in abundance calculation. [default: 60]")

    p.add_argument( '-mz','--maxZscore', metavar='<FLOAT>', type=float, default=30,
                    help="Maximum estimated z-score for the depths of the mapped region. Set to 0 to disable. [default: 30]")

    p.add_argument( '-mf','--matchFactor', metavar='<FLOAT>', type=float, default=0.5,
                    help="Minimum fraction of the read or signature fragment required to be considered a valid match. [default: 0.5]")

    p.add_argument( '-nc','--noCutoff', action="store_true",
                    help="Remove all cutoffs. This option is equivalent to use. [-mc 0 -mr 0 -ml 0 -mf 0 -mz 0]")

    p.add_argument( '-sm','--skipRemoveMultiple', action="store_true",
                    help="This option can be use to skip the step to removal of multiple hits If you are sure there are no multiple hits for the same reads in the SAM file.")

    p.add_argument( '-c','--stdout', action="store_true",
                    help="Write on standard output.")

    eg.add_argument( '-v','--version', action="store_true",
                    help="Print version number.")

    p.add_argument( '--silent', action="store_true",
                    help="Disable all messages.")

    p.add_argument( '--verbose', action="store_true",
                    help="Provide verbose messages.")

    p.add_argument( '--debug', action="store_true",
                    help="Debug mode. Provide verbose running messages and keep all temporary files.")

    args_parsed = p.parse_args(args)

    """
    Checking options
    """
    if args_parsed.version:
        print( ver )
        sys.exit(0)

    if not args_parsed.database:
        p.error( '--database option is missing.' )

    if args_parsed.input and args_parsed.sam:
        p.error( '--input and --same are incompatible options.' )

    if args_parsed.database:
        #assign default path for database name
        if "/" not in args_parsed.database and not os.path.isfile( args_parsed.database + ".mmi" ):
            bin_dir = os.path.dirname(os.path.realpath(__file__))
            args_parsed.database = bin_dir + "/database/" + args_parsed.database

    if args_parsed.database and args_parsed.database.endswith(".mmi"):
        args_parsed.database.replace('.mmi','')

    if args_parsed.database and args_parsed.input:
        if not os.path.isfile( args_parsed.database + ".mmi" ):
            p.error( 'Database index %s.mmi not found.' % args_parsed.database )

    if not args_parsed.taxInfo:
        if args_parsed.database:
            db_dir = search(r'^(.*?)[^\/]+$', args_parsed.database )
            args_parsed.taxInfo = db_dir.group(1) + "/taxonomy_db"

            if not os.path.isdir(args_parsed.taxInfo):
                bin_dir = os.path.dirname(os.path.realpath(__file__))
                args_parsed.taxInfo = bin_dir + "/database"

            if not os.path.isdir(args_parsed.taxInfo):
                args_parsed.taxInfo = db_dir.group(1)

    if not args_parsed.prefix:
        if args_parsed.input:
            name = search(r'([^\/\.]+)\..*$', args_parsed.input[0].name )
            args_parsed.prefix = name.group(1)
        elif args_parsed.sam:
            name = search(r'([^\/]+).\w+.\w+$', args_parsed.sam[0].name )
            args_parsed.prefix = name.group(1)
        else:
            args_parsed.prefix = "GOTTCHA_"

    if not args_parsed.dbLevel:
        if args_parsed.database:
            major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7, "strain":8}
            parts = args_parsed.database.split('.')
            for part in parts:
                if part in major_ranks:
                    args_parsed.dbLevel = part
                    break
        elif args_parsed.sam:
            name = search(r'\.gottcha_(\w+).sam$', args_parsed.sam[0].name )
            try:
                args_parsed.dbLevel = name.group(1)
            except:
                pass
        
        if not args_parsed.dbLevel:
            p.error( '--dbLevel is missing and cannot be auto-detected.' )

    if args_parsed.noCutoff:
        args_parsed.minCov = 0
        args_parsed.minReads = 0
        args_parsed.minLen = 0
        args_parsed.matchFactor = 0
        args_parsed.maxZscore = 0

    if args_parsed.nanopore:
        args_parsed.presetx = 'map-ont'
        args_parsed.minReads = 1
        args_parsed.matchFactor = 0
        args_parsed.maxZscore = 0

    return args_parsed

def dependency_check(cmd):
    """
    Verify that external dependencies are available in the system.
    
    Attempts to execute the specified command with --help and checks if it runs
    successfully. Exits the program if the command is not found or fails.
    
    Parameters:
        cmd (str): Command to check
        
    Returns:
        None
        
    Raises:
        SystemExit: If the command is not found or fails
    """
    try:
        subprocess.check_call([cmd, "--help"], stdout=subprocess.DEVNULL)
    except Exception as e:
        sys.stderr.write(f"[ERROR] {cmd}: {e}\n")
        sys.exit(1)

def merge_ranges(ranges):
    """
    Merge overlapping or consecutive genomic ranges.
    
    Takes a list of (start, end) tuples representing genomic ranges and merges
    any ranges that overlap or are directly adjacent (consecutive positions).
    This is used to calculate accurate linear coverage for mapped reads.
    
    Parameters:
        ranges (list): List of (start, end) tuples representing genomic ranges
        
    Returns:
        list: List of merged (start, end) tuples with no overlaps
    
    Example:
        >>> merge_ranges([(1, 5), (4, 8), (10, 12)])
        [(1, 8), (10, 12)]
    """
    # Sort ranges by start position
    sorted_ranges = sorted(ranges, key=lambda x: x[0])
    merged = []
    for current in sorted_ranges:
        if not merged:
            merged.append(current)
        else:
            last = merged[-1]
            # Check if current range overlaps or is adjacent to the last range
            if current[0] <= last[1] + 1:
                # Merge the ranges
                merged[-1] = (last[0], max(last[1], current[1]))
            else:
                merged.append(current)
    return merged

def worker(filename, chunkStart, chunkSize, matchFactor):
    """
    Process a chunk of a SAM file to extract mapping information.
    
    This function is intended to be run in parallel to process different chunks
    of a SAM file. It parses lines within the specified chunk and builds a dictionary
    of reference sequences with their mapped regions, base counts, read counts, etc.
    
    Parameters:
        filename (str): Path to the SAM file to process
        chunkStart (int): Byte position in the file where to start reading
        chunkSize (int): Number of bytes to read from the start position
        matchFactor (float): Minimum fraction required for a valid match
        
    Returns:
        dict: Dictionary with reference sequences as keys and mapping statistics as values
    """
    # processing alignments in SAM format
    f = open( filename )
    f.seek(chunkStart)
    lines = f.read(chunkSize).splitlines()
    res={}

    for line in lines:
        k, r, n, rd, rs, rq, flag, cigr, pri_aln_flag, valid_flag = parse(line, matchFactor)
        # parsed values from SAM line
        # only k, r, n, pri_aln_flag, valid_flag are used
        # k: reference name
        # r: (start, end) of mapped region
        # n: number of mismatches
        # pri_aln_flag: whether this is a primary alignment
        # valid_flag: whether this alignment meets match criteria

        if pri_aln_flag and valid_flag:
            if k in res:
                res[k]['REGIONS'] = merge_ranges(res[k]['REGIONS']+[r])
                res[k]["MB"] += r[1] - r[0] + 1
                res[k]["MR"] += 1
                res[k]["NM"] += n
            else:
                res[k]={}
                res[k]["REGIONS"] = [r]
                res[k]["MB"] = r[1] - r[0] + 1
                res[k]["MR"] = 1
                res[k]["NM"] = n
    return res

def parse(line, matchFactor):
    """
    Parse a line from a SAM file and extract relevant mapping information.
    
    Parses alignment details from a SAM format line, including reference ID, 
    match position, mismatches, sequence quality, and flags. Determines if
    the alignment is a valid match based on matchFactor criteria.
    
    Parameters:
        line (str): A line from a SAM file
        matchFactor (float): Minimum fraction required for a valid match
        
    Returns:
        tuple: (
            ref (str): Reference identifier,
            (start, end) (tuple): Mapped region coordinates,
            mismatches (int): Number of mismatches,
            read_name (str): Name of the read,
            read_seq (str): Read sequence,
            read_qual (str): Read quality string, 
            flag (str): SAM flag,
            cigar (str): CIGAR string,
            primary_alignment_flag (bool): Whether this is a primary alignment,
            valid_flag (bool): Whether this alignment meets match criteria
        )
    
    Example:
        SAM format example:
        read1   0   ABC|1|100|GCF_12345|    11  0   5S10M3S *   0   0   GGGGGCCCCCCCCCCGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
        read2   16  ABC|1|100|GCF_12345|    11  0   3S10M5S *   0   0   GGGCCCCCCCCCCGGGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
    """
    temp = line.split('\t')
    name = temp[0]
    match_len    = search(r'(\d+)M', temp[5])
    mismatch_len = search(r'NM:i:(\d+)', line)
    start = int(temp[3])
    end   = start + int(match_len.group(1)) - 1

    ref = temp[2].rstrip('|')
    ref = ref[: -2 if ref.endswith(".0") else None ]

    (acc, rstart, rend, taxid) = ref.split('|')
    rlen = int(rend)-int(rstart)+1

    # check if this is a primary alignment 256=secondary, 2048=supplementary
    primary_alignment_flag=False if int(temp[1]) & 2304 else True

    # the alignment region should cover at least matchFactor(proportion) of the read or signature fragment
    valid_flag = True # default to True

    if matchFactor > 0:
        if (int(match_len.group(1)) >= rlen * matchFactor) or (int(match_len.group(1)) >= len(temp[9])*matchFactor):
            valid_flag = True
        else:
            valid_flag = False

    return ref, (start, end), int(mismatch_len.group(1)), name, temp[9], temp[10], temp[1], temp[5], primary_alignment_flag, valid_flag

def time_spend(start):
    """
    Calculate and format elapsed time since a given start time.
    
    Parameters:
        start (float): Starting time in seconds (as returned by time.time())
        
    Returns:
        str: Formatted time string in HH:MM:SS format
    """
    done = time.time()
    elapsed = done - start
    return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

def chunkify(fname, size=1*1024*1024):
    """
    Split a file into chunks for parallel processing.
    
    Divides a file into chunks of approximately the specified size, ensuring
    that all alignments for a single read are kept in the same chunk.
    This is critical for accurate processing of multi-mapped reads.
    
    Parameters:
        fname (str): Path to the file to be chunked
        size (int): Approximate chunk size in bytes (default: 1MB)
        
    Yields:
        tuple: (chunkStart, chunkSize) where:
            - chunkStart (int): Byte position to start reading
            - chunkSize (int): Number of bytes to read from that position
    """
    fileEnd = os.path.getsize(fname)
    with open(fname, "rb") as f:
        chunkEnd = f.tell()
        while True:
            chunkStart = chunkEnd
            f.seek(size, 1)
            f.readline()
            # put all alignments of a read in the same chunck
            line = f.readline().decode('ascii')
            tmp = line.split('\t')
            if chunkEnd <= fileEnd and line:
                last = f.tell()
                line = f.readline().decode('ascii')
                while line.startswith(tmp[0]):
                    last = f.tell()
                    line = f.readline().decode('ascii')
                f.seek(last)
            # current position
            chunkEnd = f.tell()
            yield chunkStart, chunkEnd - chunkStart
            if chunkEnd > fileEnd:
                break

def process_sam_file(sam_fn, numthreads, matchFactor):
    """
    Process a SAM file using parallel execution to extract mapping information.
    
    Divides the SAM file into chunks, processes each chunk in parallel using a thread pool,
    and then merges the results. Computes the linear coverage for each reference sequence.
    
    Parameters:
        sam_fn (str): Path to the SAM file
        numthreads (int): Number of parallel processes to use
        matchFactor (float): Minimum fraction required for a valid match
        
    Returns:
        tuple: (
            result (dict): Dictionary with references as keys and mapping statistics as values,
            mapped_reads (int): Total number of reads that mapped
        )
    """
    result = gt._autoVivification()
    mapped_reads = 0

    print_message( f"Parsing SAM files with {numthreads} subprocesses...", argvs.silent, begin_t, logfile )
    pool = Pool(processes=numthreads)
    jobs = []
    results = []

    for chunkStart,chunkSize in chunkify(sam_fn):
        jobs.append( pool.apply_async(worker, (sam_fn,chunkStart,chunkSize,matchFactor)) )

    #wait for all jobs to finish
    tol_jobs = len(jobs)
    cnt=0
    for job in jobs:
        results.append( job.get() )
        cnt+=1
        if argvs.debug: print_message( f"[DEBUG] Progress: {cnt}/{tol_jobs} ({cnt/tol_jobs*100:.1f}) chunks done.", argvs.silent, begin_t, logfile )

    #clean up
    pool.close()

    print_message( "Merging results...", argvs.silent, begin_t, logfile )
    for res in results:
        for k in res:
            if k in result:
                result[k]['REGIONS'] = merge_ranges(result[k]['REGIONS']+res[k]['REGIONS'])
                result[k]["MB"] += res[k]["MB"]
                result[k]["MR"] += res[k]["MR"]
                result[k]["NM"] += res[k]["NM"]
            else:
                result[k]={}
                result[k].update(res[k])

    # convert mapped regions to linear length
    refs = result.keys()
    for k in list(refs):
        if not result[k]["MR"]:
            del result[k]
        else:
            result[k]["LL"] = sum(end - start + 1 for start, end in result[k]['REGIONS'])
            del result[k]['REGIONS']
            mapped_reads += result[k]["MR"]

    return result, mapped_reads

def is_descendant(taxid, taxid_ant):
    """
    Check if one taxid is a descendant of another in the taxonomy tree.
    
    Parameters:
        taxid (str): The taxid to check
        taxid_ant (str): The potential ancestor taxid
        
    Returns:
        bool: True if taxid is a descendant of taxid_ant, False otherwise
    """
    fullLineage = gt.taxid2fullLineage( taxid )
    if "|%s|" % taxid_ant in fullLineage:
        return True
    else:
        return False

def extract_read_from_sam(sam_fn, o, taxid, numthreads, matchFactor):
    """
    Extract reads from a SAM file that map to a specific taxid or its descendants.
    
    Processes the SAM file in parallel chunks and writes matching reads in FASTQ format
    to the provided output file.
    
    Parameters:
        sam_fn (str): Path to the SAM file
        o (file): Output file handle to write extracted reads
        taxid (str): Taxid to extract reads for (including descendants)
        numthreads (int): Number of parallel processes to use
        matchFactor (float): Minimum fraction required for a valid match
        
    Returns:
        None: Results are written directly to the output file
    """
    pool = Pool(processes=numthreads)
    jobs = []

    for chunkStart,chunkSize in chunkify(sam_fn):
        jobs.append( pool.apply_async(ReadExtractWorker, (sam_fn,chunkStart,chunkSize,taxid,matchFactor)) )

    #wait for all jobs to finish
    for job in jobs:
        outread = job.get()
        o.write(outread)
        o.flush()

    #clean up
    pool.close()

def ReadExtractWorker(filename, chunkStart, chunkSize, taxids, matchFactor):
    """Extract reads matching any of the specified taxids"""
    # Convert input to list of taxids
    taxid_list = parse_taxids(taxids)
    if not taxid_list:
        return ""
        
    readstr = ""
    with open(filename) as f:
        f.seek(chunkStart)
        lines = f.read(chunkSize).splitlines()
        
        for line in lines:
            ref, region, nm, rname, rseq, rq, flag, cigr, pri_aln_flag, valid_flag = parse(line, matchFactor)
            if not (pri_aln_flag and valid_flag): 
                continue

            acc, start, stop, t = ref.split('|')
            
            # Check if current taxid matches any requested taxids
            for taxid in taxid_list:
                if is_descendant(t, taxid):
                    if int(flag) & 16:
                        g = findall(r'\d+\w', cigr)
                        cigr = "".join(list(reversed(g)))
                        rseq = seqReverseComplement(rseq)
                        rq = rq[::-1]
                    readstr += "@%s %s:%s..%s %s\n%s\n+\n%s\n" % (rname, ref, region[0], region[1], cigr, rseq, rq)
                    break
                    
    return readstr

def seqReverseComplement(seq):
    """
    Generate the reverse complement of a DNA sequence.
    
    Creates a mapping dictionary for complementary bases and applies it to the
    reversed sequence. Handles both uppercase and lowercase nucleotides.
    
    Parameters:
        seq (str): DNA sequence string
        
    Returns:
        str: Reverse complemented DNA sequence
    """
    seq1 = 'ACGTURYSWKMBDHVNTGCAAYRSWMKVHDBNacgturyswkmbdhvntgcaayrswmkvhdbn'
    seq_dict = { seq1[i]:seq1[i+16] for i in range(64) if i < 16 or 32<=i<48 }
    return "".join([seq_dict[base] for base in reversed(seq)])

def group_refs_to_strains(r):
    """
    Group reference mapping results by strains and calculate strain-level statistics.
    
    Converts the mapping results dictionary to a pandas DataFrame and groups by
    taxonomic identifier. Calculates various statistics including total mapped bases,
    read counts, coverage, and depth of coverage.
    
    Parameters:
        r (dict): Dictionary with reference sequences as keys and mapping statistics
                 as values (output from process_sam_file)
        
    Returns:
        pandas.DataFrame: DataFrame with strain-level statistics
    """
    # covert mapping info to df
    r_df = pd.DataFrame.from_dict(r, orient='index').reset_index()
    r_df.rename(columns={"index": "RNAME"}, inplace=True)
    # retrieve sig fragment info
    r_df['RNAME'] = r_df['RNAME'].str.rstrip('|')
    r_df[['ACC','RSTART','REND','TAXID']] = r_df['RNAME'].str.split('|', expand=True)
    r_df['RSTART'] = r_df['RSTART'].astype(int)
    r_df['REND'] = r_df['REND'].astype(int)
    r_df['RLEN'] = r_df['REND']-r_df['RSTART']+1

    # group by strain
    str_df = r_df.groupby(['TAXID']).agg({
        'MB':'sum', # of mapped bases
        'MR':'sum', # of mapped reads
        'NM':'sum', # of mismatches
        'LL':'sum', # linear length
        'RLEN':'sum' # length of this signature fragments (mapped)
    }).reset_index()
    # total length of signatures
    str_df['TS'] = str_df['TAXID'].apply(lambda x: db_stats[x])
    str_df['bDOC'] = str_df['MB']/str_df['TS'] # bDOC: best Depth of Coverage of a strain
    str_df['bLC'] = str_df['LL']/str_df['TS'] # bLC:  best linear coverage of a strain
    str_df['RD'] = str_df['MB']/str_df['TS'] # roll-up DoC

    # rename columns
    str_df.rename(columns={
        "MB":   "TOTAL_BP_MAPPED",
        "MR":   "READ_COUNT",
        "NM":   "TOTAL_BP_MISMATCH",
        "LL":   "LINEAR_LEN",
        "RLEN": "MAPPED_SIG_LENGTH",
        "TS":   "TOL_SIG_LENGTH",
        "RD":   "ROLLUP_DOC",
        "bDOC": "BEST_DOC",
        "bLC":  "BEST_LINEAR_COV"
    }, inplace=True)

    str_df['ZSCORE'] = str_df.apply(lambda x: pile_lvl_zscore(x.TOTAL_BP_MAPPED, x.TOL_SIG_LENGTH, x.LINEAR_LEN), axis=1)

    return str_df

def aggregate_taxonomy(r, abu_col, tg_rank, mc, mr, ml, mz):
    """
    Aggregate strain-level results to higher taxonomic ranks.
    
    Starting from strain-level mapping data, this function rolls up statistics to
    higher taxonomic ranks (species, genus, family, etc.). It applies the specified
    cutoff criteria to filter results and marks entries that fall below these thresholds.
    
    Parameters:
        r (dict): Dictionary with reference sequences as keys and mapping stats as values
        abu_col (str): Column name to use for abundance calculations
        tg_rank (str): Target taxonomic rank
        mc (float): Minimum linear coverage threshold
        mr (int): Minimum read count threshold
        ml (int): Minimum linear length threshold
        mz (float): Maximum Z-score threshold (0 to disable)
        
    Returns:
        pandas.DataFrame: DataFrame with rolled-up taxonomy at all ranks
    """
    major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}

    # roll up references to strains
    str_df = group_refs_to_strains(r)
    # produce columns for the final report at each ranks
    rep_df = pd.DataFrame()

    # qualified strain
    qualified_idx = (str_df['LINEAR_LEN']/str_df['TOL_SIG_LENGTH'] >= mc) & \
                    (str_df['READ_COUNT'] >= mr) & \
                    (str_df['LINEAR_LEN'] >= ml)
    
    if mz > 0:
        qualified_idx &= (str_df['ZSCORE'] <= mz)

    for rank in sorted(major_ranks, key=major_ranks.__getitem__):
        try:
            str_df['LVL_NAME'] = str_df['TAXID'].apply(lambda x: gt.taxid2lineageDICT(x, True, True)[rank]['name'])
            str_df['LVL_TAXID'] = str_df['TAXID'].apply(lambda x: gt.taxid2lineageDICT(x, True, True)[rank]['taxid'])
            str_df['LEVEL'] = rank
        except Exception as e:
            logging.error(f"Error processing rank {rank}: {e}. Please verify that your taxonomy file matches the expected database.")
            sys.exit(1)

        # rollup strains that make cutoffs
        lvl_df = pd.DataFrame()
        if rank == 'strain':
            lvl_df = str_df
            lvl_df['LVL_TAXID'] = str_df['TAXID']
        else:
            lvl_df = str_df[qualified_idx].groupby(['LVL_NAME']).agg({
                'LEVEL':'first',
                'LVL_TAXID':'first',
                'TOTAL_BP_MAPPED': 'sum', 'READ_COUNT': 'sum', 'TOTAL_BP_MISMATCH': 'sum',
                'LINEAR_LEN': 'sum', 'MAPPED_SIG_LENGTH': 'sum', 'TOL_SIG_LENGTH': 'sum',
                'ROLLUP_DOC': 'sum', 'BEST_DOC': 'max', 'BEST_LINEAR_COV': 'max', 'ZSCORE': 'min',
            }).reset_index().copy()

        lvl_df['ABUNDANCE'] = lvl_df[abu_col]
        tol_abu = lvl_df[abu_col].sum()
        lvl_df['REL_ABUNDANCE'] = lvl_df[abu_col]/tol_abu

        # add NOTE if ranks is higher than target rank
        lvl_df['NOTE'] = ""
        if major_ranks[rank] > major_ranks[tg_rank]:
            lvl_df['NOTE'] = f"Not shown ({rank}-result biased); "

        # concart ranks-dataframe to the report-dataframe
        rep_df = pd.concat([rep_df, lvl_df.sort_values('ABUNDANCE', ascending=False)], sort=False)

    rep_df["LINEAR_COV"] = rep_df["LINEAR_LEN"]/rep_df["TOL_SIG_LENGTH"]
    rep_df["LINEAR_COV_MAPPED_SIG"] = rep_df["LINEAR_LEN"]/rep_df["MAPPED_SIG_LENGTH"]
    rep_df["DOC"] = rep_df["TOTAL_BP_MAPPED"]/rep_df["TOL_SIG_LENGTH"]
    rep_df["LINEAR_DOC"] = rep_df["TOTAL_BP_MAPPED"]/rep_df["LINEAR_LEN"]

    # add filtered reason to NOTE column
    filtered = (rep_df['LINEAR_COV'] < mc)
    rep_df.loc[filtered, 'NOTE'] += "Filtered out (minCov > " + rep_df.loc[filtered, 'LINEAR_COV'].astype(str) + "); "
    filtered = (rep_df['READ_COUNT'] < mr)
    rep_df.loc[filtered, 'NOTE'] += "Filtered out (minReads > " + rep_df.loc[filtered, 'READ_COUNT'].astype(str) + "); "
    filtered = (rep_df['LINEAR_LEN'] < ml)
    rep_df.loc[filtered, 'NOTE'] += "Filtered out (minLen > " + rep_df.loc[filtered, 'LINEAR_LEN'].astype(str) + "); "

    if mz > 0:
        filtered = (rep_df['ZSCORE'] > mz)
        rep_df.loc[filtered, 'NOTE'] += "Filtered out (maxZscore < " + rep_df.loc[filtered, 'ZSCORE'].astype(str) + "); "

    rep_df.drop(columns=['TAXID'], inplace=True)
    rep_df.rename(columns={"LVL_NAME": "NAME", "LVL_TAXID": "TAXID"}, inplace=True)

    logging.debug(f'rep_df: {rep_df}')

    return rep_df

def pile_lvl_zscore(tol_bp, tol_sig_len, linear_len):
    """
    Calculate Z-score for the depth of coverage of mapped regions.
    
    This determines how unusual the coverage depth is compared to expected depth
    based on a statistical model. Higher Z-scores may indicate biased mapping.
    
    Parameters:
        tol_bp (int): Total number of mapped bases
        tol_sig_len (int): Total length of the signature
        linear_len (int): Linear length (de-duplicated) covered by mappings
        
    Returns:
        float: Z-score for the depth distribution (or 0 if calculation fails)
    """
    try:
        avg_doc = tol_bp/tol_sig_len
        lin_doc = tol_bp/linear_len
        v = (linear_len*(lin_doc-avg_doc)**2 + (tol_sig_len-linear_len)*(avg_doc)**2)/tol_sig_len
        sd = math.sqrt(v)
        if sd == 0.0:
            return 0
        else:
            return (lin_doc-avg_doc)/sd
    except:
        return 0
    
def estimate_ani(tol_bp, tol_mismatch):
    """
    """
    try:
        return (1 - tol_mismatch/tol_bp)
    except:
        return 0

def generaete_taxonomy_file(rep_df, o, fullreport_o, fmt="tsv"):
    """
    Generate taxonomy result files in TSV or CSV format.
    
    Creates two files: a summary file with only qualified results, and
    a full report with all results including filtered entries.
    
    Parameters:
        rep_df (pandas.DataFrame): Taxonomy results DataFrame
        o (file): Output file handle for the summary results
        fullreport_o (str): Path for the full report file
        fmt (str): Output format, either 'tsv' or 'csv'
        
    Returns:
        bool: True if successful
    """
    # Fields for full mode
    cols = ['LEVEL', 'NAME', 'TAXID', 'READ_COUNT', 'TOTAL_BP_MAPPED',
            'TOTAL_BP_MISMATCH', 'LINEAR_LEN', 'LINEAR_DOC', 'ROLLUP_DOC', 'REL_ABUNDANCE',
            'LINEAR_COV', 'LINEAR_COV_MAPPED_SIG', 'BEST_LINEAR_COV', 'MAPPED_SIG_LENGTH', 'TOL_SIG_LENGTH',
            'ABUNDANCE', 'ZSCORE', 'NOTE']

    qualified_idx = (rep_df['NOTE']=="")
    qualified_df = rep_df.loc[qualified_idx, cols[:10]]

    sep = ',' if fmt=='csv' else '\t'
    # save full report
    rep_df[cols].to_csv(fullreport_o, index=False, sep=sep, float_format='%.6f')
    # save summary
    qualified_df.to_csv(o, index=False, sep=sep, float_format='%.6f')

    return True

def generaete_biom_file(res_df, o, tg_rank, sampleid):
    """
    Generate a BIOM format file from taxonomy results.
    
    Creates a BIOM (Biological Observation Matrix) formatted file for
    compatibility with downstream microbiome analysis tools.
    
    Parameters:
        res_df (pandas.DataFrame): Taxonomy results DataFrame
        o (file): Output file handle
        tg_rank (str): Target taxonomic rank to include in the output
        sampleid (str): Sample identifier
        
    Returns:
        bool: True if successful
        
    Raises:
        SystemExit: If the biom library version is incompatible
    """
    import numpy as np
    import biom
    from biom.table import Table
    if biom.__version__ < '2.1.7':
        sys.exit("[ERROR] Biom library requires v2.1.7 or above.\n")

    target_df = pd.DataFrame()
    target_idx = (res_df['LEVEL']==tg_rank)
    target_df = res_df.loc[target_idx, ['ABUNDANCE','TAXID']]
    target_df['LINEAGE'] = target_df['TAXID'].apply(lambda x: gt.taxid2lineage(x, True, True)).str.split('|')

    sample_ids = [sampleid]
    data = np.array(target_df['ABUNDANCE']).reshape(len(target_df), 1)
    observ_ids = target_df['TAXID']
    observ_metadata = [{'taxonomy': x} for x in target_df['LINEAGE'].tolist()]
    biom_table = Table(data, observ_ids, sample_ids, observ_metadata, table_id='GOTTCHA2')
    biom_table.to_json('GOTTCHA2', direct_io=o)

    return True

def generaete_lineage_file(target_df, o, tg_rank):
    """
    Generate a lineage file showing taxonomic paths with abundances.
    
    Creates a tab-delimited file with abundance values followed by
    the complete taxonomic lineage for each taxon.
    
    Parameters:
        target_df (pandas.DataFrame): DataFrame containing abundance and taxids
        o (str): Output file path
        tg_rank (str): Target taxonomic rank
        
    Returns:
        bool: True if successful
    """
    lineage_df = target_df['TAXID'].apply(lambda x: gt.taxid2lineage(x, True, True)).str.split('|', expand=True)
    result = pd.concat([target_df['ABUNDANCE'], lineage_df], axis=1, sort=False)
    result.to_csv(o, index=False, header=False, sep='\t', float_format='%.4f')

    return True

def readMapping(reads, db, threads, mm_penalty, presetx, samfile, logfile):
    """
    Map reads to the reference database using minimap2.
    
    Builds and executes a command to run minimap2 for read mapping, with parameters
    adjusted based on input settings. Filters the SAM output to keep only relevant
    alignments.
    
    Parameters:
        reads (list): List of input read file objects
        db (str): Path to the minimap2 database (without .mmi extension)
        threads (int): Number of threads to use
        mm_penalty (int): Mismatch penalty for alignment
        presetx (str): Minimap2 preset mode ('sr', 'map-pb', or 'map-ont')
        samfile (str): Output SAM file path
        logfile (str): Log file path
        nanopore (bool): Whether to use Nanopore-specific settings
        
    Returns:
        tuple: (
            exitcode (int): Exit code from the mapping process,
            cmd (str): Command that was executed,
            errs (str): Error output from the command
        )
    """
    input_file = " ".join([x.name for x in reads])

    # Minimap2 options for short reads: the options here is essentailly the -x 'sr' equivalent with some modifications on scoring
    sr_opts = f"--sr --frag=yes -b0 -r100 -f1000,5000 -n2 -m20 -s40 -g100 -2K50m -k24 -w12 -A1 -B{mm_penalty} -O30 -E30 -a -N20 --secondary=no --sam-hit-only"
    
    if presetx != 'sr':
        sr_opts = f"-x {presetx} -N20 --secondary=no --sam-hit-only -a"

    bash_cmd   = f"set -o pipefail; set -x;"
    mm2_cmd    = f"minimap2 {sr_opts} -t{threads} {db}.mmi {input_file}"
    filter_cmd = f"grep -v '^@'"
    cmd        = f"{bash_cmd} {mm2_cmd} 2>> {logfile} | {filter_cmd} > {samfile}"

    logging.info(f"Readmapping command: {mm2_cmd}")

    proc = subprocess.Popen(cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = proc.communicate()
    exitcode = proc.poll()

    return exitcode, mm2_cmd, errs

def remove_multiple_hits(samfile):
    """
    Removing multiple hits from the SAM file by keeping only the best alignment for each read.

    Parameters:
        samfile (str): Path to the SAM file

    Returns:
        str: Path to the temporary SAM file with only the best alignments
    """
    logging.info(f'Loading the sam file...')

    df = pd.read_csv(samfile,
                sep='\t',
                header=None,
                usecols=[0, 1, 13],
                names=['QNAME', 'FLAG', 'AS'],
                converters={
                    'AS': lambda x: x.replace('AS:i:', '')
                },
                dtype={'QNAME': 'str', 'FLAG': 'uint16'}
    )

    logging.info(f'Total alignments in SAM file: {len(df)}')

    df[['AS']] = df[['AS']].astype('int16') 

    logging.info(f'Filtering non-primary hits...')

    # for each row, if the flag bitwise AND with 256 (not primary alignment) or 2048 (supplementary), then remove them from the df
    df = df[~(df['FLAG'] & (256|2048)).astype(bool)]

    logging.info(f'After removing non-parmary hits: {len(df)}')

    logging.info(f'Identifying top score hits...')

    # if FLAG bitwise AND with 128 (second in pair), append '/2' to the QNAME
    idx = (df['FLAG'] & 128).astype(bool)
    df.loc[idx, 'QNAME'] = df.loc[idx, 'QNAME'] + '/2'

    # get the index with the best alignment score for each read
    idxmax = df.groupby('QNAME')['AS'].idxmax()

    logging.info(f'Total top score hits: {len(idxmax)}')

    # Create a set of indices for faster lookup
    idxmax_set = set(idxmax.values)

    del idxmax

    logging.info(f'Writing top score hits...')

    # Use a buffered approach for better I/O performance
    buffer_size = 100000  # Number of lines to process at once

    with open(f'{samfile}.temp', 'w') as fout, open(samfile, 'r') as fin:
        for idx, line in enumerate(fin):
            if not idx%100000:
                logging.debug(f'Processed {idx} lines...')
            
            if idx in idxmax_set:
                fout.write(line)

    logging.info(f'Done writing hits.')

    return f'{samfile}.temp'


def loadDatabaseStats(db_stats_file):
    """
    Load database signature statistics from a stats file.
    
    Reads a tab-delimited stats file containing information about
    taxonomic signatures and their lengths.
    
    Parameters:
        db_stats_file (str): Path to the database stats file
        
    Returns:
        dict: Dictionary with taxids as keys and signature lengths as values
        
    Note:
        The input stats file is an 8-column tab-delimited file with:
        1. Rank
        2. Name
        3. Taxid
        4. Superkingdom
        5. NumOfSeq
        6. Max
        7. Min
        8. TotalLength
    """
    db_stats = {}

    with open(db_stats_file) as f:
        for line in f:
            fields = line.split("\t")
            major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7, "strain":8}
            if fields[0] in major_ranks:
                db_stats[fields[2]] = int(fields[7])
            else:
                continue

    return db_stats

def print_message(msg, silent, start, logfile, errorout=0):
    """
    Print and log a timestamped message.
    
    Writes a message to the log file and optionally to stderr. Can also
    terminate the program with an error message.
    
    Parameters:
        msg (str): Message to print
        silent (bool): If True, suppress output to stderr
        start (float): Start time for timestamp calculation
        logfile (str): Path to the log file
        errorout (int): If non-zero, exit with error after printing
        
    Returns:
        None
        
    Raises:
        SystemExit: If errorout is non-zero
    """
    message = "[%s] %s\n" % (time_spend(start), msg)

    with open( logfile, "a" ) as f:
        f.write( message )
        f.close()

    if errorout:
        sys.exit( message )
    elif not silent:
        sys.stderr.write( message )

def parse_taxids(taxid_arg):
    """Parse taxids from command line arg or file"""
    if not taxid_arg:
        return []
    if taxid_arg.startswith('@'):
        # Read taxids from file
        filename = taxid_arg[1:]  # Remove @ prefix
        try:
            with open(filename) as f:
                return [x.strip() for x in f.readlines() if x.strip()]
        except IOError as e:
            sys.stderr.write(f"Error reading taxid file {filename}: {e}\n")
            sys.exit(1)
    else:
        # Parse comma-separated list
        return [x.strip() for x in taxid_arg.split(',')]

def main(args):
    """
    Main execution function for GOTTCHA2.
    """
    global argvs
    global logfile
    global begin_t
    global db_stats
    argvs = parse_params( __version__, args )
    begin_t  = time.time()
    sam_fp   = argvs.sam[0] if argvs.sam else ""
    samfile  = "%s/%s.gottcha_%s.sam" % ( argvs.outdir, argvs.prefix, argvs.dbLevel ) if not argvs.sam else sam_fp.name
    logfile  = "%s/%s.gottcha_%s.log" % ( argvs.outdir, argvs.prefix, argvs.dbLevel )

    logging_level = logging.WARNING

    if argvs.debug:
        logging_level = logging.DEBUG
    elif argvs.silent:
        logging_level = logging.FATAL
    elif argvs.verbose:
        logging_level = logging.INFO

    logging.basicConfig(
        level=logging_level,
        format='%(asctime)s [%(levelname)s] %(module)s: %(message)s',
        datefmt='%Y-%m-%d %H:%M',
    )

    #dependency check
    if sys.version_info < (3,6):
        sys.exit("[ERROR] Python 3.6 or above is required.")
    
    dependency_check("minimap2")
    # dependency_check("gawk")

    #prepare output object
    argvs.relAbu = argvs.relAbu.upper()
    outfile_full = "%s/%s.full.tsv" % (argvs.outdir, argvs.prefix)
    outfile_lineage = "%s/%s.lineage.tsv" % (argvs.outdir, argvs.prefix)

    # remove previous log file if exists
    if os.path.isfile(logfile):
        os.remove(logfile)

    out_fp = sys.stdout
    outfile = "STDOUT"

    if not argvs.stdout:
        #create output directory if not exists
        if not os.path.exists(argvs.outdir):
            os.makedirs(argvs.outdir)
        ext = "tsv"
        outfile = "%s/%s.tsv" % (argvs.outdir, argvs.prefix)
        if argvs.extract:
            outfile = "%s/%s.extract%s.fastq" % (argvs.outdir, argvs.prefix, argvs.extract)
        elif argvs.format == "csv":
            outfile = "%s/%s.csv" % (argvs.outdir, argvs.prefix)
        elif argvs.format == "biom":
            outfile = "%s/%s.biom" % (argvs.outdir, argvs.prefix)

        out_fp = open(outfile, 'w')

    print_message( f"Starting GOTTCHA (v{__version__})", argvs.silent, begin_t, logfile )
    print_message( f"Arguments and dependencies checked:", argvs.silent, begin_t, logfile )
    if argvs.input:
        print_message( f"    Input reads      : {[x.name for x in argvs.input]}",     argvs.silent, begin_t, logfile )
    print_message( f"    Input SAM file   : {samfile}",           argvs.silent, begin_t, logfile )
    print_message( f"    Database         : {argvs.database}",    argvs.silent, begin_t, logfile )
    print_message( f"    Database level   : {argvs.dbLevel}",     argvs.silent, begin_t, logfile )
    print_message( f"    Mismatch penalty : {argvs.mismatch}",    argvs.silent, begin_t, logfile )
    print_message( f"    Abundance        : {argvs.relAbu}",      argvs.silent, begin_t, logfile )
    print_message( f"    Output path      : {argvs.outdir}",      argvs.silent, begin_t, logfile )
    print_message( f"    Prefix           : {argvs.prefix}",      argvs.silent, begin_t, logfile )
    print_message( f"    Extract taxid    : {argvs.extract}",     argvs.silent, begin_t, logfile )
    print_message( f"    Threads          : {argvs.threads}",     argvs.silent, begin_t, logfile )
    print_message( f"    Minimal L_DOC    : {argvs.minCov}",      argvs.silent, begin_t, logfile )
    print_message( f"    Minimal L_LEN    : {argvs.minLen}",      argvs.silent, begin_t, logfile )
    print_message( f"    Minimal reads    : {argvs.minReads}",    argvs.silent, begin_t, logfile )
    print_message( f"    Minimal mFactor  : {argvs.matchFactor}", argvs.silent, begin_t, logfile )
    print_message( f"    Maximal zScore   : {argvs.maxZscore}",   argvs.silent, begin_t, logfile )

    #load taxonomy
    print_message( "Loading taxonomy information...", argvs.silent, begin_t, logfile )
    custom_taxa_tsv = None
    if os.path.isfile( argvs.database + ".tax.tsv" ):
        custom_taxa_tsv = argvs.database + ".tax.tsv"
    
    gt.loadTaxonomy( argvs.taxInfo, custom_taxa_tsv )
    print_message( "Done.", argvs.silent, begin_t, logfile )

    #load database stats
    print_message( "Loading database stats...", argvs.silent, begin_t, logfile )
    if os.path.isfile( argvs.database + ".stats" ):
        db_stats = loadDatabaseStats(argvs.database+".stats")
    else:
        sys.exit( "[%s] ERROR: %s not found.\n" % (time_spend(begin_t), argvs.database+".stats") )
    print_message( "Done.", argvs.silent, begin_t, logfile )

    #main process
    if argvs.input:
        print_message( "Running read-mapping...", argvs.silent, begin_t, logfile )
        exitcode, cmd, msg = readMapping( argvs.input, argvs.database, argvs.threads, argvs.mismatch, argvs.presetx, samfile, logfile)
        gc.collect()
        print_message( f"Logfile saved to {logfile}.", argvs.silent, begin_t, logfile )
        print_message( f"COMMAND: {cmd}", argvs.silent, begin_t, logfile )

        if exitcode != 0:
            sys.exit( "[%s] ERROR: error occurred while running read mapping (exit: %s, message: %s).\n" % (time_spend(begin_t), exitcode, msg) )
        else:
            print_message( f"Done mapping reads to {argvs.dbLevel} signature database.", argvs.silent, begin_t, logfile )
            print_message( f"Mapped SAM file saved to {samfile}.", argvs.silent, begin_t, logfile )
            sam_fp = open( samfile, "r" )

    # remove multiple hits
    if not argvs.skipRemoveMultiple:
        # remove multiple hits from the SAM file
        print_message( "Removing multiple hits from SAM file...", argvs.silent, begin_t, logfile )
        samfile_temp = remove_multiple_hits(samfile)
        os.rename(samfile_temp, samfile)
        gc.collect()

    if argvs.extract:
        print_message( f"Extracting reads mapped to taxid: {argvs.extract}...", argvs.silent, begin_t, logfile )
        extract_read_from_sam( os.path.abspath(samfile), out_fp, argvs.extract, argvs.threads, argvs.matchFactor )
        print_message( f"Done extracting reads to {outfile}.", argvs.silent, begin_t, logfile )
    else:
        print_message( "Loading SAM file...", argvs.silent, begin_t, logfile )
        (res, mapped_r_cnt) = process_sam_file( os.path.abspath(samfile), argvs.threads, argvs.matchFactor)
        print_message( f"Done processing SAM file. {mapped_r_cnt} qualified mapped reads loaded.", argvs.silent, begin_t, logfile )
        gc.collect()

        if mapped_r_cnt:
            res_df = aggregate_taxonomy(res, argvs.relAbu, argvs.dbLevel , argvs.minCov, argvs.minReads, argvs.minLen, argvs.maxZscore)
            print_message( "Done taxonomy rolling up.", argvs.silent, begin_t, logfile )

            if not len(res_df):
                print_message( "No qualified taxonomy profiled.", argvs.silent, begin_t, logfile )
            else:
                # generate output files
                if argvs.format == "biom":
                    generaete_biom_file(res_df, out_fp, argvs.dbLevel, argvs.prefix)
                else:
                    generaete_taxonomy_file(res_df, out_fp, outfile_full, argvs.format)
                # generate lineage file
                target_idx = (res_df['LEVEL']==argvs.dbLevel)
                target_df = res_df.loc[target_idx, ['ABUNDANCE','TAXID']]
                tax_num = len(target_df)
                if tax_num:
                     generaete_lineage_file(target_df, outfile_lineage, argvs.dbLevel)

                print_message( f"{tax_num} qualified {argvs.dbLevel} profiled; Results saved to {outfile}.", argvs.silent, begin_t, logfile )
        else:
            print_message( "GOTTCHA2 stopped.", argvs.silent, begin_t, logfile )

if __name__ == '__main__':
    main(sys.argv[1:])