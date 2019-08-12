#!/usr/bin/env python3

__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__credits__   = ["Po-E Li", "Jason Gans", "Tracey Freites", "Patrick Chain"]
__version__   = "2.1.5 BETA"
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
import taxonomy as gt
import pandas as pd
import gc
from re import search,findall
from re import compile as recompile
from multiprocessing import Pool
from itertools import chain
import math

def parse_params( ver ):
    p = ap.ArgumentParser( prog='gottcha2.py', description="""Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an
            annotation-independent and signature-based metagenomic taxonomic profiling tool
            that has significantly smaller FDR than other profiling tools. This program
            is a wrapper to map input reads to pre-computed signature databases using minimap2
            and/or to profile mapped reads in SAM format. (VERSION: %s)""" % ver)

    eg = p.add_mutually_exclusive_group( required=True )

    eg.add_argument( '-i','--input', metavar='[FASTQ]', nargs='+', type=str,
                    help="Input one or multiple FASTQ/FASTA file(s). Use space to separate multiple input files.")

    eg.add_argument( '-s','--sam', metavar='[SAMFILE]', nargs=1, type=ap.FileType('r'),
                    help="Specify the input SAM file. Use '-' for standard input.")

    p.add_argument( '-d','--database', metavar='[MINIMAP2_INDEX]', type=str, default=None,
                    help="The path of signature database. The database can be in FASTA format or minimap2 index (5 files).")

    p.add_argument( '-l','--dbLevel', metavar='[LEVEL]', type=str, default='',
                    choices=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'],
                    help="""Specify the taxonomic level of the input database. You can choose one rank from "superkingdom", "phylum", "class", "order", "family", "genus", "species" and "strain". The value will be auto-detected if the input database ended with levels (e.g. GOTTCHA_db.species).""")

    p.add_argument( '-ti','--taxInfo', metavar='[FILE]', type=str, default='',
                    help="""Specify the path of taxonomy information file (taxonomy.tsv). GOTTCHA2 will try to locate this file when user doesn't specify a path. If '--database' option is used, the program will try to find this file in the directory of specified database. If not, the 'database' directory under the location of gottcha.py will be used as default.""")

    p.add_argument( '-np','--nanopore', action="store_true",
                    help="Adjust options for Nanopore reads. The 'mismatch' option will be ignored. [-xm map-ont -mr 1]")
                    
    p.add_argument( '-pm','--mismatch', metavar='<INT>', type=int, default=10,
                    help="Mismatch penalty for the aligner. [default: 10]")

    p.add_argument( '-e','--extract', metavar='[TAXID]', type=str, default=None,
                    help="""Extract reads mapping to a specific TAXID. [default: None]""" )

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
                    help="Minimum linear coverage to be considered valid in abundance calculation [default: 0.005]")

    p.add_argument( '-mr','--minReads', metavar='<INT>', type=int, default=3,
                    help="Minimum number of reads to be considered valid in abundance calculation [default: 3]")

    p.add_argument( '-ml','--minLen', metavar='<INT>', type=int, default=60,
                    help="Minimum unique length to be considered valid in abundance calculation [default: 60]")

    p.add_argument( '-mz','--maxZscore', metavar='<FLOAT>', type=float, default=10,
                    help="Maximum estimated zscore of depths of mapped region [default: 10]")

    p.add_argument( '-nc','--noCutoff', action="store_true",
                    help="Remove all cutoffs. This option is equivalent to use [-mc 0 -mr 0 -ml 0].")

    p.add_argument( '-c','--stdout', action="store_true",
                    help="Write on standard output.")

    eg.add_argument( '-v','--version', action="store_true",
                    help="Print version number.")

    p.add_argument( '--silent', action="store_true",
                    help="Disable all messages.")

    p.add_argument( '--debug', action="store_true",
                    help="Debug mode. Provide verbose running messages and keep all temporary files.")

    args_parsed = p.parse_args()

    """
    Checking options
    """
    if args_parsed.version:
        print( ver )
        os._exit(0)

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
            args_parsed.taxInfo = db_dir.group(1)
        else:
            bin_dir = os.path.dirname(os.path.realpath(__file__))
            args_parsed.taxInfo = bin_dir + "/database"

    if not args_parsed.prefix:
        if args_parsed.input:
            name = search(r'([^\/\.]+)\..*$', args_parsed.input[0] )
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
            args_parsed.dbLevel = name.group(1)
        else:
            p.error( '--dbLevel is missing and cannot be auto-detected.' )

    if args_parsed.noCutoff:
        args_parsed.minCov = 0
        args_parsed.minReads = 0
        args_parsed.minLen = 0
        args_parsed.maxZscore = 99

    if args_parsed.nanopore:
        args_parsed.presetx = 'map-ont'
        args_parsed.minReads = 1

    return args_parsed

def dependency_check(cmd):
    try:
        subprocess.check_call([cmd, "--help"], stdout=subprocess.DEVNULL)
    except Exception as e:
        sys.stderr.write(f"[ERROR] {cmd}: {e}\n")
        sys.exit(1)

def worker(filename, chunkStart, chunkSize):
    """Make a dict out of the parsed, supplied lines"""
    # processing alignments in SAM format
    f = open( filename )
    f.seek(chunkStart)
    lines = f.read(chunkSize).splitlines()
    res={}

    for line in lines:
        k, r, m, n, rd, rs, rq, flag, cigr, pri_aln_flag, valid_flag = parse(line)
        if valid_flag:
            if k in res:
                res[k]["ML"] = res[k]["ML"] | m
                if pri_aln_flag and valid_flag:
                    res[k]["MB"] += r[1] - r[0] + 1
                    res[k]["MR"] += 1
                    res[k]["NM"] += n
            else:
                res[k]={}
                res[k]["ML"] = m
                if pri_aln_flag and valid_flag:
                    res[k]["MB"] = r[1] - r[0] + 1
                    res[k]["MR"] = 1
                    res[k]["NM"] = n
                else:
                    res[k]["MB"] = 0
                    res[k]["MR"] = 0
                    res[k]["NM"] = 0
    return res

def parse(line):
    """
    Parse SAM format
    read1   0   test    11  0   5S10M3S *   0   0   GGGGGCCCCCCCCCCGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
    read2   16  test    11  0   3S10M5S *   0   0   GGGCCCCCCCCCCGGGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
    """
    temp = line.split('\t')
    name = temp[0]
    match_len    = search(r'(\d+)M', temp[5])
    mismatch_len = search(r'NM:i:(\d+)', temp[11])
    start = int(temp[3])
    end   = start + int(match_len.group(1)) - 1

    ref = temp[2].rstrip('|')
    ref = ref[: -2 if ref.endswith(".0") else None ]

    (acc, rstart, rend, taxid) = ref.split('|')
    rlen = int(rend)-int(rstart)+1
    mask = int( "%s%s"%("1"*(end-start+1), "0"*(rlen-end)), 2)

    primary_alignment_flag=False if int(temp[1]) & 256 else True
    valid_flag=True if (int(match_len.group(1)) >= rlen*0.5) or (int(match_len.group(1)) >= len(temp[9])*0.5) else False

    return ref, [start, end], mask, int(mismatch_len.group(1)), name, temp[9], temp[10], temp[1], temp[5], primary_alignment_flag, valid_flag

def time_spend( start ):
    done = time.time()
    elapsed = done - start
    return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

def chunkify(fname, size=1*1024*1024):
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

def process_sam_file( sam_fn, numthreads, numlines ):
    result = gt._autoVivification()
    mapped_reads = 0

    #clean memory
    gc.collect()

    print_message( "Parsing SAM files with %s subprocesses..."%numthreads, argvs.silent, begin_t, logfile )
    pool = Pool(processes=numthreads)
    jobs = []
    results = []

    for chunkStart,chunkSize in chunkify(sam_fn):
        jobs.append( pool.apply_async(worker, (sam_fn,chunkStart,chunkSize)) )

    #wait for all jobs to finish
    tol_jobs = len(jobs)
    cnt=0
    for job in jobs:
        results.append( job.get() )
        cnt+=1
        if argvs.debug: print_message( "[DEBUG] Progress: %s/%s (%.1f%%) chunks done."%(cnt, tol_jobs, cnt/tol_jobs*100), argvs.silent, begin_t, logfile )

    #clean up
    pool.close()

    print_message( "Merging results...", argvs.silent, begin_t, logfile )
    for res in results:
        for k in res:
            if k in result:
                result[k]["ML"] = result[k]["ML"] | res[k]["ML"]
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
            result[k]["LL"] = 0
            mapped_reads += result[k]["MR"]
            mask = result[k]["ML"]
            del result[k]["ML"]
            bitstr = bin(mask)
            result[k]["LL"] = len(bitstr[2:].replace("0",""))

    return result, mapped_reads

def is_descendant( taxid, taxid_ant ):
    fullLineage = gt.taxid2fullLineage( taxid )
    if "|%s|" % taxid_ant in fullLineage:
        return True
    else:
        return False

def extract_read_from_sam( sam_fn, o, taxid, numthreads ):
    pool = Pool(processes=numthreads)
    jobs = []

    for chunkStart,chunkSize in chunkify(sam_fn):
        jobs.append( pool.apply_async(ReadExtractWorker, (sam_fn,chunkStart,chunkSize,taxid)) )

    #wait for all jobs to finish
    for job in jobs:
        outread = job.get()
        o.write(outread)
        o.flush()

    #clean up
    pool.close()

def ReadExtractWorker( filename, chunkStart, chunkSize, taxid ):
    # output
    readstr=""
    # processing alignments in SAM format
    f = open( filename )
    f.seek(chunkStart)
    lines = f.read(chunkSize).splitlines()
    for line in lines:
        ref, region, mask, nm, rname, rseq, rq, flag, cigr, pri_aln_flag, valid_flag = parse(line)

        if not (pri_aln_flag and valid_flag): continue

        acc, start, stop, t = ref.split('|')
        fullLineage = gt.taxid2fullLineage(t)

        if int(flag) & 16:
            g = findall(r'\d+\w', cigr)
            cigr = "".join(list(reversed(g)))
            rseq = seqReverseComplement(rseq)
            rq = rq[::-1]

        if is_descendant( t, taxid ):
            readstr += "@%s %s:%s..%s %s\n%s\n+\n%s\n" % (rname, ref, region[0], region[1], cigr, rseq, rq)
    return readstr

def seqReverseComplement( seq ):
    seq1 = 'ACGTURYSWKMBDHVNTGCAAYRSWMKVHDBNacgturyswkmbdhvntgcaayrswmkvhdbn'
    seq_dict = { seq1[i]:seq1[i+16] for i in range(64) if i < 16 or 32<=i<48 }
    return "".join([seq_dict[base] for base in reversed(seq)])

def group_refs_to_strains(r):
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
        'MB':sum, # of mapped bases
        'MR':sum, # of mapped reads
        'NM':sum, # of mismatches
        'LL':sum, # linear length
        'RLEN':sum # length of this signature fragments (mapped)
    }).reset_index()
    # total length of signatures
    str_df['TS'] = str_df['TAXID'].apply(lambda x: db_stats[x])
    str_df['bDOC'] = str_df['MB']/str_df['TS'] # bDOC: best Depth of Coverage of a strain
    str_df['bLC'] = str_df['LL']/str_df['TS'] # bLC:  best linear coverage of a strain
    str_df['RD'] = str_df['MB']/str_df['TS'] # roll-up DoC

    # rename columns
    str_df.rename(columns={
        "MB":  "TOTAL_BP_MAPPED",
        "MR":  "READ_COUNT",
        "NM":  "TOTAL_BP_MISMATCH",
        "LL":  "LINEAR_LEN",
        "RLEN":"MAPPED_SIG_LENGTH",
        "TS":  "TOL_SIG_LENGTH",
        "RD":  "ROLLUP_DOC",
        "bDOC":"BEST_DOC",
        "bLC": "BEST_LINEAR_COV"
    }, inplace=True)

    str_df['ZSCORE'] = str_df.apply(lambda x: pile_lvl_zscore(x.TOTAL_BP_MAPPED, x.TOL_SIG_LENGTH, x.LINEAR_LEN), axis=1)
    
    return str_df

def roll_up_taxonomy( r, db_stats, abu_col, tg_rank, mc, mr, ml, mz):
    """
    Take parsed SAM output and rollup to superkingdoms
    """
    major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}

    # roll up references to strains
    str_df = group_refs_to_strains(r)
    # produce columns for the final report at each ranks
    rep_df = pd.DataFrame()

    # qualified strain
    qualified_idx = (str_df['LINEAR_LEN']/str_df['TOL_SIG_LENGTH'] >= mc) & \
                    (str_df['READ_COUNT'] >= mr) & \
                    (str_df['LINEAR_LEN'] >= ml) & \
                    (str_df['ZSCORE'] <= mz)

    for rank in sorted(major_ranks, key=major_ranks.__getitem__):
        str_df['LVL_NAME'] = str_df['TAXID'].apply(lambda x: gt.taxid2lineageDICT(x, True, True)[rank]['name'])
        str_df['LVL_TAXID'] = str_df['TAXID'].apply(lambda x: gt.taxid2lineageDICT(x, True, True)[rank]['taxid'])
        str_df['LEVEL'] = rank

        # rollup strains that make cutoffs
        lvl_df = pd.DataFrame()
        if rank == 'strain':
            lvl_df = str_df
        else:
            lvl_df = str_df[qualified_idx].groupby(['LVL_NAME']).agg({
                'LEVEL':'first',
                'LVL_TAXID':'first',
                'TOTAL_BP_MAPPED': sum, 'READ_COUNT': sum, 'TOTAL_BP_MISMATCH': sum,
                'LINEAR_LEN': sum, 'MAPPED_SIG_LENGTH': sum, 'TOL_SIG_LENGTH': sum, 
                'ROLLUP_DOC': sum, 'BEST_DOC': max, 'BEST_LINEAR_COV': max, 'ZSCORE': min,
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
    filtered = (rep_df['ZSCORE'] > mz)
    rep_df.loc[filtered, 'NOTE'] += "Filtered out (minLen > " + rep_df.loc[filtered, 'ZSCORE'].astype(str) + "); "

    rep_df.drop(columns=['TAXID'], inplace=True)
    rep_df.rename(columns={"LVL_NAME": "NAME", "LVL_TAXID": "TAXID"}, inplace=True)

    return rep_df

def pile_lvl_zscore(tol_bp, genome_size, linear_len):
    try:
        doc = tol_bp/genome_size
        lin_doc = tol_bp/linear_len
        v = (linear_len*(lin_doc-doc)**2 + (genome_size-linear_len)*(doc)**2)/genome_size
        sd = math.sqrt(v)
        return (lin_doc-doc)/sd
    except:
        return 99

def generaete_taxonomy_file(rep_df, o, fullreport_o, fmt="tsv"):
    """
    output result in tsv or csv format
    """
    # Fields for full mode
    cols = ['LEVEL', 'NAME', 'TAXID', 'READ_COUNT', 'TOTAL_BP_MAPPED',
            'TOTAL_BP_MISMATCH', 'LINEAR_LEN', 'LINEAR_DOC', 'ROLLUP_DOC', 'REL_ABUNDANCE',
            'LINEAR_COV', 'LINEAR_COV_MAPPED_SIG', 'BEST_LINEAR_COV', 'MAPPED_SIG_LENGTH', 'TOL_SIG_LENGTH',
            'ABUNDANCE', 'ZSCORE', 'NOTE'
            #'PILE_NORM_D','PILE_POIS_D','PILE_N_POIS_D'
    ]

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
    output result in biom format
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

def generaete_lineage_file(res_df, o, tg_rank):
    """
    output abundance-lineage in tsv format
    """
    import csv
    target_df = pd.DataFrame()
    target_idx = (res_df['LEVEL']==tg_rank)
    target_df = res_df.loc[target_idx, ['ABUNDANCE','TAXID']]
    lineage_df = target_df['TAXID'].apply(lambda x: gt.taxid2lineage(x, True, True)).str.split('|', expand=True)
    result = pd.concat([target_df['ABUNDANCE'], lineage_df], axis=1, sort=False)
    result.to_csv(o, index=False, header=False, sep='\t', float_format='%.4f')

    return True

def readMapping(reads, db, threads, mm_penalty, presetx, samfile, logfile, nanopore):
    """
    mapping reads to database
    """
    input_file = " ".join(reads)
    
    sr_opts = f"-x {presetx} -k24 -A1 -B{mm_penalty} -O30 -E30 -a -N1 -n1 -p1 -m24 -s30"
    if nanopore:
        sr_opts = f"-x {presetx} -a"

    bash_cmd   = "set -o pipefail; set -x;"
    mm2_cmd    = f"minimap2 {sr_opts} -t{threads} {db}.mmi {input_file}"
    filter_cmd = "gawk -F\\\\t '!/^@/ && !and($2,4) && !and($2,2048) { if(r!=$1.and($2,64)){r=$1.and($2,64); s=$14} if($14>=s){print} }'"
    cmd        = "%s %s 2>> %s | %s > %s"%(bash_cmd, mm2_cmd, logfile, filter_cmd, samfile)

    proc = subprocess.Popen( cmd, shell=True, executable='/bin/bash', stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    outs, errs = proc.communicate()
    exitcode = proc.poll()

    return exitcode, mm2_cmd, errs

def loadDatabaseStats(db_stats_file):
	"""
	loading database stats from db_path.stats

	The input stats file is a 8 column tab delimited text file:
	1. Rank
	2. Name
	3. Taxid
	4. Superkingdom
	5. NumOfSeq
	6. Max
	7. Min
	8. TotalLength

	We will only save stain level info with their taxid (2) and total signature length (8).
	"""
	db_stats = {}

	with open(db_stats_file) as f:
		for line in f:
			fields = line.split("\t")
			if fields[0] == "strain": #or fields[0] == "species":
				db_stats[fields[2]] = int(fields[7])
			else:
				continue

	return db_stats

def print_message(msg, silent, start, logfile, errorout=0):
    message = "[%s] %s\n" % (time_spend(start), msg)

    with open( logfile, "a" ) as f:
        f.write( message )
        f.close()

    if errorout:
        sys.exit( message )
    elif not silent:
        sys.stderr.write( message )

if __name__ == '__main__':
    argvs    = parse_params( __version__ )
    begin_t  = time.time()
    sam_fp   = argvs.sam[0] if argvs.sam else ""
    samfile  = "%s/%s.gottcha_%s.sam" % ( argvs.outdir, argvs.prefix, argvs.dbLevel ) if not argvs.sam else sam_fp.name
    logfile  = "%s/%s.gottcha_%s.log" % ( argvs.outdir, argvs.prefix, argvs.dbLevel )
    lines_per_process = 10000

    #dependency check
    if sys.version_info < (3,4):
        sys.exit("[ERROR] Python 3.4 or above is required.")
    dependency_check("minimap2")
    dependency_check("gawk")

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

    print_message( "Starting GOTTCHA (v%s)" % __version__, argvs.silent, begin_t, logfile )
    print_message( "Arguments and dependencies checked:", argvs.silent, begin_t, logfile )
    print_message( "    Input reads      : %s" % argvs.input,     argvs.silent, begin_t, logfile )
    print_message( "    Input SAM file   : %s" % samfile,         argvs.silent, begin_t, logfile )
    print_message( "    Database         : %s" % argvs.database,  argvs.silent, begin_t, logfile )
    print_message( "    Database level   : %s" % argvs.dbLevel,   argvs.silent, begin_t, logfile )
    print_message( "    Mismatch penalty : %s" % argvs.mismatch,  argvs.silent, begin_t, logfile )
    print_message( "    Abundance        : %s" % argvs.relAbu,    argvs.silent, begin_t, logfile )
    print_message( "    Output path      : %s" % argvs.outdir,    argvs.silent, begin_t, logfile )
    print_message( "    Prefix           : %s" % argvs.prefix,    argvs.silent, begin_t, logfile )
    print_message( "    Extract taxid    : %s" % argvs.extract,   argvs.silent, begin_t, logfile )
    print_message( "    Threads          : %d" % argvs.threads,   argvs.silent, begin_t, logfile )
    print_message( "    Minimal L_DOC    : %s" % argvs.minCov,    argvs.silent, begin_t, logfile )
    print_message( "    Minimal L_LEN    : %s" % argvs.minLen,    argvs.silent, begin_t, logfile )
    print_message( "    Minimal reads    : %s" % argvs.minReads,  argvs.silent, begin_t, logfile )

    #load taxonomy
    print_message( "Loading taxonomy information...", argvs.silent, begin_t, logfile )
    custom_taxa_tsv = None
    if os.path.isfile( argvs.database + ".tax.tsv" ):
        custom_taxa_tsv = argvs.database+".tax.tsv"
    gt.loadTaxonomy( argvs.taxInfo, custom_taxa_tsv )
    print_message( "Done.", argvs.silent, begin_t, logfile )

    #load database stats
    print_message( "Loading database stats...", argvs.silent, begin_t, logfile )
    if os.path.isfile( argvs.database + ".stats" ):
        db_stats = loadDatabaseStats(argvs.database+".stats")
    else:
        sys.exit( "[%s] ERROR: %s not found.\n" % (time_spend(begin_t), argvs.database+".stats") )
    print_message( "Done.", argvs.silent, begin_t, logfile )

    if argvs.input:
        print_message( "Running read-mapping...", argvs.silent, begin_t, logfile )
        exitcode, cmd, msg = readMapping( argvs.input, argvs.database, argvs.threads, argvs.mismatch, argvs.presetx, samfile, logfile, argvs.nanopore )
        print_message( "Logfile saved to %s." % logfile, argvs.silent, begin_t, logfile )
        print_message( "COMMAND: %s" % cmd, argvs.silent, begin_t, logfile )
        if exitcode != 0:
            sys.exit( "[%s] ERROR: error occurred while running read mapping (exit: %s, message: %s).\n" % (time_spend(begin_t), exitcode, msg) )
        else:
            print_message( "Done mapping reads to %s signature database." % argvs.dbLevel, argvs.silent, begin_t, logfile )
            print_message( "Mapped SAM file saved to %s." % samfile, argvs.silent, begin_t, logfile )
            sam_fp = open( samfile, "r" )
    
    if argvs.extract:
        extract_read_from_sam( os.path.abspath(samfile), out_fp, argvs.extract, argvs.threads )
        print_message( "Done extracting reads to %s." % outfile, argvs.silent, begin_t, logfile )
    else:
        (res, mapped_r_cnt) = process_sam_file( os.path.abspath(samfile), argvs.threads, lines_per_process)
        print_message( "Done processing SAM file. %s reads mapped." % mapped_r_cnt, argvs.silent, begin_t, logfile )
        res_df = roll_up_taxonomy( res, db_stats, argvs.relAbu, argvs.dbLevel , argvs.minCov, argvs.minReads, argvs.minLen, argvs.maxZscore)
        print_message( "Done taxonomy rolling up.", argvs.silent, begin_t, logfile )

        # generate output files
        if argvs.format == "biom":
            generaete_biom_file(res_df, out_fp, argvs.dbLevel, argvs.prefix)
        else:
            generaete_taxonomy_file(res_df, out_fp, outfile_full, argvs.format)
        # generate lineage file
        generaete_lineage_file(res_df, outfile_lineage, argvs.dbLevel)

        print_message( "Done taxonomy profiling; Results saved to %s." %(outfile), argvs.silent, begin_t, logfile )