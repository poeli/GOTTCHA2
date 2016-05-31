#!/usr/bin/env python

"""
[AUTHOR]
Po-E (Paul) Li
Bioinformatics Team
Los Alamos National Laboratory
2016/05/31

[COPYRIGHT]
Copyright (2014). Los Alamos National Security, LLC. This material was produced
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
import sys, time, os.path as op
import gottcha_taxonomy as gt
from subprocess import getstatusoutput
from re import search
from multiprocessing import Pool
from itertools import chain

def parse_params( ver ):
	p = ap.ArgumentParser( prog='gottcha.py', description="""Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an
			annotation-independent and signature-based metagenomic taxonomic profiling tool
			that has significantly smaller FDR than other profiling tools. This program
			is a wrapper to map input reads to pre-computed signature databases using BWA-MEM
			and/or to profile mapped reads in SAM format. (VERSION: %s)""" % ver)

	eg = p.add_mutually_exclusive_group( required=True )

	eg.add_argument( '-i','--input', metavar='[FASTQ]', nargs='+', type=str,
	  				help="Input one or multiple FASTQ file(s). Use space to separate multiple input files.")

	eg.add_argument( '-s','--sam', metavar='[SAMFILE]', nargs=1, type=ap.FileType('r'),
					help="Specify the input SAM file. Use '-' for standard input.")

	p.add_argument( '-d','--database', metavar='[BWA_INDEX]', type=str, default=None,
	                help="The path of signature database. The database can be in FASTA format or BWA index (5 files).")

	p.add_argument( '-l','--dbLevel', metavar='[LEVEL]', type=str, default='',
	                choices=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'],
	                help="""Specify the taxonomic level of the input database. You can choose one rank from "superkingdom", "phylum", "class", "order", "family", "genus", "species" and "strain". The value will be auto-detected if the input database ended with levels (e.g. GOTTCHA_db.species).""")

	p.add_argument( '-pm','--mismatch', metavar='<INT>', type=int, default=5,
					help="Mismatch penalty for BWA-MEM (pass to option -B while BWA-MEM is running). You can use 99 for not allowing mismatch in alignments (except for extreme cases). [default: 5]")

	p.add_argument( '-m','--mode', type=str, default='summary',
	                choices=['summary', 'full', 'tree', 'class', 'extract'],
					help="""You can specify one of the following output modes:
                            "summary" : report a summary of profiling result;
                            "full"    : other than a summary result, this mode will report unfiltered profiling results with more detail;
                            "tree"    : report results with lineage of taxonomy;
                            "class"   : output results of classified reads;
                            "extract" : extract mapped reads;
						  Note that only results/reads belongs to descendants of TAXID will be reported/extracted if option [--taxonomy TAXID] is specified. [default: summary]""" )

	p.add_argument( '-x','--taxonomy', metavar='[TAXID]', type=str,
	        		help="Specify a NCBI taxonomy ID. The program  will only report/extract the taxonomy you specified.")

	p.add_argument( '-r','--relAbu', metavar='[FIELD]', type=str, default='LINEAR_DOC',
					choices=['LINEAR_LENGTH','TOTAL_BP_MAPPED','READ_COUNT','LINEAR_DOC'],
					help='The field will be used to calculate relative abundance. You can specify one of the following fields: "LINEAR_LENGTH", "TOTAL_BP_MAPPED", "READ_COUNT" and "LINEAR_DOC". [default: LINEAR_DOC]')

	p.add_argument( '-t','--threads', metavar='<INT>', type=int, default=1,
					help="Number of threads [default: 1]")

	p.add_argument( '-o','--outdir', metavar='[DIR]', type=str, default='.',
					help="Output directory [default: .]")

	p.add_argument( '-p','--prefix', metavar='<STR>', type=str, required=False,
					help="Prefix of the output file [default: <INPUT_FILE_PREFIX>]")

	p.add_argument( '-mc','--minCov', metavar='<FLOAT>', type=float, default=0.005,
					help="Minimum linear coverage to be considered valid in abundance calculation [default: 0.005]")

	p.add_argument( '-mr','--minReads', metavar='<INT>', type=int, default=3,
					help="Minimum number of reads to be considered valid in abundance calculation [default: 3]")

	p.add_argument( '-ml','--minLen', metavar='<INT>', type=int, default=60,
					help="Minimum unique length to be considered valid in abundance calculation [default: 60]")

	p.add_argument( '-c','--stdout', action="store_true",
					help="Write on standard output.")

	p.add_argument( '-v','--verbose', action="store_true",
	                help="Verbose output")

	args_parsed = p.parse_args()

	"""
	Checking options
	"""
	if args_parsed.input and not args_parsed.database:
		p.error( '--database option is missing.' )

	if args_parsed.input and args_parsed.sam:
		p.error( '--input and --same are incompatible options.' )

	if not args_parsed.prefix:
		if args_parsed.input:
			name = search('([^\/\.]+)\..*$', args_parsed.input[0] )
			args_parsed.prefix = name.group(1)
		elif args_parsed.sam:
			name = search('([^\/]+).\w+.\w+$', args_parsed.sam[0].name )
			args_parsed.prefix = name.group(1)
		else:
			args_parsed.prefix = "GOTTCHA_"

	if not args_parsed.dbLevel:
		if args_parsed.database:
			name = search('\.(\w+).\w+$', args_parsed.database )
			args_parsed.dbLevel = name.group(1)
		elif args_parsed.sam:
			name = search('\.gottcha_(\w+).sam$', args_parsed.sam[0].name )
			args_parsed.dbLevel = name.group(1)
		else:
			p.error( '--dbLevel is missing and cannot be auto-detected.' )

	return args_parsed

def join_ranges(data, offset=0):
	"""
	Merge location regions
	"""
	flatten = chain.from_iterable
	LEFT, RIGHT = 1, -1
	data = sorted(flatten(((start, LEFT), (stop + offset, RIGHT)) for start, stop in data))
	c = 0
	for value, label in data:
		if c == 0:
			x = value
		c += label
		if c == 0:
			yield [x, value - offset]

def worker(lines):
	"""Make a dict out of the parsed, supplied lines"""
	res = gt._autoVivification()

	for line in lines:
		k, r, n, rd, rs, rq, flag, cigr = parse(line)
		if k in res:
			res[k]["ML"] = list( join_ranges( res[k]["ML"] + [r] ) )
			res[k]["MB"] += r[1] - r[0] + 1
			res[k]["MR"] += 1
			res[k]["NM"] += n
		else:
			res[k]["ML"] = [r]
			res[k]["MB"] = r[1] - r[0] + 1
			res[k]["MR"] = 1
			res[k]["NM"] = n

	return res

def parse(line):
	"""
	Parse SAM format
	read1   0   test    11  0   5S10M3S *   0   0   GGGGGCCCCCCCCCCGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
	read2   16  test    11  0   3S10M5S *   0   0   GGGCCCCCCCCCCGGGGG  HHHHHHHHHHHHHHHHHH  NM:i:0  MD:Z:10 AS:i:10 XS:i:0
	"""
	temp = line.split('\t')

	name = temp[0]
	match_len    = search('(\d+)M', temp[5])
	mismatch_len = search('NM:i:(\d+)', temp[11])
	start = int(temp[3])
	end   = start + int(match_len.group(1)) - 1;

	return temp[2], [start, end], int(mismatch_len.group(1)), name, temp[9], temp[10], temp[1], temp[5]

def timeSpend( start ):
	done = time.time()
	elapsed = done - start
	return time.strftime( "%H:%M:%S", time.gmtime(elapsed) )

def processSAMfile( sam_fp, numthreads, numlines ):
	result = gt._autoVivification()
	mapped_reads = 0

	lines = sam_fp.readlines();
	if len(lines) > numlines and numthreads > 2:
		pool = Pool(processes=numthreads)
		result_list = pool.map(worker,
		    (lines[line:line+numlines] for line in range(0,len(lines),numlines) ) )

		for res in result_list:
			for k in res:
				if k in result:
					result[k]["ML"] = list( join_ranges( result[k]["ML"] + res[k]["ML"] ) )
					result[k]["MB"] += res[k]["MB"]
					result[k]["MR"] += res[k]["MR"]
					result[k]["NM"] += res[k]["NM"]
				else:
					result[k].update(res[k])

	else:
		result = worker(lines)

	# convert mapped regions to linear length
	for k in result:
		result[k]["LL"] = 0
		mapped_reads += result[k]["MR"]
		for cr in result[k]["ML"]:
			result[k]["LL"] += cr[1]-cr[0]

	return result, mapped_reads

def isDescendant( taxid, taxid_ant ):
	fullLineage = gt.taxid2fullLineage( taxid )
	if "|%s|" % taxid_ant in fullLineage:
		return True
	else:
		return False

def processSAMfileReadClass( f, o, tg_rank, taxid_fi ):
	for line in f:
		ref, region, nm, rname, rseq, rq, flag, cigr = parse(line)
		acc, start, stop, tid = ref.split('|')

		if taxid_fi:
			if not isDescendant( tid, taxid_fi ):
				continue

		o.write( "%s\t%s\t%s\t%s\t%s\n" % (
		    rname,
		    gt.taxid2taxidOnRank(tid, tg_rank),
			"%s:%s..%s" % (ref, region[0], region[1]),
			cigr,
		    gt.taxid2nameOnRank( tid, tg_rank )
		))

def processSAMfileReadExtract( f, o, taxid ):
	for line in f:
		ref, region, nm, rname, rseq, rq, flag, cigr = parse(line)
		acc, start, stop, t = ref.split('|')
		fullLineage = gt.taxid2fullLineage(t)

		if int(flag) & 16:
			rseq = seqReverseComplement(rseq)
			rq = rq[::-1]

		if isDescendant( t, taxid ):
			o.write( "@%s %s:%s..%s %s\n%s\n+\n%s\n" % (rname, ref, region[0], region[1], cigr, rseq, rq) )

def seqReverseComplement( seq ):
	for base in seq:
		if base not in 'ACGTURYSWKMBDHVNacgturyswkmbdhvn':
			print("[ERROR] NOT a DNA sequence")
			return None
	seq1 = 'ACGTURYSWKMBDHVNTGCAAYRSWMKVHDBNacgturyswkmbdhvntgcaayrswmkvhdbn'
	seq_dict = { seq1[i]:seq1[i+16] for i in range(64) if i < 16 or 32<=i<48 }
	return "".join([seq_dict[base] for base in reversed(seq)])

def taxonomyRollUp( r ):
	"""
	Take parsed SAM output and rollup to superkingdoms
	"""
	res_rollup = gt._autoVivification()
	res_tree = gt._autoVivification()

	for ref in r:
		(acc, start, stop, t) = ref.split('|')

		tree = gt.taxid2fullLinkDict( t )

		for pid, tid in tree.items():
			res_tree[pid][tid] = 1
			if tid in res_rollup:
				res_rollup[tid]["ML"] += ";%s:%s" %  ( ref, ",".join("..".join(map(str,l)) for l in r[ref]["ML"]) )
				res_rollup[tid]["MB"] += r[ref]["MB"]
				res_rollup[tid]["MR"] += r[ref]["MR"]
				res_rollup[tid]["NM"] += r[ref]["NM"]
				res_rollup[tid]["LL"] += r[ref]["LL"]
				res_rollup[tid]["SL"] += int(stop) - int(start) + 1
			else:
				res_rollup[tid]["ML"] = "%s:%s" %  ( ref, ",".join("..".join(map(str,l)) for l in r[ref]["ML"]) )
				res_rollup[tid]["MB"] = r[ref]["MB"]
				res_rollup[tid]["MR"] = r[ref]["MR"]
				res_rollup[tid]["NM"] = r[ref]["NM"]
				res_rollup[tid]["LL"] = r[ref]["LL"]
				res_rollup[tid]["SL"] = int(stop) - int(start) + 1

	return res_rollup, res_tree

def outputResultsAsTree( tid, res_tree, res_rollup, indent, taxid_fi, o, mc, mr, ml ):
	"""
	iterate taxonomy tree and print results recursively
	"""
	if int(tid) == 1:
		o.write( "%s%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % ( "", "NAME", "LEVEL", "READ_COUNT", "TOTAL_BP_MAPPED", "TOTAL_BP_MISMATCH", "LINEAR_LENGTH", "LINEAR_DOC" ) )
	else:
		if mc <= int(res_rollup[tid]["MB"])/int(res_rollup[tid]["LL"]) and mr <= int(res_rollup[tid]["MR"]) and ml <= res_rollup[tid]["LL"]:
			o.write( "%s%s\t%s\t%s\t%s\t%s\t%s\t%.4f\n" % (
				indent,
				gt.taxid2name(tid),
				gt.taxid2rank(tid),
				res_rollup[tid]["MR"],
				res_rollup[tid]["MB"],
				res_rollup[tid]["NM"],
				res_rollup[tid]["LL"],
				int(res_rollup[tid]["MB"])/int(res_rollup[tid]["LL"])
				)
			)

	if len( res_tree[tid] ):
		indent += "    "
		for cid in res_tree[tid]:
			if taxid_fi:
				if not isDescendant( tid, taxid_fi ):
					continue

			outputResultsAsTree( cid, res_tree, res_rollup, indent, taxid_fi, o, mc, mr, ml )

def outputResultsAsRanks( res_rollup, o, tg_rank, relAbu, mode, mc, mr, ml ):
	output = gt._autoVivification()
	major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}

	for tid in res_rollup:
		rank = gt.taxid2rank(tid)
		if rank in major_ranks and major_ranks[rank] <= major_ranks[tg_rank]:
			if relAbu == "LINEAR_LENGTH":
				abundance = int(res_rollup[tid]["LL"])
			elif relAbu == "TOTAL_BP_MAPPED":
				abundance = int(res_rollup[tid]["MB"])
			elif relAbu == "READ_COUNT":
				abundance = int(res_rollup[tid]["MR"])
			else:
				abundance = int(res_rollup[tid]["MB"])/int(res_rollup[tid]["LL"])

			output[rank]["RES"][tid] = abundance

			if "TOT_ABU" in output[rank]:
				output[rank]["TOT_ABU"] += abundance
			else:
				output[rank]["TOT_ABU"] = abundance

	add_field = "\tLINEAR_COV\tNOTE\tMAPPED_REGION" if mode == "full" else ""

	o.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n" %
			( "LEVEL", "NAME", "TAXID", "READ_COUNT", "TOTAL_BP_MAPPED",
			  "TOTAL_BP_MISMATCH", "LINEAR_LENGTH", "LINEAR_DOC", "REL_ABUNDANCE", add_field ) )

	for rank in sorted( major_ranks, key=major_ranks.__getitem__ ):
		if major_ranks[rank] > major_ranks[tg_rank]: break
		taxas = output[rank]["RES"]
		for tid in sorted( taxas, key=taxas.__getitem__, reverse=True):
			mc_note = "Filtered out for minCov. "   if mc > res_rollup[tid]["LL"]/res_rollup[tid]["SL"] else ""
			mr_note = "Filtered out for minReads. " if mr > int(res_rollup[tid]["MR"]) else ""
			ml_note = "Filtered out for minLen. "   if ml > int(res_rollup[tid]["LL"]) else ""

			add_field = "\t%.4f\t%s%s%s\t%s" % (
			    (res_rollup[tid]["LL"]/res_rollup[tid]["SL"]),
				mc_note,
				mr_note,
				ml_note,
				res_rollup[tid]["ML"]
			) if mode == "full" else ""

			if mc_note and mr_note and ml_note and mode=="summary":
				continue

			o.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\t%.4f%s\n" %
				( rank,
				  gt.taxid2name(tid),
				  tid,
				  res_rollup[tid]["MR"],
				  res_rollup[tid]["MB"],
				  res_rollup[tid]["NM"],
				  res_rollup[tid]["LL"],
				  int(res_rollup[tid]["MB"])/int(res_rollup[tid]["LL"]),
				  int(res_rollup[tid]["MB"])/int(res_rollup[tid]["LL"])/output[rank]["TOT_ABU"],
				  add_field
				)
			)

def readMapping( reads, db, threads, mm_penalty, samfile, logfile ):
	"""
	mapping reads to database
	"""
	input_file = " ".join(reads)
	bwa_cmd = "bwa mem -k30 -T30 -A1 -B%s -O99 -E99 -L0 -P -S -t%s %s %s | samtools view -F4 > %s 2> %s" % ( mm_penalty, threads, db, input_file, samfile, logfile )
	exitcode, msg = getstatusoutput( bwa_cmd )
	return exitcode, bwa_cmd

if __name__ == '__main__':
	verison  = "2.0 BETA"
	argvs    = parse_params( verison )
	start    = time.time()
	numlines = 10000
	sam_fp   = argvs.sam[0] if argvs.sam else ""
	samfile  = "%s/%s.gottcha_%s.sam" % ( argvs.outdir, argvs.prefix, argvs.dbLevel ) if not argvs.sam else sam_fp.name
	logfile  = "%s/%s.gottcha_%s.log" % ( argvs.outdir, argvs.prefix, argvs.dbLevel )

	sys.stderr.write( "[%s] Starting GOTTCHA (v%s)\n" % (timeSpend(start), verison) )

	#prepare output object
	out_fp = sys.stdout
	outfile = "STDOUT"
	if not argvs.stdout:
		ext = "fastq" if argvs.mode == "reads" else "tsv"
		tg_taxid = argvs.taxonomy if argvs.taxonomy else ""
		outfile = "%s/%s.%s%s.%s" % ( argvs.outdir, argvs.prefix, argvs.mode, tg_taxid, ext)
		out_fp = open( outfile, 'w')

	verbose_msg = """[%s] Arguments checked:
           Input reads    : %s
           Input SAM file : %s
           Database       : %s
           Database level : %s
           Output path    : %s
           Prefix         : %s
           Mode           : %s
           Specific taxid : %s
           Threads        : %d
           Minimal L_DOC  : %s
           Minimal L_LEN  : %s
           Minimal reads  : %s\n""" % (
		timeSpend(start), argvs.input, samfile, argvs.database, argvs.dbLevel,
		argvs.outdir, argvs.prefix, argvs.mode, argvs.taxonomy, argvs.threads, argvs.minCov, argvs.minLen, argvs.minReads
	) if argvs.verbose else ""
	sys.stderr.write( verbose_msg )

	#load taxonomy
	gt.loadTaxonomy()
	verbose_msg = "[%s] Taxonomy information loaded.\n" % timeSpend(start) if argvs.verbose else ""
	sys.stderr.write( verbose_msg )

	#load strain
	gt.loadStrainName()
	verbose_msg = "[%s] Non-standard groups/strains information loaded.\n" % timeSpend(start) if argvs.verbose else ""
	sys.stderr.write( verbose_msg )

	if argvs.input:
		verbose_msg = "[%s] Running read-mapping...\n" % timeSpend(start) if argvs.verbose else ""
		sys.stderr.write( verbose_msg )
		exitcode, cmd = readMapping( argvs.input, argvs.database, argvs.threads, argvs.mismatch, samfile, logfile )
		verbose_msg  = "[%s] Logfile saved to %s.\n" % (timeSpend(start), logfile) if argvs.verbose else ""
		verbose_msg += "[%s] COMMAND: %s\n" % (timeSpend(start), cmd) if argvs.verbose else ""
		sys.stderr.write( verbose_msg )

		if exitcode:
			sys.exit( "[%s] ERROR: error occurred while running read mapping (exit code: %s).\n" % (timeSpend(start), exitcode) )
		else:
			sys.stderr.write( "[%s] Done mapping reads to %s signature database.\n" % (timeSpend(start), argvs.dbLevel) )
			sam_fp = open( samfile, "r" )

	if argvs.mode == 'class':
		processSAMfileReadClass( sam_fp, out_fp, argvs.dbLevel, argvs.taxonomy )
		sys.stderr.write( "[%s] Done classifying reads. Results saved to %s.\n" % (timeSpend(start), outfile) )

	elif argvs.mode == 'extract':
		processSAMfileReadExtract( sam_fp, out_fp, argvs.taxonomy )
		sys.stderr.write( "[%s] Done extracting reads to %s.\n" % timeSpend(start), outfile )

	else:
		(res, mapped_r_cnt) = processSAMfile( sam_fp, argvs.threads, numlines)
		sys.stderr.write( "[%s] Done processing SAM file. %s reads mapped.\n" % (timeSpend(start), mapped_r_cnt) )

		(res_rollup, res_tree) = taxonomyRollUp( res )
		sys.stderr.write( "[%s] Done taxonomy rolling up.\n" % timeSpend(start) )

		if argvs.mode == 'summary' or argvs.mode == 'full':
			outputResultsAsRanks( res_rollup, out_fp, argvs.dbLevel, argvs.relAbu, argvs.mode, argvs.minCov, argvs.minReads, argvs.minLen )
		elif argvs.mode == 'tree':
			outputResultsAsTree( "1", res_tree, res_rollup, "", argvs.taxonomy, out_fp, argvs.minCov, argvs.minReads, argvs.minLen )

		sys.stderr.write( "[%s] Done taxonomy preofiling; %s results saved to %s.\n" % (timeSpend(start), argvs.mode, outfile) )
