#!/usr/bin/env python3

__author__    = "Po-E (Paul) Li, Bioscience Division, Los Alamos National Laboratory"
__credits__   = ["Po-E Li", "Jason Gans", "Tracey Freites", "Patrick Chain"]
__version__   = "2.1.4 BETA"
__date__      = "2018/10/07"
__copyright__ = """
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
import sys, os, time, subprocess
import taxonomy as gt
import gc
from re import search,findall
from re import compile as recompile
from multiprocessing import Pool
from itertools import chain

def parse_params( ver ):
	p = ap.ArgumentParser( prog='gottcha.py', description="""Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an
			annotation-independent and signature-based metagenomic taxonomic profiling tool
			that has significantly smaller FDR than other profiling tools. This program
			is a wrapper to map input reads to pre-computed signature databases using minimap2
			and/or to profile mapped reads in SAM format. (VERSION: %s)""" % ver)

	eg = p.add_mutually_exclusive_group( required=True )

	eg.add_argument( '-i','--input', metavar='[FASTQ]', nargs='+', type=str,
					help="Input one or multiple FASTQ file(s). Use space to separate multiple input files.")

	eg.add_argument( '-s','--sam', metavar='[SAMFILE]', nargs=1, type=ap.FileType('r'),
					help="Specify the input SAM file. Use '-' for standard input.")

	p.add_argument( '-d','--database', metavar='[MINIMAP2_INDEX]', type=str, default=None,
					help="The path of signature database. The database can be in FASTA format or minimap2 index (5 files).")

	p.add_argument( '-l','--dbLevel', metavar='[LEVEL]', type=str, default='',
					choices=['superkingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain'],
					help="""Specify the taxonomic level of the input database. You can choose one rank from "superkingdom", "phylum", "class", "order", "family", "genus", "species" and "strain". The value will be auto-detected if the input database ended with levels (e.g. GOTTCHA_db.species).""")

	p.add_argument( '-ti','--taxInfo', metavar='[FILE]', type=str, default='',
					help="""Specify the path of taxonomy information file (taxonomy.tsv). GOTTCHA2 will try to locate this file when user doesn't specify a path. If '--database' option is used, the program will try to find this file in the directory of specified database. If not, the 'database' directory under the location of gottcha.py will be used as default.""")

	p.add_argument( '-pm','--mismatch', metavar='<INT>', type=int, default=5,
					help="Mismatch penalty for BWA-MEM (pass to option -B while BWA-MEM is running). You can use 99 for not allowing mismatch in alignments (except for extreme cases). [default: 5]")

	p.add_argument( '-m','--mode', type=str, default='summary',
					choices=['summary', 'full', 'class', 'extract', 'lineage'],
					help="""You can specify one of the following output modes:
							"summary" : report a summary of profiling result;
							"full"    : other than a summary result, this mode will report unfiltered profiling results with more detail;
							"class"   : output results of classified reads;
							"extract" : extract mapped reads;
							"lineage" : output abundance and lineage in a line;
						  Note that only results/reads belongs to descendants of TAXID will be reported/extracted if option [--taxonomy TAXID] is specified. [default: summary]""" )

	p.add_argument( '-x','--taxonomy', metavar='[TAXID]', type=str,
					help="Specify a NCBI taxonomy ID. The program  will only report/extract the taxonomy you specified.")

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

	p.add_argument( '-mh','--minMLRL', metavar='<INT>', type=float, default=1,
					help="Minimum mean linear read length [default: 1]")

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
			db_dir = search( '^(.*?)[^\/]+$', args_parsed.database )
			args_parsed.taxInfo = db_dir.group(1)
		else:
			bin_dir = os.path.dirname(os.path.realpath(__file__))
			args_parsed.taxInfo = bin_dir + "/database"

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
			major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7, "strain":8}
			parts = args_parsed.database.split('.')
			for part in parts:
				if part in major_ranks:
					args_parsed.dbLevel = part
					break
		elif args_parsed.sam:
			name = search('\.gottcha_(\w+).sam$', args_parsed.sam[0].name )
			args_parsed.dbLevel = name.group(1)
		else:
			p.error( '--dbLevel is missing and cannot be auto-detected.' )

	if args_parsed.noCutoff:
		#p.error( 'conflict options: cutoff(s) are specified with --noCutoff option.' )
		args_parsed.minCov = 0
		args_parsed.minReads = 0
		args_parsed.minLen = 0

	return args_parsed

def dependency_check(cmd):
	proc = subprocess.Popen("which " + cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
	outs, errs = proc.communicate()
	return outs.decode().rstrip() if proc.returncode == 0 else False

def worker(filename, chunkStart, chunkSize):
	"""Make a dict out of the parsed, supplied lines"""
	# processing alignments in SAM format
	f = open( filename )
	f.seek(chunkStart)
	lines = f.read(chunkSize).splitlines()
	res={}

	for line in lines:
		k, r, m, n, rd, rs, rq, flag, cigr, pri_aln_flag, valid_flag = parse(line)
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
	match_len    = search('(\d+)M', temp[5])
	mismatch_len = search('NM:i:(\d+)', temp[11])
	start = int(temp[3])
	end   = start + int(match_len.group(1)) - 1

	ref = temp[2].rstrip('|')
	ref = ref[: -2 if ref.endswith(".0") else None ]

	(acc, rstart, rend, taxid) = ref.split('|')
	rlen = int(rend)-int(rstart)+1
	mask = int( "%s%s%s"%("0"*(start-1), "1"*(end-start+1), "0"*(rlen-end)), 2)

	primary_alignment_flag=False if int(temp[1]) & 256 else True
	valid_flag=True if (int(match_len.group(1)) >= rlen*0.5) or (int(match_len.group(1)) >= len(temp[9])*0.5) else False

	return ref, [start, end], mask, int(mismatch_len.group(1)), name, temp[9], temp[10], temp[1], temp[5], primary_alignment_flag, valid_flag

def timeSpend( start ):
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

def processSAMfile( sam_fn, numthreads, numlines ):
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
			p = recompile('1+')
			bitstr = bin(mask).replace('0b','')
			iterator = p.finditer(bitstr)
			for match in iterator:
				r = match.span()
				result[k]["LL"] += r[1]-r[0]

	return result, mapped_reads

def isDescendant( taxid, taxid_ant ):
	fullLineage = gt.taxid2fullLineage( taxid )
	if "|%s|" % taxid_ant in fullLineage:
		return True
	else:
		return False

def processSAMfileReadClass( f, o, tg_rank, taxid_fi ):
	for line in f:
		ref, region, mask, nm, rname, rseq, rq, flag, cigr, pri_aln_flag, valid_flag = parse(line)
		acc, start, stop, tid = ref.split('|')

		if taxid_fi:
			if not isDescendant( tid, taxid_fi ):
				continue

		if pri_aln_flag and valid_flag:
			o.write( "%s\t%s\t%s\t%s\t%s\n" % (
				rname,
				gt.taxid2taxidOnRank(tid, tg_rank),
				"%s:%s..%s" % (ref, region[0], region[1]),
				cigr,
				gt.taxid2nameOnRank( tid, tg_rank )
			))

def processSAMfileReadExtract( sam_fn, o, taxid, numthreads ):
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

		if not (pri_aln_flag and valid_flag): next

		acc, start, stop, t = ref.split('|')
		fullLineage = gt.taxid2fullLineage(t)

		if int(flag) & 16:
			g = findall("\d+\w", cigr)
			cigr = "".join(list(reversed(g)))
			rseq = seqReverseComplement(rseq)
			rq = rq[::-1]

		if isDescendant( t, taxid ):
			readstr += "@%s %s:%s..%s %s\n%s\n+\n%s\n" % (rname, ref, region[0], region[1], cigr, rseq, rq)
	return readstr

def seqReverseComplement( seq ):
	for base in seq:
		if base not in 'ACGTURYSWKMBDHVNacgturyswkmbdhvn':
			print("[ERROR] NOT a DNA sequence")
			return None
	seq1 = 'ACGTURYSWKMBDHVNTGCAAYRSWMKVHDBNacgturyswkmbdhvntgcaayrswmkvhdbn'
	seq_dict = { seq1[i]:seq1[i+16] for i in range(64) if i < 16 or 32<=i<48 }
	return "".join([seq_dict[base] for base in reversed(seq)])

def taxonomyRollUp( r, db_stats, relAbu, mc, mr, ml, mh ):
	"""
	Take parsed SAM output and rollup to superkingdoms
	"""
	res_rollup = gt._autoVivification()
	res_tree = gt._autoVivification()
	major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7}

	# rollup to strain first
	for ref in r:
		(acc, start, stop, stid) = ref.split('|')
		if stid in res_rollup:
			# ML: mapped region
			# MB: # of mapped bases
			# MR: # of mapped reads
			# NM: # of mismatches
			# LL: linear length
			# SL: length of this signature fragments (mapped)
			# TS: length of total signature fragments for a strain (mapped + unmapped)
			#res_rollup[stid]["ML"] += ";%s:%s" %  ( ref, ",".join("..".join(map(str,l)) for l in r[ref]["ML"]) )
			res_rollup[stid]["MB"] += r[ref]["MB"]
			res_rollup[stid]["MR"] += r[ref]["MR"]
			res_rollup[stid]["NM"] += r[ref]["NM"]
			res_rollup[stid]["LL"] += r[ref]["LL"]
			res_rollup[stid]["SL"] += int(stop) - int(start) + 1
		else:
			#res_rollup[stid]["ML"] = "%s:%s" %  ( ref, ",".join("..".join(map(str,l)) for l in r[ref]["ML"]) )
			res_rollup[stid]["MB"] = r[ref]["MB"]
			res_rollup[stid]["MR"] = r[ref]["MR"]
			res_rollup[stid]["NM"] = r[ref]["NM"]
			res_rollup[stid]["LL"] = r[ref]["LL"]
			res_rollup[stid]["SL"] = int(stop) - int(start) + 1
			res_rollup[stid]["TS"] = db_stats[stid]

	# get all strain tax id
	allStrTaxid = list(res_rollup)

	# Calculating DOC, LC and CC for strains
	# These calculations need to be done before rollup step because
	# it's possible that a strain's parent is a strain (no rank) as well
	for stid in allStrTaxid:
		res_rollup[stid]["bDOC"] = res_rollup[stid]["MB"]/db_stats[stid]
		res_rollup[stid]["bLC"]  = res_rollup[stid]["LL"]/db_stats[stid]
		res_rollup[stid]["RD"]   = res_rollup[stid]["MB"]/db_stats[stid]

	# roll strain results to upper levels
	for stid in allStrTaxid:
		# apply cutoffs strain level and rollup to higher levels		
		if mc > res_rollup[stid]["LL"]/db_stats[stid] or \
			mr > res_rollup[stid]["MR"] or \
			ml > res_rollup[stid]["LL"] or \
			mh > res_rollup[stid]["LL"]/res_rollup[stid]["MR"]:
			continue

		tree = gt.taxid2fullLinkDict( stid )

		for pid, tid in tree.items():
			res_tree[pid][tid] = 1
			if tid == stid: # skip strain id, rollup only
				continue
			if not gt.taxid2rank(tid) in major_ranks:
				continue
			if tid in res_rollup:
				# bDOC: best Depth of Coverage of a strain
				# bLC:  best linear coverage of a strain
				#res_rollup[tid]["ML"]   += ";%s" % res_rollup[stid]["ML"]
				res_rollup[tid]["MB"]   += res_rollup[stid]["MB"]
				res_rollup[tid]["MR"]   += res_rollup[stid]["MR"]
				res_rollup[tid]["NM"]   += res_rollup[stid]["NM"]
				res_rollup[tid]["LL"]   += res_rollup[stid]["LL"]
				res_rollup[tid]["SL"]   += res_rollup[stid]["SL"]
				res_rollup[tid]["TS"]   += res_rollup[stid]["TS"]
				res_rollup[tid]["RD"]   += res_rollup[stid]["RD"]
				res_rollup[tid]["bDOC"]  = res_rollup[stid]["bDOC"] if res_rollup[stid]["bDOC"] > res_rollup[tid]["bDOC"] else res_rollup[tid]["bDOC"]
				res_rollup[tid]["bLC"]   = res_rollup[stid]["bLC"] if res_rollup[stid]["bLC"] > res_rollup[tid]["bLC"] else res_rollup[tid]["bLC"]
			else:
				#res_rollup[tid]["ML"]    = res_rollup[stid]["ML"]
				res_rollup[tid]["MB"]    = res_rollup[stid]["MB"]
				res_rollup[tid]["MR"]    = res_rollup[stid]["MR"]
				res_rollup[tid]["NM"]    = res_rollup[stid]["NM"]
				res_rollup[tid]["LL"]    = res_rollup[stid]["LL"]
				res_rollup[tid]["SL"]    = res_rollup[stid]["SL"]
				res_rollup[tid]["TS"]    = res_rollup[stid]["TS"]
				res_rollup[tid]["RD"]    = res_rollup[stid]["RD"]
				res_rollup[tid]["bDOC"]  = res_rollup[stid]["bDOC"]
				res_rollup[tid]["bLC"]   = res_rollup[stid]["bLC"]

	#add abundance to res_rollup
	for tid in res_rollup:
		if relAbu == "LINEAR_LENGTH":
			res_rollup[tid]["ABU"] = res_rollup[tid]["LL"]
		elif relAbu == "TOTAL_BP_MAPPED":
			res_rollup[tid]["ABU"] = res_rollup[tid]["MB"]
		elif relAbu == "READ_COUNT":
			res_rollup[tid]["ABU"] = res_rollup[tid]["MR"]
		elif relAbu == "LINEAR_DOC":
			res_rollup[tid]["ABU"] = res_rollup[tid]["MB"]/res_rollup[tid]["LL"]
		else:
			res_rollup[tid]["ABU"] = res_rollup[tid]["RD"]

	return res_rollup, res_tree

def outputResultsAsRanks( res_rollup, o, tg_rank, mode, mc, mr, ml, mh ):
	output = gt._autoVivification()
	major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7,"strain":8}

	# init total abundance
	tol_abu = {}
	tol_abu["ROLLUP_DOC"] = 0
	tol_abu["LINEAR_DOC"] = 0
	tol_abu["READ_COUNT"] = 0
	tol_abu["TOTAL_BP_MAPPED"] = 0
	tol_abu["ABU"] = 0

	# calculate total abundances and prepare dictionary using ranks as keys
	for tid in res_rollup:
		rank = gt.taxid2rank(tid)
		if rank == "superkingdom":
			tol_abu["ROLLUP_DOC"]      += res_rollup[tid]["RD"]
			tol_abu["READ_COUNT"]      += res_rollup[tid]["MR"]
			tol_abu["TOTAL_BP_MAPPED"] += res_rollup[tid]["MB"]
			tol_abu["ABU"]             += res_rollup[tid]["ABU"]

		if rank in major_ranks:
			if not rank in output:
				output[rank] = []
			output[rank].append(tid)

	# Fields for full mode
	add_field = "\t" + "\t".join([
			"LINEAR_COV",
			"LINEAR_COV_MAPPED_SIG",
			"BEST_LINEAR_COV",
			"DOC",
			"BEST_DOC",
			"MAPPED_SIG_LENGTH",
			"TOL_SIG_LENGTH",
			"ABUNDANCE",
			"REL_ABU_ROLLUP_DOC",
			"REL_ABU_READ_COUNT",
			"REL_ABU_TOL_BP_MAPPED",
			"MLRL",
			"NOTE" ]) if mode == "full" else ""

	# essential fields
	o.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s%s\n" % (
			"LEVEL",
			"NAME",
			"TAXID",
			"READ_COUNT",
			"TOTAL_BP_MAPPED",
			"TOTAL_BP_MISMATCH",
			"LINEAR_LENGTH",
			"LINEAR_DOC",
			"ROLLUP_DOC",
			"REL_ABUNDANCE", add_field ) )

	for rank in sorted( major_ranks, key=major_ranks.__getitem__ ):
		if major_ranks[rank] > major_ranks[tg_rank] and mode == "summary":
			break

		for tid in sorted( output[rank], key=lambda tid: res_rollup[tid]["ABU"], reverse=True):
			note = ""
			note += "Filtered out (minCov > %.2f); "%(res_rollup[tid]["LL"]/db_stats[tid]) if rank == "strain" and tid in db_stats and mc > res_rollup[tid]["LL"]/db_stats[tid] else ""
			note += "Filtered out (minReads > %s); "%res_rollup[tid]["MR"] if mr > int(res_rollup[tid]["MR"]) else ""
			note += "Filtered out (minLen > %s); "%res_rollup[tid]["LL"] if ml > int(res_rollup[tid]["LL"]) else ""
			note += "Filtered out (minMLRL > %.2f); "%(res_rollup[tid]["LL"]/res_rollup[tid]["MR"]) if mh > (res_rollup[tid]["LL"]/res_rollup[tid]["MR"]) else ""
			note += "Not shown (%s-result biased); "%rank if major_ranks[rank] > major_ranks[tg_rank] else ""

			# additional fileds for full mode
			add_field = "\t%.4f\t%.4f\t%.4f\t%.4f\t%.4f\t%s\t%s\t%.2f\t%.4f\t%.4f\t%.4f\t%.4f\t%s" % (
				res_rollup[tid]["LL"]/res_rollup[tid]["TS"],                                              # LINEAR_COV
				res_rollup[tid]["LL"]/res_rollup[tid]["SL"],                                              # LINEAR_COV_MAPPED_SIG
				res_rollup[tid]["bLC"],                                                                   # BEST_LINEAR_COV
				res_rollup[tid]["MB"]/res_rollup[tid]["TS"],                                              # DOC
				res_rollup[tid]["bDOC"],                                                                  # BEST_DOC
				res_rollup[tid]["SL"],                                                                    # MAPPED_SIG_LENGTH
				res_rollup[tid]["TS"],                                                                    # TOL_SIG_LENGTH
				res_rollup[tid]["ABU"],                                                                   # ABUNDANCE
				res_rollup[tid]["RD"]/tol_abu["ROLLUP_DOC"] if tol_abu["ROLLUP_DOC"] else 0,              # REL_ABU_ROLLUP_DOC
				res_rollup[tid]["MR"]/tol_abu["READ_COUNT"] if tol_abu["READ_COUNT"] else 0,              # REL_ABU_READ_COUNT
				res_rollup[tid]["MB"]/tol_abu["TOTAL_BP_MAPPED"] if tol_abu["TOTAL_BP_MAPPED"] else 0,    # REL_ABU_TOL_BP_MAPPED
				res_rollup[tid]["LL"]/res_rollup[tid]["MR"],                                              # MLRL
				note,                                                                                     # NOTE
				#res_rollup[tid]["ML"]
			) if mode == "full" else ""

			if note and mode=="summary": continue

			#relative abundance
			o.write( "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%.4f\t%.4f\t%.4f%s\n" %
				(   rank,
					gt.taxid2name(tid),
					tid,
					res_rollup[tid]["MR"],
					res_rollup[tid]["MB"],
					res_rollup[tid]["NM"],
					res_rollup[tid]["LL"],
					res_rollup[tid]["MB"]/res_rollup[tid]["LL"],
					res_rollup[tid]["RD"],
					res_rollup[tid]["ABU"]/tol_abu["ABU"] if tol_abu["ABU"] else 0,
					add_field
				)
			)

def outputResultsAsLineage( res_rollup, o, tg_rank, mode, mc, mr, ml, mh ):
	for tid in res_rollup:
		rank = gt.taxid2rank(tid)

		if rank != tg_rank or ( mh > res_rollup[tid]["LL"]/res_rollup[tid]["MR"] or mc > res_rollup[tid]["LL"]/res_rollup[tid]["SL"] or mr > int(res_rollup[tid]["MR"]) or ml > int(res_rollup[tid]["LL"]) ):
			continue

		o.write( "%s\t%s\n" %
			( res_rollup[tid]["ABU"],
			'\t'.join( gt.taxid2lineage(tid).split('|') )
			)
		)

def readMapping( reads, db, threads, mm_penalty, presetx, samfile, logfile ):
	"""
	mapping reads to database
	"""
	input_file = " ".join(reads)

	sr_opts = "-k24 -A1 -B%s -O30 -E30 -a -n1 -m24 -s30"%mm_penalty

	bash_cmd   = "set -o pipefail; set -x;"
	mm2_cmd    = "minimap2 -x %s %s --second=yes -t%s %s.mmi %s" % (presetx, sr_opts, threads, db, input_file )
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
	message = "[%s] %s\n" % (timeSpend(start), msg)

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
	if sys.version_info < (3,0):
		sys.exit("[ERROR] Python 3.0 or above is required.")

	if not dependency_check("minimap2"):
		sys.exit("[ERROR] Executable minimap2 not found.")

	if not dependency_check("gawk"):
		sys.exit("[ERROR] Executable gawk not found.")

	#prepare output object
	out_fp = sys.stdout
	outfile = "STDOUT"
	argvs.relAbu = argvs.relAbu.upper()

	# remove previous log file if exists
	if os.path.isfile(logfile):
		os.remove(logfile)

	if not argvs.stdout:
		#create output directory if not exists
		if not os.path.exists(argvs.outdir):
			os.makedirs(argvs.outdir)
		ext = "fastq" if argvs.mode == "extract" else "tsv"
		tg_taxid = argvs.taxonomy if argvs.taxonomy else ""
		outfile = "%s/%s.%s%s.%s" % ( argvs.outdir, argvs.prefix, argvs.mode, tg_taxid, ext)
		out_fp = open( outfile, 'w')

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
	print_message( "    Mode             : %s" % argvs.mode,      argvs.silent, begin_t, logfile )
	print_message( "    Specific taxid   : %s" % argvs.taxonomy,  argvs.silent, begin_t, logfile )
	print_message( "    Threads          : %d" % argvs.threads,   argvs.silent, begin_t, logfile )
	print_message( "    Minimal L_DOC    : %s" % argvs.minCov,    argvs.silent, begin_t, logfile )
	print_message( "    Minimal L_LEN    : %s" % argvs.minLen,    argvs.silent, begin_t, logfile )
	print_message( "    Minimal reads    : %s" % argvs.minReads,  argvs.silent, begin_t, logfile )
	print_message( "    Minimal MLHL     : %s" % argvs.minMLRL,  argvs.silent, begin_t, logfile )
	print_message( "    Minimap2 path    : %s" % dependency_check("minimap2"),       argvs.silent, begin_t, logfile )

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
		db_stats = loadDatabaseStats( argvs.database + ".stats" )
	else:
		sys.exit( "[%s] ERROR: %s not found.\n" % (timeSpend(begin_t), argvs.database+".stats") )
	print_message( "Done.", argvs.silent, begin_t, logfile )

	if argvs.input:
		print_message( "Running read-mapping...", argvs.silent, begin_t, logfile )
		exitcode, cmd, msg = readMapping( argvs.input, argvs.database, argvs.threads, argvs.mismatch, argvs.presetx, samfile, logfile )
		print_message( "Logfile saved to %s." % logfile, argvs.silent, begin_t, logfile )
		#print_message( "COMMAND: %s" % cmd, argvs.silent, begin_t, logfile )

		if exitcode != 0:
			sys.exit( "[%s] ERROR: error occurred while running read mapping (exit: %s, message: %s).\n" % (timeSpend(begin_t), exitcode, msg) )
		else:
			print_message( "Done mapping reads to %s signature database." % argvs.dbLevel, argvs.silent, begin_t, logfile )
			print_message( "Mapped SAM file saved to %s." % samfile, argvs.silent, begin_t, logfile )
			sam_fp = open( samfile, "r" )

	if argvs.mode == 'class':
		processSAMfileReadClass( sam_fp, out_fp, argvs.dbLevel, argvs.taxonomy )
		print_message( "Done classifying reads. Results printed to %s." % outfile, argvs.silent, begin_t, logfile )

	elif argvs.mode == 'extract':
		processSAMfileReadExtract( os.path.abspath(samfile), out_fp, argvs.taxonomy, argvs.threads )
		print_message( "Done extracting reads to %s." % outfile, argvs.silent, begin_t, logfile )

	else:
		(res, mapped_r_cnt) = processSAMfile( os.path.abspath(samfile), argvs.threads, lines_per_process)
		print_message( "Done processing SAM file. %s reads mapped." % mapped_r_cnt, argvs.silent, begin_t, logfile )

		(res_rollup, res_tree) = taxonomyRollUp( res, db_stats, argvs.relAbu, argvs.minCov, argvs.minReads, argvs.minLen, argvs.minMLRL )
		print_message( "Done taxonomy rolling up.", argvs.silent, begin_t, logfile )

		if argvs.mode == 'summary' or argvs.mode == 'full':
			outputResultsAsRanks( res_rollup, out_fp, argvs.dbLevel, argvs.mode, argvs.minCov, argvs.minReads, argvs.minLen, argvs.minMLRL )
		#elif argvs.mode == 'tree':
		#	outputResultsAsTree( "1", res_tree, res_rollup, "", argvs.dbLevel, 0, argvs.taxonomy, out_fp, argvs.minCov, argvs.minReads, argvs.minLen )
		elif argvs.mode == 'lineage':
			outputResultsAsLineage( res_rollup, out_fp, argvs.dbLevel, argvs.mode, argvs.minCov, argvs.minReads, argvs.minLen, argvs.minMLRL )

		print_message( "Done taxonomy profiling; %s results printed to %s." % (argvs.mode, outfile), argvs.silent, begin_t, logfile )
