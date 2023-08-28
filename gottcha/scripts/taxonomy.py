#!/usr/bin/env python

# Po-E (Paul) Li
# B-11, Los Alamos National Lab
# Date: 05/15/2016

import sys
import io
import os.path
import json
import gzip
import subprocess
import fileinput
import tarfile
import requests

####################
# Global variables #
####################

lib_path = os.path.dirname(os.path.realpath(__file__))
taxonomy_dir = lib_path + "/taxonomy_db"
DEBUG = True

taxDepths      = {}
taxParents     = {}
taxRanks       = {}
taxNames       = {}
taxMerged      = {}
taxNumChilds   = {}
accTid         = {}
tidLineage     = {}
tidLineageDict = {}

major_level_to_abbr = {
	'superkingdom' : 'k',
	'phylum'       : 'p',
	'class'        : 'c',
	'order'        : 'o',
	'family'       : 'f',
	'genus'        : 'g',
	'species'      : 's'
}
abbr_to_major_level = {
	'k'            : 'superkingdom',
	'p'            : 'phylum',
	'c'            : 'class',
	'o'            : 'order',
	'f'            : 'family',
	'g'            : 'genus',
	's'            : 'species'
}

####################
#      Methods     #
####################

def taxidStatus( taxID ):
	if taxID in taxMerged:
		return taxMerged[taxID]

	if taxID in taxNames and taxID in taxNames and taxID in taxRanks:
		if '.' in taxID:
			return "valid custom"
		return "valid"
	else:
		return "invalid"

def acc2taxid( acc ):
	_checkTaxonomy()

	accession2taxid_file = f"{taxonomy_dir}/accession2taxid.tsv"
	#remove version number#
	acc = acc.split('.')[0]

	if DEBUG: sys.stderr.write( f"[INFO] acc2taxid from file: {accession2taxid_file}\n" )

	if not acc in accTid:
		with open( accession2taxid_file ) as f:
			f.seek(0, 2)
			start = 0
			end = f.tell()
			accCur = ""
			
			if DEBUG: sys.stderr.write( f"[INFO] acc2taxid from file: {accession2taxid_file}\n" )
			
			while( acc != accCur and start < end ):
				
				posNew = (end+start)/2
				
				f.seek( posNew )
		
				if posNew != start: f.readline()

				line = f.readline()	
				
				if DEBUG: sys.stderr.write( "[INFO] start: %15d, posNew: %15d, end: %15d, line: %s" % (start, posNew, end, line) )
				if line :
					(accNew, tid) = line.split('\t')
				else:
					break

				if acc > accNew and accCur != accNew and accNew:
					if accNew: posNew = f.tell()
					start = posNew
					if start >= end: end = start+1
				else:
					end = posNew
				
				accCur = accNew

			f.close()

			if accCur == acc:
				accTid[acc] = tid.strip()
			else:
				accTid[acc] = ""

	tid = _checkTaxonomy(accTid[acc])

	return tid

def taxid2rank( taxID, guess_strain=True ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"

	if taxID == '1':
		return "root"

	if taxRanks[taxID] == "no rank" and guess_strain:
		# a leaf taxonomy is a strain
		if taxidIsLeaf(taxID):
			return "strain"
		# if not
		else:
			nmtid = taxid2nearestMajorTaxid(taxID)
			nmrank = _getTaxRank(nmtid)
			if nmrank == "species":
				return "species - others"
			else:
				return "others"
	
	return taxRanks[taxID]

def taxid2name( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown":
		return "unknown"
	else:
		return _getTaxName(taxID)

def taxid2depth( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown":
		return "unknown"
	else:
		return _getTaxDepth(taxID)

def taxid2type( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"

	origID = taxID
	lastID = taxID
	taxID = taxParents[taxID]

	while taxID != '1' and taxRanks[taxID] != 'species':
		lastID = taxID
		taxID = taxParents[taxID]

	if taxRanks[taxID] != 'species':
		taxID = 0
	else:
		taxID = lastID
		if taxID == origID: taxID = 0

	return taxID

def taxid2parent( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"

	taxID = taxParents[taxID]
	while taxID != '1' and taxRanks[taxID] == 'no rank':
		taxID = taxParents[taxID]

	return taxID

def taxid2nameOnRank( taxID, target_rank=None ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"

	if taxID == 1: return "root"
	if target_rank == "root": return "root"

	rank = _getTaxRank(taxID)
	name = _getTaxName(taxID)

	if target_rank == "strain" and taxidIsLeaf(taxID):
		return name

	while taxID:
		if rank.upper() == target_rank.upper(): return name
		if name == 'root': break
		taxID = _getTaxParent(taxID)
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)

	return ""

def taxid2taxidOnRank( taxID, target_rank=None ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"

	rank = _getTaxRank(taxID)
	name = _getTaxName(taxID)

	if target_rank == rank or ( target_rank == 'strain' and rank == 'no rank'): return taxID
	if target_rank == "root": return 1

	while taxID:
		if rank.upper() == target_rank.upper(): return taxID
		if name == 'root': break

		taxID = _getTaxParent(taxID)
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)

	return ""

def taxidIsLeaf( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return False
	if not taxID in taxNumChilds:
		return True
	else:
		return False

def taxid2fullLineage( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"
	fullLineage = ""

	while taxID != '1':
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)
		if not name: break
		fullLineage += "%s|%s|%s|"%(rank,taxID,name)
		taxID = taxParents[taxID]

	return fullLineage

def taxid2fullLinkDict( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"
	link = {}

	while taxID != '1':
		name = _getTaxName(taxID)
		if not name: break

		parID = taxParents[taxID]
		link[parID] = taxID
		taxID = parID

	return link

def taxid2nearestMajorTaxid( taxID ):
	taxID = _checkTaxonomy( taxID )
	if taxID == "unknown": return "unknown"
	ptid = _getTaxParent( taxID )
	while ptid != '1':
		tmp = taxid2rank( ptid )
		if tmp in major_level_to_abbr:
			return ptid
		else:
			ptid = _getTaxParent( ptid )

	return "1"

def taxid2lineage( tid, print_all_rank=True, print_strain=False, replace_space2underscore=True, output_type="auto"):
	return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)

def taxid2lineageDICT( tid, print_all_rank=True, print_strain=False, replace_space2underscore=False, output_type="DICT" ):
	return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type )

def taxid2lineageTEXT( tid, print_all_rank=True, print_strain=False, replace_space2underscore=True, output_type="DICT"):
	lineage = _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)
	texts = []
	for rank in major_level_to_abbr:
		if rank in lineage:
			texts.append( f"{major_level_to_abbr[rank]}__{lineage[rank]['name']}" ) 
	
	return ";".join(texts).replace(" ","_")

def _taxid2lineage(tid, print_all_rank, print_strain, replace_space2underscore, output_type):
	taxID = _checkTaxonomy( tid )
	if taxID == "unknown": return "unknown"

	if output_type == "DICT":
		if taxID in tidLineageDict: return tidLineageDict[taxID]
	else:
		if taxID in tidLineage: return tidLineage[taxID]

	info = _autoVivification()
	lineage = []

	level = {
		'k' : '',
		'p' : '',
		'c' : '',
		'o' : '',
		'f' : '',
		'g' : '',
		's' : ''
	}

	rank = taxid2rank(taxID)
	orig_rank = rank
	name = _getTaxName(taxID)
	str_name = name
	if replace_space2underscore: str_name.replace(" ", "_")

	while taxID:
		if rank in major_level_to_abbr:
			if replace_space2underscore: name.replace(" ", "_")
			level[major_level_to_abbr[rank]] = name

			#for output JSON
			info[rank]["name"] = name
			info[rank]["taxid"] = taxID

		taxID = _getTaxParent(taxID)
		rank = _getTaxRank(taxID)
		name = _getTaxName(taxID)

		if name == 'root': break

	# try to get the closest "no_rank" taxa to "type" representing subtype/group (mainly for virus)
	typeTID = taxid2type(tid)
	if typeTID:
		info["type"]["name"]  = _getTaxName(typeTID)
		info["type"]["taxid"] = typeTID

	last = str_name

	ranks = ['s','g','f','o','c','p','k']
	idx = 0
	
	# input taxid is a major rank
	if orig_rank in major_level_to_abbr:
		idx = ranks.index( major_level_to_abbr[orig_rank] )
	# if not, find the next major rank
	else:
		nmtid = taxid2nearestMajorTaxid( tid )
		nmrank = taxid2rank( nmtid )
		if nmrank == "root":
			idx = 7
		else:
			idx = ranks.index( major_level_to_abbr[nmrank] )

	for lvl in ranks[idx:]:
		if print_all_rank == 0:
			if not level[lvl]: continue

		if not level[lvl]:
			level[lvl] = "%s - no_%s_rank"%(last,lvl)
			info[abbr_to_major_level[lvl]]["name"]  = "%s - no_%s_rank"%(last,lvl)
			info[abbr_to_major_level[lvl]]["taxid"] = 0

		last=level[lvl]
		#lineage.append( "%s__%s"%(lvl, level[lvl]) )
		lineage.append( level[lvl] )

	lineage.reverse()

	if print_strain:
		if orig_rank == "strain":
			#lineage.append( "n__%s"%(str_name) )
			lineage.append( "%s"%(str_name) )
			info["strain"]["name"]  = str_name
			info["strain"]["taxid"] = tid

	if output_type == "DICT":
		tidLineageDict[tid] = info
		return info
	else:
		tidLineage[tid] = "|".join(lineage)
		return "|".join(lineage)

def _getTaxDepth( taxID ):
	return taxDepths[taxID]

def _getTaxName( taxID ):
	return taxNames[taxID]

def _getTaxParent( taxID ):
	return taxParents[taxID]

def _getTaxRank( taxID ):
	return taxRanks[taxID]

def lca_taxid(taxids):
	""" lca_taxid
	Return lowest common ancestor (LCA) taxid of input taxids
	"""
	ranks = ['strain','species','genus','family','order','class','phylum','superkingdom']

	merged_dict = _autoVivification()
	for tid in taxids:
		lineage = taxid2lineageDICT(tid, 1, 1)
		for r in ranks:
			if not r in lineage:
				ttid = "0"
			else:
				ttid = lineage[r]['taxid']

			if ttid in merged_dict[r]:
				merged_dict[r][ttid] += 1
			else:
				merged_dict[r][ttid] = 1

	for r in ranks:
		if len(merged_dict[r]) == 1:
			for ttid in merged_dict[r]:
				# skip if no tid in this rank
				if ttid=="0":
					continue
				return ttid

	return '1'

def loadTaxonomy( dbpath=None, cus_taxonomy_file=None, auto_download=True, debug=False ):
	global taxonomy_dir
	global DEBUG

	DEBUG = debug

	if dbpath:
		taxonomy_dir = dbpath

	if DEBUG: sys.stderr.write( f"[INFO] Taxonomy directory: {taxonomy_dir}\n" )

	#NCBI ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz
	taxdump_tgz_file = taxonomy_dir+"/taxdump.tar.gz"

	#raw taxonomy dmp files from NCBI
	names_dmp_file = taxonomy_dir+"/names.dmp"
	nodes_dmp_file = taxonomy_dir+"/nodes.dmp"
	merged_dmp_file = taxonomy_dir+"/merged.dmp"

	#parsed taxonomy tsv file
	taxonomy_file = taxonomy_dir+"/taxonomy.tsv"
	merged_taxonomy_file = taxonomy_dir+"/taxonomy.merged.tsv"

	#custom taxonomy file
	if not cus_taxonomy_file:
		cus_taxonomy_file = taxonomy_dir+"/taxonomy.custom.tsv"

	# checking if taxonomy files downloaded
	if not os.path.isfile( taxdump_tgz_file ) \
	  and not os.path.isfile( merged_taxonomy_file ) \
	  and not (os.path.isfile( names_dmp_file ) and os.path.isfile( nodes_dmp_file )) \
	  and not os.path.isfile( cus_taxonomy_file ):
		if DEBUG: sys.stderr.write( "[INFO] Local taxonomy files not found.\n" )
		if auto_download:
			url = 'http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
			if DEBUG: sys.stderr.write( f"[INFO] Auto downloading taxanomy from {url}...\n" )
			# download taxonomy file if auto_download enabled
			r = requests.get(url)
			if not os.path.exists( taxonomy_dir ):
				os.makedirs( taxonomy_dir )
			with open(f'{taxonomy_dir}/taxdump.tar.gz', 'wb') as f:
				f.write(r.content)
			if os.path.getsize( f'{taxonomy_dir}/taxdump.tar.gz' ):
				if DEBUG: sys.stderr.write( f"[INFO] Saved to {taxonomy_dir}/taxdump.tar.gz.\n" )
			else:
				_die( "[ERROR] Failed to download taxonomy files.\n" )	
		else:
			if DEBUG: sys.stderr.write( f"[INFO] Auto-download is off.\n" )
			_die( "[ERROR] No available taxonomy files.\n" )

	# try to load taxonomy from taxonomy.tsv
	if os.path.isfile( taxdump_tgz_file ):
		try:
			if DEBUG: sys.stderr.write( f"[INFO] Open taxonomy file: {taxdump_tgz_file}\n" )
			tar = tarfile.open(taxdump_tgz_file, "r:gz")
			
			# read name from names.dmp
			if DEBUG: sys.stderr.write( "[INFO] Extract taxonomy names file: names.dmp\n" )
			member = tar.getmember("names.dmp")
			f = tar.extractfile(member)
			for line in f.readlines():
				tid, name, tmp, nametype = line.decode('utf8').rstrip('\r\n').split('\t|\t')
				if not nametype.startswith("scientific name"):
					continue
				taxNames[tid] = name
			f.close()
			
			# read taxonomy info from nodes.dmp
			if DEBUG: sys.stderr.write( "[INFO] Extract taxonomy nodes file: nodes.dmp\n" )
			member = tar.getmember("nodes.dmp")
			f = tar.extractfile(member)
			for line in f.readlines():
				fields = line.decode('utf8').rstrip('\r\n').split('\t|\t')
				tid = fields[0]
				parent = fields[1]
				taxParents[tid] = parent
				taxDepths[tid] = taxDepths[parent]+1 if parent in taxDepths else 0 # could have potiential bug if child node is parsed before parent node.
				taxRanks[tid] = fields[2]
				if parent in taxNumChilds:
					taxNumChilds[parent] += 1
				else:
					taxNumChilds[parent] = 1
			f.close()

			if DEBUG: sys.stderr.write( "[INFO] Extract merged taxonomy node file: merged.dmp\n")
			member = tar.getmember("merged.dmp")
			f = tar.extractfile(member)
			for line in f.readlines():
				fields = line.decode('utf8').rstrip('\r\n').replace("\t","").split('|')
				mtid = fields[0]
				tid = fields[1]
				taxMerged[mtid] = tid
			f.close()

		except IOError:
			_die( "Failed to load taxonomy from %s\n"%taxdump_tgz_file )
	elif os.path.isfile( names_dmp_file ):
		try:
			# read name from names.dmp
			if DEBUG: sys.stderr.write( f"[INFO] Open taxonomy name file: {names_dmp_file}\n"  )
			with open(names_dmp_file) as f:
				for line in f:
					tid, name, tmp, nametype = line.rstrip('\r\n').split('\t|\t')
					if not nametype.startswith("scientific name"):
						continue
					taxNames[tid] = name
				f.close()
				if DEBUG: sys.stderr.write( f"[INFO] Done parsing taxonomy name file.\n" )	

			# read taxonomy info from nodes.dmp
			if DEBUG: sys.stderr.write( f"[INFO] Open taxonomy node file: {nodes_dmp_file}\n" )
			with open(nodes_dmp_file) as f:
				for line in f:
					fields = line.rstrip('\r\n').split('\t|\t')
					tid = fields[0]
					parent = fields[1]
					taxParents[tid] = parent
					taxDepths[tid] = taxDepths[parent]+1 if parent in taxDepths else 0 # could have potiential bug if child node is parsed before parent node.
					taxRanks[tid] = fields[2]
					if parent in taxNumChilds:
						taxNumChilds[parent] += 1
					else:
						taxNumChilds[parent] = 1
				f.close()
				if DEBUG: sys.stderr.write( f"[INFO] Done parsing taxonomy node file.\n" )

			if os.path.isfile( merged_dmp_file ):
				if DEBUG: sys.stderr.write( f"[INFO] Open merged taxonomy node file: {merged_dmp_file}\n" )
				with open(merged_dmp_file) as f:
					for line in f:
						line = line.rstrip('\r\n')
						if not line: continue
						fields = line.replace("\t","").split('|')
						mtid = fields[0]
						tid = fields[1]
						taxMerged[mtid] = tid
					f.close()
					if DEBUG: sys.stderr.write( f"[INFO] Done parsing merged taxonomy file.\n" )
		except IOError:
			_die( "Failed to open taxonomy files (taxonomy.tsv, nodes.dmp and names.dmp).\n" )
	elif os.path.isfile( taxonomy_file ):
		if DEBUG: sys.stderr.write( "[INFO] Open taxonomy file: %s\n"% taxonomy_file )
		try:
			with open(taxonomy_file) as f:
				for line in f:
					tid, depth, parent, rank, name = line.rstrip('\r\n').split('\t')
					taxParents[tid] = parent
					taxDepths[tid] = depth
					taxRanks[tid] = rank
					taxNames[tid] = name
					if parent in taxNumChilds:
						taxNumChilds[parent] += 1
					else:
						taxNumChilds[parent] = 1
				f.close()
				if DEBUG: sys.stderr.write( f"[INFO] Done parsing taxonomy file.\n" )

			#try to load merged taxids
			if os.path.isfile( merged_taxonomy_file ):
				if DEBUG: sys.stderr.write( "[INFO] Open merged taxonomy node file: %s\n"% merged_taxonomy_file )
				with open(merged_taxonomy_file) as f:
					for line in f:
						line = line.rstrip('\r\n')
						if not line: continue
						mtid, tid = line.split('\t')
						taxMerged[mtid] = tid
					f.close()
					if DEBUG: sys.stderr.write( f"[INFO] Done parsing merged taxonomy file.\n" )
		except IOError:
			_die( "Failed to open taxonomy file: %s.\n" % taxonomy_file )

	# try to load custom taxonomy from taxonomy.custom.tsv
	if os.path.isfile( cus_taxonomy_file ):
		if DEBUG: sys.stderr.write( "[INFO] Open custom taxonomy node file: %s\n"% cus_taxonomy_file)
		try:
			with open(cus_taxonomy_file) as f:
				for line in f:
					line = line.rstrip('\r\n')
					if not line: continue
					tid, depth, parent, rank, name = line.split('\t')
					taxParents[tid] = parent
					taxDepths[tid] = depth
					taxRanks[tid] = rank
					taxNames[tid] = name
					if parent in taxNumChilds:
						taxNumChilds[parent] += 1
					else:
						taxNumChilds[parent] = 1
				f.close()
				if DEBUG: sys.stderr.write( f"[INFO] Done parsing custom taxonomy file.\n" )
		except IOError:
			_die( "Failed to open custom taxonomy file: %s.\n" % cus_taxonomy_file )

	if DEBUG: sys.stderr.write( "[INFO] Done parsing taxonomy files (%d taxons loaded)\n" % len(taxParents) )

##########################
##  Internal functions  ##
##########################

class _autoVivification(dict):
	"""Implementation of perl's autovivification feature."""
	def __getitem__(self, item):
		try:
			return dict.__getitem__(self, item)
		except KeyError:
			value = self[item] = type(self)()
			return value

def _die( msg ):
	sys.exit(msg)

def _checkTaxonomy(taxID="", acc=""):
	if not len(taxParents):
		_die("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.\n")

	if taxID:
		if taxID in taxMerged:
			taxID = taxMerged[taxID]

	if taxID in taxNames and taxID in taxParents and taxID in taxParents:
		return taxID
	else:
		return "unknown"

if __name__ == '__main__':
	#loading taxonomy
	loadTaxonomy( sys.argv[1] if len(sys.argv)>1 else None, debug=True)

	inid = 0
	try:
		inid = input("\nEnter acc/taxid: ")
	except:
		inid = 0

	while(inid):
		if inid[0] in "1234567890":
			taxid = inid
		else:
			taxid = acc2taxid( inid )
			print( "acc2taxid( %s ) => %s"   % (inid, taxid) )

		if taxid:
			print( "taxid2name( %s )                 => %s" % (taxid, taxid2name(taxid)) )
			print( "taxid2rank( %s )                 => %s" % (taxid, taxid2rank(taxid)) )
			print( "taxid2type( %s )                 => %s" % (taxid, taxid2type(taxid)) )
			print( "taxid2depth( %s )                => %s" % (taxid, taxid2depth(taxid)) )
			print( "taxid2parent( %s )               => %s" % (taxid, taxid2parent(taxid)) )
			print( "taxidIsLeaf( %s )                => %s" % (taxid, taxidIsLeaf(taxid)) )
			print( "taxid2nearestMajorTaxid( %s )    => %s" % (taxid, taxid2nearestMajorTaxid(taxid)) )
			print( "taxid2nameOnRank( %s, 'genus')   => %s" % (taxid, taxid2nameOnRank(taxid, "genus")) )
			print( "taxid2taxidOnRank( %s, 'genus')  => %s" % (taxid, taxid2taxidOnRank(taxid, "genus")) )
			print( "taxid2nameOnRank( %s, 'phylum')  => %s" % (taxid, taxid2nameOnRank(taxid, "phylum")) )
			print( "taxid2taxidOnRank( %s, 'phylum') => %s" % (taxid, taxid2taxidOnRank(taxid, "phylum")) )
			print( "taxid2lineage( %s )              => %s" % (taxid, taxid2lineage(taxid)) )
			print( "taxid2lineageDICT( %s, 1, 1 )    => %s" % (taxid, taxid2lineageDICT(taxid,1,1)) )
			print( "taxid2lineageTEXT( %s, 1, 1 )    => %s" % (taxid, taxid2lineageTEXT(taxid,1,1)) )
			print( "taxid2fullLineage( %s )          => %s" % (taxid, taxid2fullLineage(taxid)) )
			print( "taxid2fullLinkDict( %s )         => %s" % (taxid, taxid2fullLinkDict(taxid)) )
		else:
			print( "No taxid found.\n" )

		try:
			inid = input("\nEnter acc/taxid: ")
		except:
			inid = 0
