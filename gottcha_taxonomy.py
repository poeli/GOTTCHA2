#!/usr/bin/env python

# Po-E (Paul) Li
# B-11, Los Alamos National Lab
# Date: 05/15/2016
# Last update: 05/15/2016

import sys
import json

# functions:
#    taxid2rank
#    taxid2name
#    taxid2nameOnRank
#    taxid2taxidOnRank
#    taxid2lineage
#    taxid2lineageJSON
#    taxid2fullLineage
#    taxidIsLeaf
#    getTaxRank
#    getTaxName
#    getTaxDepth
#    getTaxParent
#    loadTaxonomy
# Other function:
#    _autoVivification()

####################
# Global variables #
####################

libPath = "/panfs/biopan01/scratch-218817/opt/src/KronaTools-2.6.1/lib"
taxonomyDir = libPath + "/../taxonomy"
DEBUG=0

taxDepths  = {}
taxParents = {}
taxRanks   = {}
taxNames   = {}
taxLeaves  = {}

####################
#      Methods     #
####################

def taxid2rank( taxID ):
	_checkTaxonomy();
	return getTaxRank(taxID)

def taxid2name( taxID ):
	_checkTaxonomy()
	return getTaxName(taxID)

def taxid2nameOnRank( taxID, r ):
	_checkTaxonomy()

	if taxID == 1: return "root"
	if r == "root": return "root"

	rank = getTaxRank(taxID)
	name = getTaxName(taxID)

	if r.startswith("strain"): return name

	while taxID:
		if rank.upper() == r.upper(): return name
		if name == 'root': break
		taxID = getTaxParent(taxID)
		rank = getTaxRank(taxID)
		name = getTaxName(taxID)

	return ""

def taxid2taxidOnRank( taxID, r ):
	_checkTaxonomy()
	rank = getTaxRank(taxID)
	name = getTaxName(taxID)

	if r == rank or ( r == 'strain' and rank == 'no rank'): return taxID
	if r == "root": return 1

	while taxID:
		if rank.upper() == r.upper(): return taxID
		if name == 'root': break

		taxID = getTaxParent(taxID)
		rank = getTaxRank(taxID)
		name = getTaxName(taxID)

	return ""

def taxidIsLeaf( taxID ):
	if taxID in taxLeaves:
		return True
	else:
		return False

def taxid2fullLineage( taxID ):
	_checkTaxonomy()
	fullLineage = ""

	while int(taxID) > 1:
		rank = getTaxRank(taxID)
		name = getTaxName(taxID)
		if not name: break
		fullLineage += "%s|%s|%s|"%(rank,taxID,name)
		taxID = taxParents[taxID]

	return fullLineage

def taxid2fullLinkDict( taxID ):
	_checkTaxonomy()
	fullLineage = ""
	link = {}

	while int(taxID) > 1:
		rank = getTaxRank(taxID)
		name = getTaxName(taxID)
		if not name: break

		parID = taxParents[taxID]
		link[parID] = taxID
		taxID = parID

	return link

def taxid2lineage( tid, print_all_rank=1, print_strain=0, replace_space2underscore=1, output_type="auto"):
	return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type)

def taxid2lineageDICT( tid, print_all_rank=1, print_strain=0, replace_space2underscore=0, output_type="DICT" ):
	return _taxid2lineage( tid, print_all_rank, print_strain, replace_space2underscore, output_type )

def _taxid2lineage(tid, print_all_rank, print_strain, replace_space2underscore, output_type):
	_checkTaxonomy()

	info = _autoVivification()
	lineage = []
	taxID = tid

	major_level = {
		'superkingdom' : 'k',
		'phylum'       : 'p',
		'class'        : 'c',
		'order'        : 'o',
		'family'       : 'f',
		'genus'        : 'g',
		'species'      : 's',
		'k'            : 'superkingdom',
		'p'            : 'phylum',
		'c'            : 'class',
		'o'            : 'order',
		'f'            : 'family',
		'g'            : 'genus',
		's'            : 'species'
	}

	level = {
		'k' : '',
		'p' : '',
		'c' : '',
		'o' : '',
		'f' : '',
		'g' : '',
		's' : ''
	}

	rank = getTaxRank(taxID)
	name = getTaxName(taxID)
	str_name = name
	if replace_space2underscore: str_name.replace(" ", "_")

	while taxID:
		if rank in major_level:
			if replace_space2underscore: name.replace(" ", "_")
			level[major_level[rank]] = name

			#for output JSON
			info[rank]["name"] = name
			info[rank]["taxid"] = taxID

		taxID = getTaxParent(taxID)
		rank = getTaxRank(taxID)
		name = getTaxName(taxID)

		if name == 'root': break

	# try to get the closest "no_rank" taxa to "type" representing subtype/group (mainly for virus)
	typeTID = getTaxType(tid)
	if typeTID:
		info["type"]["name"]  = getTaxName(typeTID)
		info["type"]["taxid"] = typeTID

	last = str_name
	for lvl in ['s','g','f','o','c','p','k']:
		if print_all_rank == 0:
			if not level[lvl]: continue

		if not level[lvl]:
			level[lvl] = "%s - no_%s_rank"%(last,lvl)
			info[major_level[lvl]]["name"]  = "%s - no_%s_rank"%(last,lvl)
			info[major_level[lvl]]["taxid"] = 0

		last=level[lvl]
		lineage.append( "%s__%s"%(lvl, level[lvl]) )

	lineage.reverse()

	if print_strain:
		lineage.append( "n__%s"%(str_name) )
		info["strain"]["name"]  = str_name
		info["strain"]["taxid"] = tid

	if output_type == "DICT":
		return info
	else:
		return "|".join(lineage)

def loadStrainName():
	return ""

def getTaxDepth( taxID ):
	return taxDepths[taxID]

def getTaxName( taxID ):
	return taxNames[taxID]

def getTaxParent( taxID ):
	taxID = taxParents[taxID]
	while int(taxID) > 1 and taxRanks[taxID] == 'no rank':
		taxID = taxParents[taxID]

	return taxID

def getTaxType( taxID ):
	origID = taxID
	lastID = taxID
	taxID = taxParents[taxID]

	while int(taxID) > 1 and taxRanks[taxID] != 'species':
		lastID = taxID
		taxID = taxParents[taxID]

	if taxRanks[taxID] != 'species':
		taxID = 0
	else:
		taxID = lastID
		if taxID == origID: taxID = 0

	return taxID

def getTaxRank( taxID ):
	return taxRanks[taxID]

def loadTaxonomy( custom_taxonomy_file="", taxonomy_file = taxonomyDir+"/taxonomy.tab" ):
	if DEBUG: sys.stderr.write( "[INFO] Open taxonomy file: %s\n"%(taxonomy_file) )

	with open(taxonomy_file, 'r') as f:
		for line in f:
			tid, depth, parent, rank, name = line.rstrip('\r\n').split('\t')
			taxParents[tid] = parent
			taxDepths[tid] = depth
			taxRanks[tid] = rank
			taxNames[tid] = name
			taxLeaves[tid] = 1;
			if parent in taxLeaves: del taxLeaves[parent]
		f.close()

	if custom_taxonomy_file:
		with open(custom_taxonomy_file, 'r') as f:
			for line in f:
				tid, name = line.rstrip('\r\n').split('\t')
				if "." in t:
					parent, sid = tid.split('.')
					taxParents[tid] = parent
					taxRanks[tid] = "strain"
					taxNames[tid] = name
					taxLeaves[tid] = 1;
		f.close()

	if DEBUG: sys.stderr.write( "[INFO] Done parsing taxonomy.tab (%d taxons loaded)\n" % len(taxParents) )

	if taxParents["2"] == "1":
		_die( "Local taxonomy database is out of date. Update using updateTaxonomy.sh." )

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

def _checkTaxonomy():
	if not len(taxParents):
		_die("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.\n")

if __name__ == '__main__':
	loadTaxonomy()
	taxid = input("Enter TaxID => ")
	while taxid:
		taxid.rstrip('\r\n')
		print( "Name: %s" % taxid2name(taxid) )
		print( "Rank: %s" % taxid2rank(taxid) )
		print( "Is leaf: %s" % taxidIsLeaf(taxid) )
		print( "Gneus: %s" % taxid2nameOnRank(taxid, "genus") )
		print( "Gneus TaxID: %s" % taxid2taxidOnRank(taxid, "genus") )
		print( "Lineage: %s\n" % taxid2lineage(taxid) )
		print( "LineageJSON: %s\n" % taxid2lineageJSON(taxid) )
		print( "Full lineage: %s\n" % taxid2fullLineage(taxid) )
		print( taxid2fullLinkDict(taxid) )
		taxid = input("Enter TaxID => ")
