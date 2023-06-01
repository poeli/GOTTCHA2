#!/usr/bin/env python

# Po-E (Paul) Li
# B-11, Los Alamos National Lab
# Date: 05/15/2016
# Latest Update: 05/05/2023

import sys
import os
import tarfile
import requests
import logging
from . import __version__
from typing import Union, Optional, Literal

logger = logging.getLogger(__name__)

# Set default path of `taxonomy_db/` and `major_level_to_abbr.json`:
# The default `taxonomy_db/` path is your the location
lib_path = os.path.dirname(os.path.realpath(__file__))
taxonomy_dir = f"{lib_path}/taxonomy_db"
abbr_json_path = f"{taxonomy_dir}/major_level_to_abbr.json"

if os.path.isdir( "./taxonomy_db" ):
    taxonomy_dir = "./taxonomy_db"
    # if there is a new `major_level_to_abbr.json` in in the taxonomy_dir, use the json file.
    # Otherwise, use the default one comes with this package
    if os.path.isfile( f"{taxonomy_dir}/major_level_to_abbr.json" ):
        abbr_json_path = f"{taxonomy_dir}/major_level_to_abbr.json"
    else:
        pass
elif os.path.isdir( os.getcwd()+"/taxonomy_db" ):
    taxonomy_dir = os.getcwd()+"/taxonomy_db"

# init global dir
taxDepths      = {}
taxParents     = {}
taxRanks       = {}
taxNames       = {}
taxMerged      = {}
taxNumChilds   = {}
accTid         = {}
tidLineage     = {}
tidLineageDict = {}
nameTid        = {}
major_level_to_abbr = {}
abbr_to_major_level = {}
df_names = None

# --- helper functions ---
def _getTaxDepth(tid: str) -> str:
    """Get the depth of a taxonomy ID [warning: only support Kraken taxa inputs]"""
    if tid in taxMerged: tid = taxMerged[tid]
    return taxDepths[tid]

def _getTaxName(tid: str) -> str:
    """Get the name of a taxonomy ID"""
    if tid in taxMerged: tid = taxMerged[tid]
    return taxNames[tid]

def _getTaxParent(tid: str) -> str:
    """Get the parent taxid of a taxonomy ID"""
    if tid in taxMerged: tid = taxMerged[tid]
    return taxParents[tid]

def _getTaxRank(tid: str) -> str:
    """Get the rank of a taxonomy ID"""
    if tid in taxMerged: tid = taxMerged[tid]
    return taxRanks[tid]

def _die(msg: str) -> str:
    sys.exit(msg)

def _checkTaxonomy(tid: Union[int, str]):
    """Check if a taxonomy ID is present in the taxonomy database"""

    if not len(taxParents):
        logger.fatal("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.")
        _die("Taxonomy not loaded. \"loadTaxonomy()\" must be called first.")

    if tid:
        # tid must be in string type
        tid = str(tid)

        # convert to merged tid first if needs
        if tid in taxMerged:
            logger.info( f"Merged tid found: {tid} -> {taxMerged[tid]}." )
            tid = taxMerged[tid]

        if (tid in taxNames) and (tid in taxParents):
            return tid
        else:
            return "unknown"

def _taxid2fullLink(tid: Union[int, str]) -> dict:
    """
    Returns a dictionary containing the full lineage of the target taxon (e.g. {`[parent_tid]`: `tid`,...}).

    Args:
        tid (Union[int, str]): Taxonomy ID of the target taxon.

    Returns:
        dict: A dictionary containing the full lineage of the target taxon.

    """
    tid = _checkTaxonomy(tid)
    if tid == "unknown": return {}
    link = _autoVivification()

    while tid != '1':
        name = _getTaxName(tid)
        if not name: break
        
        parID = _getTaxParent(tid)
        link[parID] = tid
        tid = parID
    
    return link


def _taxid2lineage(tid: Union[int, str], all_major_rank: bool=True, print_strain: bool=True, 
                   space2underscore: bool=False, guess_type: bool=False):
    """Get the lineage of a taxonomy ID as a dictionary"""

    tid = _checkTaxonomy( tid )
    if tid == "unknown": return {}
    # if tid in tidLineageDict: return tidLineageDict[tid]

    info = _autoVivification()
    level = {abbr: '' for abbr in abbr_to_major_level}

    rank = _getTaxRank(tid)
    orig_rank = rank
    name = _getTaxName(tid)
    str_name = name
    if space2underscore: str_name = str_name.replace(" ", "_")

    while tid:
        if rank in major_level_to_abbr:
            if space2underscore: name = name.replace(" ", "_")
            level[major_level_to_abbr[rank]] = name

            #for output JSON
            info[rank]["name"] = name
            info[rank]["taxid"] = tid

        tid = _getTaxParent(tid)
        rank = _getTaxRank(tid)
        name = _getTaxName(tid)

        if name == 'root': break

    # try to get the closest "no_rank" taxa to "type" representing subtype/group (mainly for virus)
    if guess_type==True:
        typeTID = taxid2type(tid)
        if typeTID:
            info["type"]["name"]  = _getTaxName(typeTID)
            info["type"]["taxid"] = typeTID

    last = str_name

    ranks = list(abbr_to_major_level.keys())
    ranks.reverse()
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
        if all_major_rank == False:
            if not level[lvl]: continue

        if not level[lvl]:
            level[lvl] = "%s - no_%s_rank"%(last,lvl)
            info[abbr_to_major_level[lvl]]["name"]  = "%s - no_%s_rank"%(last,lvl)
            info[abbr_to_major_level[lvl]]["taxid"] = 0

        last=level[lvl]

    if print_strain==True:
        if orig_rank == "strain":
            info["strain"]["name"]  = str_name
            info["strain"]["taxid"] = tid

    tidLineageDict[tid] = info
    return info

def _loadAbbrJson(abbr_json_path: str) -> None:
    """Load abbreviations for major ranks from a JSON file"""

    import json
    global major_level_to_abbr, abbr_to_major_level

    # Opening JSON file
    with open(abbr_json_path) as f:    
        major_level_to_abbr = json.load(f)
        f.close

    if len(major_level_to_abbr):
        abbr_to_major_level = {v: k for k, v in major_level_to_abbr.items()}
    else:
        logger.fatal( f"None of the major level to aberration loaded from {abbr_json_path}." )
        _die(f"[ERROR] None of the major level to aberration loaded from {abbr_json_path}.")

class _autoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

# --- main functions ---

def taxid2rank(tid: Union[int, str], guess_strain: bool=True) -> str:
    """
    Get the taxonomic rank of a given taxonomic ID.
    
    Args:
        tid (Union[int, str]): Taxonomic ID.
        guess_strain (bool, optional): Whether to guess the strain when the rank, strain, is not available. 
            Defaults to True.
    
    Returns:
        str: The taxonomic rank.
    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown": return "unknown"

    if tid == '1':
        return "root"

    if taxRanks[tid] == "no rank" and guess_strain:
        # a leaf taxonomy is a strain
        if taxidIsLeaf(tid):
            return "strain"
        # if not
        else:
            nmtid = taxid2nearestMajorTaxid(tid)
            nmrank = _getTaxRank(nmtid)
            if nmrank == "species":
                return "species - others"
            else:
                return "others"
    
    return taxRanks[tid]

def taxid2name(tid: Union[int, str]) -> str:
    """
    Get the taxonomic name of a given taxonomic ID.
    
    Args:
        tid (Union[int, str]): Taxonomic ID.
    
    Returns:
        str: The taxonomic name.
    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown":
        return "unknown"
    else:
        return _getTaxName(tid)

def taxid2depth(tid: Union[int, str]) -> int:
    """
    Get the depth of a given taxonomic ID.
    
    Args:
        tid (Union[int, str]): Taxonomic ID.
    
    Returns:
        int: The depth of the taxonomic ID in the taxonomy tree.
    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown":
        return "unknown"
    else:
        return _getTaxDepth(tid)

def taxid2type(tid: Union[int, str]):
    """
    Guessing the `type` of a given taxonomic ID. It could be rank of `group`, `serotype`...etc.
    
    Args:
        tid (Union[int, str]): Taxonomic ID.
    
    Returns:
        Union[int, str]: The type of the taxonomic ID.
    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown": return "unknown"

    origID = tid
    lastID = tid
    tid = taxParents[tid]

    while tid != '1' and taxRanks[tid] != 'species':
        lastID = tid
        tid = taxParents[tid]

    if taxRanks[tid] != 'species':
        tid = 0
    else:
        tid = lastID
        if tid == origID: tid = 0

    return tid

def taxid2parent(tid: Union[int, str], norank: bool=False) -> str:
    """
    Get the parent of a given taxonomic ID.
    
    Args:
        tid (Union[int, str]): Taxonomic ID.
        norank (bool): Find the next parent if it's 'no-rank'.
    
    Returns:
        str: The taxonomic ID of the parent.
    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown": return "unknown"

    tid = _getTaxParent(tid)

    if not norank:
        while tid != '1' and (taxRanks[tid] == 'no rank'):
            tid = _getTaxParent(tid)

    return tid

def name2taxid(name, rank=None, superkingdom=None, fuzzy=True, cutoff=0.7, max_matches=3, reset=False, expand=True) -> list:
    """
    Get the taxonomic ID of a given taxonomic name.
    
    Args:
        name (str): Taxonomic scientific name.
        rank (str, optional): The expected rank of the taxonomic name.
        superkingdom (str, optional): The expected superkingdom of the taxonomic name.
        fuzzy (bool, optional): Whether to allow fuzzy search. Defaults to True.
        cutoff (float, optional): Similarity cutoff for difflib.get_close_matches(). 
            Only apply to `expand` mode. Cutoff will set to 1 if `fuzzy` set to False. Defaults to 0.7.
        max_matches (int, optional): Reporting max number of taxid. Defaults to 3.
        reset (bool, optional): Mapping results are cached. Whether to clean up previous searches. 
            Defaults to False.
        expand (bool, optional): Search the entire 'names.dmp' if True, otherwise search sientific names only. 
            Defaults to False.
    
    Returns:
        list: The list of matched taxonomic ID.
    """
    global nameTid, df_names, taxonomy_dir
    import pandas as pd

    if not fuzzy: cutoff=1

    # if expand is True, loading names.dmp
    names_dmp_file = taxonomy_dir+"/names.dmp"

    if df_names is None and expand and os.path.isfile( names_dmp_file ):
        df_names = pd.read_csv(names_dmp_file, 
                         sep='\t\|\t', 
                         engine='python', 
                         header=None, 
                         names=['taxid', 'name', 'annot', 'type'], 
                         usecols=['taxid','name'],
                         index_col='name')
    
    if not name in nameTid or reset:
        matched_taxid = []

        if expand:
            import difflib
            matches = difflib.get_close_matches(name, df_names.index, max_matches, cutoff)
            logger.debug(f'{name}: {matches}')
            df_temp = df_names.loc[matches,:]

            if rank:
                df_temp['rank'] = df_temp.taxid.apply(taxid2rank)
                idx = df_temp['rank']==rank
                df_temp = df_temp[idx]
            
            if superkingdom:
                df_temp['sk'] = df_temp.taxid.apply(lambda x: taxid2nameOnRank(x, 'superkingdom'))
                idx = df_temp['sk']==superkingdom
                df_temp = df_temp[idx]
            
            nameTid[name] = df_temp.head(max_matches).taxid.to_list()
            return nameTid[name]

        else: # searching scientific names only
            for taxid in taxNames:
                if fuzzy==True:
                    if not name in taxNames[taxid]:
                        continue
                else:
                    if name!=taxNames[taxid]:
                        continue
                
                if rank:
                    if _getTaxRank(taxid)==rank:
                        matched_taxid.append(taxid)
                else:
                    matched_taxid.append(taxid)

                # return when the first match found
                if len(matched_taxid)==max_matches:
                    nameTid[name] = matched_taxid
                    return nameTid[name]
            
            nameTid[name] = matched_taxid
            return nameTid[name][:max_matches]
    else:
        return nameTid[name][:max_matches]

def taxid2nameOnRank(tid: Union[int, str], target_rank=None) -> str:
    """
    Get the taxonomic name of a given taxonomic ID at a specific rank.
    
    Args:
        tid (Union[int, str]): Taxonomic ID.
        target_rank (str, optional): Target rank. Defaults to None.
    
    Returns:
        str: The taxonomic name at the target rank.
    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown": return "unknown"
    if tid == 1: return "root"
    if target_rank == "root": return "root"

    if tid:
        rank = _getTaxRank(tid)
        name = _getTaxName(tid)

        if target_rank == "strain" and taxidIsLeaf(tid):
            return name

        while tid:
            if rank.upper() == target_rank.upper(): return name
            if name == 'root': break
            tid = _getTaxParent(tid)
            rank = _getTaxRank(tid)
            name = _getTaxName(tid)
    else:
        return ""

def taxid2taxidOnRank(tid: Union[int, str], target_rank=None ) -> str:
    """
    Returns the taxonomy ID of the nearest parent taxon at the specified rank.

    Args:
        tid (Union[int, str]): Taxonomy ID of the target taxon.
        target_rank (str): The target rank to search for. Defaults to None.

    Returns:
        str: Taxonomy ID of the nearest parent taxon at the specified rank.

    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown": return "unknown"

    if tid:
        rank = _getTaxRank(tid)
        name = _getTaxName(tid)

        if target_rank == rank or ( target_rank == 'strain' and rank == 'no rank'): return tid
        if target_rank == "root": return 1

        while tid:
            if rank.upper() == target_rank.upper(): return tid
            if name == 'root': break

            tid = _getTaxParent(tid)
            rank = _getTaxRank(tid)
            name = _getTaxName(tid)
    else:
        return ""

def taxidIsLeaf(tid: Union[int, str]) -> bool:
    """
    Checks if the taxonomy ID corresponds to a leaf node in the taxonomic tree.

    Args:
        tid (Union[int, str]): Taxonomy ID of the target taxon.

    Returns:
        bool: True if the taxon is a leaf node, False otherwise.

    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown": return False
    if not tid in taxNumChilds:
        return True
    else:
        return False


def taxid2fullLineage(tid: Union[int, str], sep: str='|', use_rank_abbr=False, space2underscore=True) -> str:
    """
    Returns the full lineage of the target taxon in a specified format.

    Args:
        tid (Union[int, str]): Taxonomy ID of the target taxon.
        sep (str): Separator used to separate fields in the output. Defaults to '|'.
        use_rank_abbr (bool): If True, abbreviated rank names are used. Defaults to False.
        space2underscore (bool): If True, spaces in the output are replaced by underscores. Defaults to True.

    Returns:
        str: Full lineage of the target taxon in the specified format.

    """
    link = _taxid2fullLink(tid)
    texts = []
    if len(link):
        for p_taxID in link:
            tid = link[p_taxID]
            rank = _getTaxRank(tid)
            name = _getTaxName(tid)

            if use_rank_abbr and (rank in major_level_to_abbr):
                rank =  major_level_to_abbr[rank]

            if sep == ';':
                texts.append(f"{rank}__{name}")
            else:
                texts.append(f"{rank}|{tid}|{name}")
    
    texts.reverse()

    if space2underscore:
        return sep.join(texts).replace(' ', '_')
    else:
        return sep.join(texts)

def taxid2fullLinkDict(tid: Union[int, str]) -> str:
    """
    Returns a dictionary containing the full lineage of the target taxon.

    Args:
        tid (Union[int, str]): Taxonomy ID of the target taxon.

    Returns:
        dict: A dictionary containing the full lineage of the target taxon.

    """
    return _taxid2fullLink(tid)

def taxid2nearestMajorTaxid(tid: Union[int, str]) -> str:
    """
    Returns the taxonomy ID of the nearest parent taxon at a major rank.

    Args:
        tid (Union[int, str]): Taxonomy ID of the target taxon.

    Returns:
        str: Taxonomy ID of the nearest parent taxon at a major rank.

    """

    tid = _checkTaxonomy(tid)
    if tid == "unknown": return "unknown"
    ptid = _getTaxParent(tid)
    while ptid != '1':
        tmp = _getTaxRank(ptid)
        if tmp in major_level_to_abbr:
            return ptid
        else:
            ptid = _getTaxParent(ptid)

    return "1"

def taxid2lineage(tid: Union[int, str], all_major_rank=True, print_strain=False, space2underscore=False, sep="|") -> str:
    """
    Returns the taxonomic lineage for a given taxonomic identifier (tid) as a formatted string.

    Parameters:
        tid (Union[int, str]): A taxonomic identifier, which can be either an integer or a string.
        all_major_rank (bool): If True, all major taxonomic ranks will be included in the lineage; if False, only the lowest common ancestor will be included. Default is True.
        print_strain (bool): If True, strain information will be included in the lineage; if False, it will be omitted. Default is False.
        space2underscore (bool): If True, spaces in the taxonomic names will be replaced with underscores; if False, they will be left as spaces. Default is False.
        sep (str): The separator used to join the taxonomic ranks, taxids, and names in the returned string. Default is "|".

    Returns:
        str: A formatted string containing the taxonomic lineage information, with each rank separated by the specified separator (default is "|").

    """
    
    lineage = _taxid2lineage( tid, all_major_rank, print_strain, space2underscore)
    texts = []
    for rank in major_level_to_abbr:
        if rank in lineage:
            if print_strain==False and rank=="strain":
                continue
            if sep == ";":
                texts.append( f"{major_level_to_abbr[rank]}__{lineage[rank]['name']}" )
            else:
                texts.append( f"{rank}|{lineage[rank]['taxid']}|{lineage[rank]['name']}" ) 
    
    if space2underscore:
        return sep.join(texts).replace(' ', '_')
    else:
        return sep.join(texts) 

def taxid2lineageDICT(tid: Union[int, str], all_major_rank=True, print_strain=True, space2underscore=False, guess_type=False):
    return _taxid2lineage( tid, all_major_rank, print_strain, space2underscore, guess_type)

def lca_taxid(taxids: list) -> str:
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

def acc2taxid_raw(acc: str, accession2taxid_file: Optional[str] = None) -> str:
    """
    Get the taxonomy ID for a given accession from NCBI accession2taxid tsv file.

    Args:
        acc (str): The accession number to look up.
        accession2taxid_file (str, optional): The path to the accession to taxonomy ID mapping file. If not specified, the default file in the taxonomy directory will be used.

    Returns:
        str: The taxonomy ID for the given accession.
    """

    # Remove version number
    acc = acc.split('.')[0]

    if not acc in accTid:
        logger.info( f"acc2taxid from file: {accession2taxid_file}" )
        with open( accession2taxid_file ) as f:
            f.seek(0, 2)
            start = 0
            end = f.tell()
            accCur = ""
            
            while( acc != accCur and start < end ):
                posNew = (end+start)/2
                f.seek( posNew )
                if posNew != start: f.readline()
                line = f.readline()    
                
                logger.debug( "start: %15d, posNew: %15d, end: %15d, line: %s" % (start, posNew, end, line) )
                if line :
                    (accNew, accNewVer, tid, gi) = line.split('\t')
                else:
                    break

                logger.debug( f'[acc, accNew, accCur]=[{acc}, {accNew}, {accCur}]')
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

    return accTid[acc]

def acc2taxid(acc: str, type: Optional[str] = 'nuc') -> str:
    """
    Get the taxonomy ID for a given accession.

    Args:
        acc (str): The accession number to look up.
        type (str, optional): Type of the acession number, either nuc, prot, or pdb. Default is 'nuc'.

    Returns:
        str: The taxonomy ID for the given accession.
    """
    global taxonomy_dir
    import glob
    acc2taxid_files = []
    
    # preparing accession2taxid files
    if type == 'nuc':
        acc2taxid_files = [
            f'{taxonomy_dir}/accession2taxid/nucl_gb.accession2taxid',
            f'{taxonomy_dir}/accession2taxid/nucl_wgs.accession2taxid.EXTRA',
            f'{taxonomy_dir}/accession2taxid/nucl_wgs.accession2taxid',
            f'{taxonomy_dir}/accession2taxid/dead_nucl.accession2taxid',
            f'{taxonomy_dir}/accession2taxid/dead_wgs.accession2taxid',
        ]
    elif type == 'prot':
        acc2taxid_files = [
            f'{taxonomy_dir}/accession2taxid/prot.accession2taxid.FULL',
            f'{taxonomy_dir}/accession2taxid/dead_prot.accession2taxid',
        ]
    elif type == 'pdb':
        acc2taxid_files = [
            f'{taxonomy_dir}/accession2taxid/pdb.accession2taxid'
        ]

    # check if accession2taxid files exist
    for acc2taxid_file in acc2taxid_files:
        if not os.path.isfile(acc2taxid_file):
            acc2taxid_files.remove(acc2taxid_file)

    # download accession2taxid files if not exist
    if len(acc2taxid_files) == 0:
        logger.info( f"NCBI accession2taxid data not found." )
        NCBITaxonomyDownload(accession2taxid=True)

    for acc2taxid_file in acc2taxid_files:
        taxid = acc2taxid_raw(acc, accession2taxid_file=acc2taxid_file)
        if taxid: return taxid

    return ""

def taxid2decendentOnRank(tid: Union[int, str], target_rank=None) -> list:
    """
    Return a list of taxids for all descendants of the given taxid at the specified target rank.

    Arguments:
    - tid: The taxid to retrieve the descendants for. Can be either an integer or a string.
    - target_rank: The rank at which to retrieve the descendants. Defaults to None, which means all descendants will be returned.

    Returns:
    - A list of taxids for all descendants of the given taxid at the specified target rank.
    """
    tid = _checkTaxonomy(tid)
    decd_tids = [k for k, v in taxParents.items() if v == tid]

    tids = []

    if not target_rank:
        return decd_tids
    else:
        for tid in decd_tids:
            if target_rank == _getTaxRank(tid):
                tids.append(tid)
            else:
                if taxidIsLeaf(tid):
                    continue
                else:
                    tids.extend(taxid2decendentOnRank(tid, target_rank))

        return tids

def loadTaxonomy(dbpath: Optional[str] = None,
                 cus_taxonomy_file: Optional[str] = None, 
                 cus_taxonomy_format: str = 'tsv',
                 auto_download: bool = True) -> None:
    """
    Load taxonomy files into memory for use in subsequent conversions.

    Args:
        dbpath (str, optional): Path to store and load the taxonomy files. Defaults to None.
        cus_taxonomy_file (str, optional): Path to a custom taxonomy file. Defaults to None.
        cus_taxonomy_format (str, optional): Format of the custom taxonomy file, one of ['tsv','mgnify_lineage','gtdb_taxonomy','gtdb_metadata']. Defaults to 'tsv'.
        auto_download (bool, optional): If True, automatically download the taxonomy files if they are not found locally. Defaults to True.

    Returns:
        None
    """
    global taxonomy_dir, abbr_json_path

    if dbpath:
        taxonomy_dir = dbpath

    logger.debug( f"v{__version__}" )

    logger.debug( f"Taxonomy directory: {taxonomy_dir}" )

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path)

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

    # checking if taxonomy files provided
    if not os.path.isfile( taxdump_tgz_file ) \
        and not os.path.isfile( merged_taxonomy_file ) \
        and not os.path.isfile( taxonomy_file ) \
        and not (os.path.isfile( names_dmp_file ) and os.path.isfile( nodes_dmp_file )) \
        and not os.path.isfile( cus_taxonomy_file ):
        
        logger.info( f"No taxonomy files not found." )
        if auto_download:
            NCBITaxonomyDownload(taxonomy_dir)
        else:
            logger.info( f"Auto-download is off." )
            logger.fatal( f"No available taxonomy files." )
            _die( "[ERROR] No available taxonomy files." )

    # try to load taxonomy from taxonomy.tsv
    if os.path.isfile( nodes_dmp_file ) and  os.path.isfile( names_dmp_file ):
        loadNCBITaxonomy(taxdump_tgz_file, names_dmp_file, nodes_dmp_file, merged_dmp_file)
    elif os.path.isfile(taxdump_tgz_file):
        loadNCBITaxonomy(taxdump_tgz_file, names_dmp_file, nodes_dmp_file, merged_dmp_file)

    if os.path.isfile(taxonomy_file):
        logger.info( "Open taxonomy file: %s"% taxonomy_file )
        loadTaxonomyTSV(taxonomy_file)

    # try to load custom taxonomy from taxonomy.custom.tsv
    if os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format=='tsv'):
        logger.info( "Open custom taxonomy node file (tsv format): %s"% cus_taxonomy_file)
        loadTaxonomyTSV(cus_taxonomy_file)
    # try to load custom taxonomy from lineage file
    elif os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format=='mgnify_lineage'):
        logger.info( "Open custom taxonomy node file (lineage format): %s"% cus_taxonomy_file)
        loadMgnifyTaxonomy(cus_taxonomy_file)
    # try to load custom taxonomy from GTDB file
    elif os.path.isfile(cus_taxonomy_file) and (cus_taxonomy_format in ['gtdb_taxonomy','gtdb_metadata']):
        loadGTDBTaxonomy(cus_taxonomy_file, cus_taxonomy_format)
    elif os.path.isfile(cus_taxonomy_file):
        logger.fatal( f"invalid cus_taxonomy_format: {cus_taxonomy_format}" )
        _die(f"[ERROR] Invalid cus_taxonomy_format: {cus_taxonomy_format}")

def NCBITaxonomyDownload(dir=None, taxdump=True, accession2taxid=False):
    global taxonomy_dir

    if not dir:
        dir = taxonomy_dir
        logger.info( f"No destination taxonomy dir input. Use default: {dir}..." )

    if not os.path.exists( dir ):
        os.makedirs( dir )
        logger.info( f"Taxonomy dir doesn't exist. Make dir: {dir}..." )

    if taxdump:
        url = 'http://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz'
        logger.info( f"Auto downloading taxanomy from {url}..." )
        # download taxonomy file if auto_download enabled
        r = requests.get(url)
        
        taxdump_tgz_file = f'{dir}/taxdump.tar.gz'

        with open(taxdump_tgz_file, 'wb') as f:
            f.write(r.content)
        if os.path.getsize( taxdump_tgz_file ):
            logger.info( f"Saved to {taxdump_tgz_file}." )
        else:
            logger.fatal( f"Failed to download or save taxonomy files." )
            _die( "[ERROR] Failed to download or save taxonomy files." )    

        # extract
        tax_tar = tarfile.open(taxdump_tgz_file, "r:gz")
        logger.info( f"Extracting nodes.dmp..." )
        tax_tar.extract('nodes.dmp', dir)
        logger.info( f"Extracting names.dmp..." )
        tax_tar.extract('names.dmp', dir)
        logger.info( f"Extracting merged.dmp..." )
        tax_tar.extract('merged.dmp', dir)
        tax_tar.close()
        # delete taxdump_tgz_file
        os.remove(taxdump_tgz_file)

    if accession2taxid:
        import subprocess
        url = "rsync://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid"
        
        logger.info( f"Auto ryncing accession2taxid data from {url}..." )
        subprocess.call(["rsync", "-auvh", "--exclude='prot.accession2taxid.gz*'", "--exclude='prot.accession2taxid.FULL.*.gz*'", url, dir])

        logger.info( f"Decompressing accession2taxid data..." )
        subprocess.call(["gzip", "-d", f"{dir}/accession2taxid/*.gz"])


def loadTaxonomyTSV(tsv_taxonomy_file):

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path)

    try:
        with open(tsv_taxonomy_file) as f:
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
            logger.info( f"Done parsing custom tsv taxonomy file." )
    except IOError:
        _die( "Failed to open custom tsv taxonomy file: %s." % tsv_taxonomy_file )

def loadNCBITaxonomy(taxdump_tgz_file: Optional[str] = None, 
                     names_dmp_file: Optional[str] = None, 
                     nodes_dmp_file: Optional[str] = None, 
                     merged_dmp_file: Optional[str] = None):

    # loading major levels from json file
    _loadAbbrJson(abbr_json_path)

    # try to load taxonomy from taxonomy.tsv
    if os.path.isfile( nodes_dmp_file ) and  os.path.isfile( names_dmp_file ):
        try:
            # read name from names.dmp
            logger.info( f"Open taxonomy name file: {names_dmp_file}" )
            with open(names_dmp_file) as f:
                for line in f:
                    tid, name, tmp, nametype = line.rstrip('\r\n').split('\t|\t')
                    if not nametype.startswith("scientific name"):
                        continue
                    taxNames[tid] = name
                f.close()
                logger.info( f"Done parsing taxonomy name file." )    

            # read taxonomy info from nodes.dmp
            logger.info( f"Open taxonomy node file: {nodes_dmp_file}" )
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
                logger.info( f"Done parsing taxonomy node file." )
        except IOError:
            _die( "Failed to open taxonomy files (taxonomy.tsv, nodes.dmp and names.dmp)." )
    elif os.path.isfile( taxdump_tgz_file ):
        try:
            logger.info( f"Open taxonomy file: {taxdump_tgz_file}" )
            tar = tarfile.open(taxdump_tgz_file, "r:gz")
            
            # read name from names.dmp
            logger.info( "Extract taxonomy names file: names.dmp" )
            member = tar.getmember("names.dmp")
            f = tar.extractfile(member)
            for line in f.readlines():
                tid, name, tmp, nametype = line.decode('utf8').rstrip('\r\n').split('\t|\t')
                if not nametype.startswith("scientific name"):
                    continue
                taxNames[tid] = name
            f.close()
            
            # read taxonomy info from nodes.dmp
            logger.info( "Extract taxonomy nodes file: nodes.dmp" )
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

            # read taxonomy info from merged.dmp
            logger.info( "Extract taxonomy merged file: merged.dmp" )
            member = tar.getmember("merged.dmp")
            f = tar.extractfile(member)
            for line in f.readlines():
                fields = line.decode('utf8').rstrip('\r\n').split('\t|')
                taxMerged[fields[0]] = fields[1].strip('\t')
            f.close()

        except IOError:
            _die( "Failed to load taxonomy from %s"%taxdump_tgz_file )
    
    #try to load merged taxids
    if os.path.isfile( merged_dmp_file ):
        logger.info( "Open merged taxonomy node file: %s"% merged_dmp_file )
        with open(merged_dmp_file) as f:
            for line in f:
                fields = line.rstrip('\r\n').split('\t|')
                taxMerged[fields[0]] = fields[1].strip('\t')
            f.close()
            logger.info( f"Done parsing merged taxonomy file." )
    
def loadMgnifyTaxonomy(mgnify_taxonomy_file=None):
    """
    loadMgnifyTaxonomy()
    """

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path)

    # try to load custom taxonomy from lineage file
    if os.path.isfile(mgnify_taxonomy_file):
        logger.info( "Open custom taxonomy node file (lineage format): %s"% mgnify_taxonomy_file)
        try:
            with open(mgnify_taxonomy_file) as f:
                for line in f:
                    line = line.rstrip('\r\n')
                    if not line: continue
                    if line.startswith('#'): continue
                    if not line.startswith('sk__'):
                        logger.warn( f"A text line of lineage has to start with 'sk__'...skipped: {line}" )
                        continue

                    temp = line.split(';')
                    p_name = ''
                    rank = ''
                    name = ''

                    for i in range(1, len(temp)+1):
                        # this taxa
                        (rank_abbr, name) = temp[-i].split('__')

                        import re
                        re_taxa = re.compile("^([^_]+)__(.*)$")
                        rank_abbr, name = re_taxa.match(temp[-i]).groups()

                        # for na taxon (no_{rank_abbr}_rank)
                        if name=="": name = p_name

                        # paranet taxa
                        try:
                            p_rank_abbr, p_name = re_taxa.match(temp[-(i+1)]).groups()
                            if p_name=="":
                                p_name = f'{name} - no_{rank_abbr}_rank'
                        except:
                            # for the superkingdom rank, assign parant taxid to 1 (root)
                            p_name = '1'
                            if not '1' in taxRanks: taxRanks['1'] = 'root'
                            if not '1' in taxNames: taxNames['1'] = 'root'

                        if rank_abbr in abbr_to_major_level:
                            rank = abbr_to_major_level[rank_abbr]
                        else:
                            rank = rank_abbr
                            
                        tid = name
                        taxParents[tid] = p_name
                        taxRanks[tid] = rank
                        taxNames[tid] = name
                        if p_name in taxNumChilds:
                            taxNumChilds[p_name] += 1
                        else:
                            taxNumChilds[p_name] = 1
                f.close()
                logger.info( f"Done parsing custom taxonomy file." )
        except IOError:
            _die( "Failed to open custom taxonomy file: %s." % mgnify_taxonomy_file )

    logger.info( f"Done parsing taxonomy files (total {len(taxParents)} taxa loaded)" )

def loadGTDBTaxonomy(gtdb_taxonomy_file=None, gtdb_taxonomy_format="gtdb_metadata"):
    """
    loadGTDBTaxonomy()
    """

    # loading major levels to json file
    _loadAbbrJson(abbr_json_path, type='gtdb')

    # try to load custom taxonomy from GTDB file
    if os.path.isfile(gtdb_taxonomy_file) and (gtdb_taxonomy_format in ['gtdb_taxonomy','gtdb_metadata']):
        logger.info( f"Open custom taxonomy node file ({gtdb_taxonomy_format}): %s"% gtdb_taxonomy_file)
        try:
            with open(gtdb_taxonomy_file) as f:
                for line in f:
                    line = line.rstrip('\r\n')
                    if not line: continue
                    if line.startswith('#'): continue
                    if line.startswith('accession'): continue

                    acc, lineage = '', ''

                    if gtdb_taxonomy_format=='gtdb_taxonomy':
                        try:
                            acc, lineage = line.split('\t')
                        except:
                            logger.fatal( f"Incorrect GTDB taxonomy .tsv format: {gtdb_taxonomy_file}" )
                            _die( f"[ERROR] 2 columns are required for GTDB taxonomy .tsv format: {gtdb_taxonomy_file}" )
                        lineage = f'{lineage};x__{acc}'
                    elif gtdb_taxonomy_format=='gtdb_metadata':
                        temp = line.split('\t')
                        if len(temp)!=110:
                            logger.fatal( f"Incorrect GTDB metadata .tsv format: {gtdb_taxonomy_file}" )
                            _die( f"[ERROR] 110 columns are required for GTDB metadata .tsv format: {gtdb_taxonomy_file}" )
                        acc = temp[0]
                        # col_17: 'gtdb_taxonomy'; col_63: 'ncbi_organism_name'
                        lineage = f'{temp[16]};x__{temp[62]}'
                    else:
                        logger.fatal( f"Incorrect format: {gtdb_taxonomy_format}: {gtdb_taxonomy_file}" )
                        _die( f"[ERROR] Incorrect format: {gtdb_taxonomy_format}: {gtdb_taxonomy_file}" )

                    temp = lineage.split(';')
                    p_name = ''
                    rank = ''
                    name = ''

                    for i in range(1, len(temp)+1):
                        # this taxa
                        import re
                        re_taxa = re.compile("^([^_]+)__(.*)$")
                        rank_abbr, name = re_taxa.match(temp[-i]).groups()
                        if i==1 and rank_abbr=='x':
                            name = f'{name} ({acc})'
                        # for na taxon (no_{rank_abbr}_rank)
                        if name=="": name = p_name

                        # paranet taxa
                        try:
                            p_rank_abbr, p_name = re_taxa.match(temp[-(i+1)]).groups()
                            if p_name=="":
                                p_name = f'{name} - no_{p_rank_abbr}_rank'
                        except:
                            # for the *first* taxa in lineage line (usually superkingdom), assign parant taxid to 1 (root)
                            p_name = '1'
                            if not '1' in taxRanks: taxRanks['1'] = 'root'
                            if not '1' in taxNames: taxNames['1'] = 'root'

                        if rank_abbr=='d':
                            rank = 'superkingdom'
                        elif rank_abbr=='x':
                            rank = 'strain'
                        if rank_abbr in abbr_to_major_level:
                            rank = abbr_to_major_level[rank_abbr]
                        else:
                            rank = rank_abbr
                            
                        tid = acc.split('.')[0] if i==1 and rank=='strain' else name
                        taxParents[tid] = p_name
                        taxRanks[tid] = rank
                        taxNames[tid] = name
                        if p_name in taxNumChilds:
                            taxNumChilds[p_name] += 1
                        else:
                            taxNumChilds[p_name] = 1
                f.close()
                logger.info( f"Done parsing custom taxonomy file." )
        except IOError:
            _die( "Failed to open custom taxonomy file: %s." % gtdb_taxonomy_file )

    logger.info( f"Done parsing taxonomy files (total {len(taxParents)} taxa loaded)" )