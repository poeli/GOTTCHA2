#!/usr/bin/env python3
import sys
import io
import taxonomy as t
import requests
import os
import tarfile


#GLOBAL DICTS
#key: gtdb reference number, value: ncbi accession
gtdb_tax = {}
#key: ncbi accession, gtdb reference number
ncbi_tax = {}
#ranks list
ranks = ['strain','species','genus','family','order','class','phylum','superkingdom']
#rank dictionary with letters for keys
rankDict = {'s':'species','g':'genus','f':'family','o':'order','c':'class','p':'phylum','d':'superkingdom'}
#depthDict
depthDict = {'s':7,'g':6,'f':5,'o':4,'c':3,'p':2,'d':1}
#  key gtdb accession, value ncbi taxid
ncbi_taxid = {}
# key gtdb accession, value ncbi spec taxid
ncbi_spec_taxid = {}
#Class for each Taxon in the taxonomy tree
class Taxon():
    def __init__(self, assigned_id, gtdb_id, parent):
         #assignments
         self.id = assigned_id
         self.gtdb_id = gtdb_id
         self.parent = parent
         splitName = gtdb_id.split("__")
         #map to ncbi accession if exists
         try:
             self.ncbi_id = gtdb_tax[gtdb_id]
         except:
             self.ncbi_id = None

        #not a leaf
         if len(splitName) > 1:
             self.name = splitName[1]
             self.rank = rankDict[splitName[0]]
             self.depth = depthDict[splitName[0]]
        #a leaf
         else:
             self.name = gtdb_id
             self.rank = "strain"
             self.id = gtdb_id
             self.depth = 8

#Taxonomy tree
class Graph():
    def __init__(self):
         self.dictionary = {}
         self.count = 1
    def add_node(self, tid, parent):
        taxon = Taxon(self.count, tid, parent)
        self.count += 1
        if tid not in self.dictionary.keys():
             self.dictionary[tid] = []
        return taxon
    def add_edge(self, node1, node2):
         self.dictionary[node1].append(node2)
    def get_dictionary(self):
        return self.dictionary
    def print_graph(self):
        self.print_graph_iter('root')
    def print_graph_iter(self, item):
        print(item + ": " + str(self.dictionary[item]))
        items = self.dictionary[item]
        for i in items:
             self.print_graph_iter(i)


# key: gtdb reference number, value: Node
nodes = {}
#graph
graph = Graph()
#
nodes_ncbi = {}
#
graph_ncbi = Graph()

#assembly accession to rank
def taxid2rank(tid):
    return nodes[tid].rank

#assembly accession to name
def taxid2name(tid):
    return nodes[nodes[tid].parent].name

def taxid2lineageDICT(tid):
    return taxid2lineage(tid)
#returns taxonomy in the format of a dictionary
def taxid2lineage(tid):
    ret = {}
    gtdb_id =  ncbi_tax[tid]
    node = nodes[gtdb_id]
    parent = node.parent
    taxid = node.id
    i = 0
    ret[ranks[i]] = {}
    ret[ranks[i]]['name'] = taxid
    ret[ranks[i]]['taxid'] = taxid
    while parent != 'root':
        i += 1
        node = nodes[parent]
        parent = node.parent
        name = node.name
        taxid = node.id
        ret[ranks[i]] = {}
        ret[ranks[i]]['name'] = name
        ret[ranks[i]]['taxid'] = taxid
    return ret

def taxid2fullLineage( taxID ):
	n = nodes[taxID]
	fullLineage = ""

	while taxID != '1':
		rank = taxid2ranks(taxID)
		name = taxid2name(taxID)
		if not name: break
		fullLineage += "%s|%s|%s|"%(rank,taxID,name)
		taxID = n.parent

	return fullLineage

#loading wol metadata format (unused method)
def loadWOL(path):
    with open(path) as f:
        header = f.readline().rstrip('\r\n').split('\t')
        ncbi_id = header.index('assembly_accession')
        gtdb_id = header.index('gtdb_id')
        for line in f:
            line = line.rstrip('\r\n').split('\t')
            ncbi_tax[line[ncbi_id]] = line[gtdb_id]

#loading gtdb metadata and taxonomy from a directory path
def loadGTDB(path):
    if path is None:
        url = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/'
        bac_taxonomy = url + 'bac120_taxonomy.tsv'
        bac_metadata = url + 'bac120_metadata.tar.gz'
        arch_taxonomy = url + 'ar122_taxonomy.tsv'
        arch_metadata = url + 'ar122_metadata.tar.gz'
        sys.stderr.write( f"[INFO] Auto downloading taxonomy from {url}...\n" )
        # download taxonomy file if auto_download enabled
        r = requests.get(bac_taxonomy, stream=True)
        open('bac120_taxonomy.tsv', 'wb').write(r.content)
        r=requests.get(bac_metadata, stream=True)
        open('bac120_metadata.tar.gz', 'wb').write(r.content)
        r=requests.get(arch_taxonomy, stream=True)
        open('ar122_taxonomy.tsv', 'wb').write(r.content)
        r=requests.get(arch_metadata, stream=True)
        open('ar122_metadata.tar.gz', 'wb').write(r.content)
        #extract metadatafiles
        tar = tarfile.open('bac120_metadata.tar.gz','r:gz')
        tar.extractall()
        tar.close()
        os.rename('bac120_metadata_r202.tsv','bac120_metadata.tsv')
        tar = tarfile.open('ar122_metadata.tar.gz','r:gz')
        tar.extractall()
        tar.close()
        os.rename('ar122_metadata_r202.tsv','ar122_metadata.tsv')
        path=""


    loadGTDBMetadata(str(path) + 'bac120_metadata.tsv')
    loadGTDBMetadata(str(path) + 'ar122_metadata.tsv')
    loadGTDBtaxonomy(str(path) + 'bac120_taxonomy.tsv')
    loadGTDBtaxonomy(str(path) + 'ar122_taxonomy.tsv')

#Load NCBI Taxonomies from GTDB Metadata file
def loadNCBI(path):
    if path is None:
        url = 'https://data.ace.uq.edu.au/public/gtdb/data/releases/latest/'
        bac_metadata = url + 'bac120_metadata.tar.gz'
        arch_metadata = url + 'ar122_metadata.tar.gz'
        sys.stderr.write( f"[INFO] Auto downloading taxanomy from {url}...\n" )
        # download taxonomy file if auto_download enabled
        wget.download(bac_metadata)
        wget.download(arch_metadata)
    loadNCBITaxonomy(str(path) + 'bac120_metadata.tsv')
    loadNCBITaxonomy(str(path) + 'ar122_metadata.tsv')

def loadNCBITaxonomy(metadata):
    if not os.path.isfile(metadata):
        raise Exception(metadata + "File Not Found")
    graph_ncbi.add_node('root', None)
    with open(metadata, encoding="utf-8") as f:
        header = f.readline().rstrip('\r\n').split('\t')
        gtdb_id = header.index("gtdb_genome_representative")
        ncbi_tax = header.index("ncbi_taxonomy")
        for line in f:
            line = line.rstrip('\r\n').split('\t')
            ncbi_id = '_'.join(line[gtdb_id].split("_")[1:])
            tax = line[ncbi_tax].split(";")
            tax.append(line[gtdb_id])
            parent = 'root'
            for i,c in enumerate(tax):
                taxon = graph_ncbi.add_node(c,parent)
                parent = c
                nodes_ncbi[c] = taxon
                if i == 0:
                    graph_ncbi.add_edge('root',c)
                else:
                    graph_ncbi.add_edge(tax[i-1],c)
                if not line:
                    continue

def gtdb2CustomDB(p):
    cus_taxonomy_file = open(p,"w")
    for tid in ncbi_tax:
        gtdb_id =  ncbi_tax[tid]
        node = nodes[gtdb_id]
        parent = node.parent
        taxid = node.id
        cus_taxonomy_file.write('\t'.join((str(taxid),str(node.depth),str(nodes[node.parent].id),node.rank,node.name)))
        while parent != 'root':
            node = nodes[parent]
            parent = node.parent
            cus_taxonomy_file.write(taxid+"\t"+node.depth+"\t"+nodes[node.parent].id+"\t"+node.rank+"\t"+node.name+"\n")
    cus_taxonomy_file.close()

#load metadata file
def loadGTDBMetadata(metadata):
    if not os.path.isfile(metadata):
        raise Exception(metadata + " File Not Found")
    with open(metadata, encoding="utf-8") as f:
        header = f.readline().rstrip('\r\n').split('\t')
        gtdb_id = header.index("gtdb_genome_representative")
        tax = header.index("ncbi_taxid")
        spec_tax = header.index("ncbi_species_taxid")
        for line in f:
            line = line.rstrip('\r\n').split('\t')
            ncbi_id = '_'.join(line[gtdb_id].split("_")[1:])
            gtdb_tax[line[gtdb_id]] = ncbi_id
            ncbi_tax[ncbi_id] = line[gtdb_id]
            ncbi_taxid[ncbi_id] = line[tax]
            ncbi_spec_taxid[ncbi_id] = line[spec_tax]

#load taxonomy file
def loadGTDBtaxonomy(taxonomy):
    if not os.path.isfile(taxonomy):
        raise Exception(taxonomy + "File Not Found")
    graph.add_node('root', None)
    with open(taxonomy, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip('\r\n')
            gtdb_id, taxa = line.split('\t',1)
            tax = taxa.split(';')
            tax.append(gtdb_id)
            parent = 'root'
            for i,c in enumerate(tax):
                taxon = graph.add_node(c, parent)
                parent = c
                nodes[c] = taxon
                if i == 0:
                     graph.add_edge('root', c)
                else:
                     graph.add_edge(tax[i-1], c)
            if not line:
                continue

def taxid2lineageDEFAULT(taxid):
    ret = None
    try:
        ret = taxid2lineage(taxid)
    except:
        try:
            ret = t.taxid2lineageDICT(taxid)
            if ret == "unknown":
                ret = t.taxid2lineageDICT(ncbi_spec_taxid[taxid])
        except:
            raise Exception('Key Error')
    print(taxid)
    print(ret)
    return ret

if __name__ == '__main__':
    loadGTDB(sys.argv[1] if len(sys.argv)>1 else None)
    loadNCBI(sys.argv[1] if len(sys.argv)>1 else None)
