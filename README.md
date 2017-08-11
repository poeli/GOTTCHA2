# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)

GOTTCHA is an application of a novel, gene-independent and signature-based metagenomic
taxonomic profiling method with significantly smaller false discovery rates (FDR) that is 
laptop deployable. Our algorithm was tested and validated on twenty synthetic and mock 
datasets ranging in community composition and complexity, was applied successfully to data
generated from spiked environmental and clinical samples, and robustly demonstrates 
superior performance compared with other available tools.

-------------------------------------------------------------------
## DISCLOSURE

GOTTCHAv2 is currently under development in ALPHA stage. Databases for v1 will not be compatible with v2.

-------------------------------------------------------------------
## SYSTEM REQUIREMENT

Python 3.0+ is required. Linux (2.6 kernel or later) or Mac (OSX 10.6 Snow Leopard or later) operating system with minimal 8 GB of RAM is recommended.

-------------------------------------------------------------------
## RESULT
GOTTCHA reports profiling results in plain text table (*.gottcha.tsv) by default. The tsv file will list the organism(s) at all taxonomic levels from STRAIN to PHYLUM, following by other information listed below. The linear depth of coverage (LINEAR_DOC) is used to calculate relative abundance of each organism or taxonomic name in the sample.

| COLUMN | NAME                  | DESCRIPTION                                                                | EQUIVALENT                                                 | 
|--------|-----------------------|----------------------------------------------------------------------------|------------------------------------------------------------| 
| 1      | LEVEL                 | Taxonomic rank                                                             |                                                            | 
| 2      | NAME                  | Taxonomic name                                                             |                                                            | 
| 3      | TAXID                 | Taxonomic ID                                                               |                                                            | 
| 4      | READ_COUNT            | Number of mapped reads                                                     |                                                            | 
| 5      | TOTAL_BP_MAPPED       | Total bases of mapped reads                                                |                                                            | 
| 6      | TOTAL_BP_MISMATCH     | Total mismatch bases of mapped reads                                       |                                                            | 
| 7      | LINEAR_LENGTH         | Number of non-overlapping bases covering the signatures                    |                                                            | 
| 8      | LINEAR_DOC            | Linear depth-of-coverage                                                   | TOTAL_BP_MAPPED / LINEAR_LENGTH                            | 
| 9      | REL_ABUNDANCE         | Normalized abundance                                                       | ABUNDANCE / ∑ ABUNDANCE (for each level)                   | 
| 10     | LINEAR_COV            | Proportion of covered signatures to total signatures of mapped organism(s) | LINEAR_LENGTH / SIG_LENGTH_TOL                             | 
| 11     | LINEAR_COV_MAPPED_SIG | Proportion of covered signatures to mapped signatures                      | LINEAR_LENGTH / SIG_LENGTH_MAPPED                          | 
| 12     | BEST_LINEAR_COV       | Best linear coverage of corresponding taxons                               |                                                            | 
| 13     | DOC                   | Average depth-of-coverage                                                  | TOTAL_BP_MAPPED / SIG_LENGTH_TOL                           | 
| 14     | BEST_DOC              | Best DOC of corresponding taxons                                           |                                                            | 
| 15     | SIG_LENGTH_TOL        | Length of all signatures in mapped organism(s)                             |                                                            | 
| 16     | SIG_LENGTH_MAPPED     | Length of signatures in mapped signature fragment(s)                       |                                                            | 
| 17     | COPY                  | Mimic cell copy number                                                     | ∑ DOC                                                      | 
| 18     | ABUNDANCE             | abundance of the taxon                                                     | [READ_COUNT|TOTAL_BP_MAPPED|LINEAR_LENGTH|LINEAR_DOC|COPY] | 
| 19     | NOTE                  | Information                                                                |                                                            | 

-------------------------------------------------------------------
## INSTRUCTIONS

```python
usage: gottcha.py [-h] (-i [FASTQ] [[FASTQ] ...] | -s [SAMFILE])
                  [-d [BWA_INDEX]] [-l [LEVEL]] [-pm <INT>]
                  [-m {summary,full,tree,class,extract}] [-x [TAXID]]
                  [-r [FIELD]] [-t <INT>] [-o [DIR]] [-p <STR>] [-mc <FLOAT>]
                  [-mr <INT>] [-ml <INT>] [-c] [-v]
```
Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an annotation-
independent and signature-based metagenomic taxonomic profiling tool that has
significantly smaller FDR than other profiling tools. This program is a
wrapper to map input reads to pre-computed signature databases using BWA-MEM
and/or to profile mapped reads in SAM format. (VERSION: 2.0 Beta)

optional arguments:

```
-i [FASTQ] [[FASTQ] ...], --input [FASTQ] [[FASTQ] ...]
```
Input one or multiple FASTQ file(s). Use space to separate multiple input files.
```
-s [SAMFILE], --sam [SAMFILE]
```
Specify the input SAM file. Use '-' for standard input.
```
-d [BWA_INDEX], --database [BWA_INDEX]
```
The path of signature database. The database can be in FASTA format or BWA index (5 files).
```
-l [LEVEL], --dbLevel [LEVEL]
```
Specify the taxonomic level of the input database. You can choose one rank from "superkingdom", "phylum","class", "order", "family", "genus", "species" and "strain". The value will be auto-detected if the input database ended with levels (e.g. GOTTCHA_db.species).
```
-pm <INT>, --mismatch <INT>
```
Mismatch penalty for BWA-MEM (pass to option -B while BWA-MEM is running). You can use 99 for not allowing mismatch in alignments (except for extreme cases). [default: 5]
```
-m {summary,full,tree,class,extract}, --mode {summary,full,tree,class,extract}
```
You can specify one of the following output modes:

* "summary" : report a summary of profiling result;
* "full" : other than a summary result, this mode will report unfiltered profiling results with more detail;
* "tree" : report results with lineage of taxonomy;
* "class" : output results of classified reads;
* "extract" : extract mapped reads; Note that only results/reads belongs to descendants of
TAXID will be reported/extracted if option [--taxonomy TAXID] is specified. [default: summary]

```
-x [TAXID], --taxonomy [TAXID]
```
Specify a NCBI taxonomy ID. The program will only report/extract the taxonomy you specified.
```
-r [FIELD], --relAbu [FIELD]
```
The field will be used to calculate relative abundance. You can specify one of the following fields: "LINEAR_LENGTH", "TOTAL_BP_MAPPED", "READ_COUNT" and "LINEAR_DOC". [default: LINEAR_DOC]
```
-t <INT>, --threads <INT>
```
Number of threads [default: 1]
```
-o [DIR], --outdir [DIR]
```
Output directory [default: .]
```
-p <STR>, --prefix <STR>
```
Prefix of the output file [default:<INPUT_FILE_PREFIX>]
```
-mc <FLOAT>, --minCov <FLOAT>
```
Minimum linear coverage to be considered valid in abundance calculation [default: 0.005]
```
-mr <INT>, --minReads <INT>
```
Minimum number of reads to be considered valid in abundance calculation [default: 3]
```
-ml <INT>, --minLen <INT>
```
Minimum unique length to be considered valid in abundance calculation [default: 60]
```
-c, --stdout
```
Write on standard output.
```
-v, --verbose
```
Enable verbose output
```
-h, --help
```
show this help message and exit
