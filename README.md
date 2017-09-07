# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)

GOTTCHA is an application of a novel, gene-independent and signature-based metagenomic
taxonomic profiling method with significantly smaller false discovery rates (FDR) that is 
laptop deployable. Our algorithm was tested and validated on twenty synthetic and mock 
datasets ranging in community composition and complexity, was applied successfully to data
generated from spiked environmental and clinical samples, and robustly demonstrates 
superior performance compared with other available tools.

GOTTCHAv2 is currently under development in BETA stage. Pre-built databases for v1 are incompatible with v2.

-------------------------------------------------------------------
## DISCLOSURE and COPYRIGHT

Copyright (2017).  Los Alamos National Security, LLC. This material was produced
 under U.S. Government contract DE-AC52-06NA25396 for Los Alamos National Labora
tory (LANL), which is operated by Los Alamos National Security, LLC for the U.S.
 Department of Energy. The U.S. Government has rights to use, reproduce, and dis
tribute this software.  NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY,
 LLC MAKES ANY WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE US
E OF THIS SOFTWARE.  If software is modified to produce derivative works, such m
odified software should be clearly marked, so as not to confuse it with the vers
ion available from LANL.

Additionally, this program is free software; you can redistribute it and/or modi
fy it under the terms of the GNU General Public License as published by the Free
 Software Foundation; either version 2 of the License, or (at your option) any l
ater version. Accordingly, this program is distributed in the hope that it will 
be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHA
NTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public Licens
e for more details.

-------------------------------------------------------------------
## SYSTEM REQUIREMENT

Python 3.0+ is required. Linux (2.6 kernel or later) or Mac (OSX 10.6 Snow Leopard or later) operating system with minimal 8 GB of RAM is recommended.

-------------------------------------------------------------------
## QUICK START

0. All available pre-computed databases for RefSeq-Release81 can be found at [LANL's FTP site](ftp://ftp.lanl.gov/public/genome/GOTTCHA2/RefSeq-Release81/).

1. Downloading GOTTCHA2 using git clone and cd to the cloned directory: 

        git clone https://gitlab.com/poeli/GOTTCHA2.git
        cd GOTTCHA2

2. Downloading the database for bacterial species-level identification and untar: 

        wget ftp://ftp.lanl.gov/public/genome/GOTTCHA2/RefSeq-Release81/RefSeq-Release81.Bacteria.species.fna.tar
        tar -xf RefSeq-Release81.Bacteria.species.fna.tar -C database

3. Downloading taxonomy information and untar:

        wget ftp://ftp.lanl.gov/public/genome/GOTTCHA2/RefSeq-Release81/taxonomy.tar
        tar -xf taxonomy.tar -C database

4. Run GOTTCHA2:

        ./gottcha.py -d database/RefSeq-Release81.Bacteria.species.fna -t <THREADS> -i <FASTQ>

-------------------------------------------------------------------
## RESULT

GOTTCHA2 provides 5 different outputs that can be specified using `--mode` option (-m [summary|full|tree|class|extract|lineage]):
- "summary" : report a summary of profiling result (10 columns) in taxonomic ranks breakdown 
- "full" : other than a summary result, a full report including unfiltered profiling results and 12 additional columns
- "tree" : summary report with taxonomic lineage breakdown
- "class" : output read classifications
- "extract" : extract reads
- "lineage" : output lineage and abundance of a taxonomy per line

#### summary and full reports

A full GOTTCHA2 report has 22 columns in tab-delimited format. The summary report is a brief version that has the first 10 columns and qualified taxonomies. The report lists profiling results at taxonomic rank breakdown from superkingdom to strain, following by other information listed below. The rollup depth of coverage (ROLLUP_DOC) is used to calculate relative abundance (column 10) by default, as well as other relative abundance calculations (column 19-21).

| COLUMN | NAME                  | DESCRIPTION                                                                | NOTE                                                       | 
|--------|-----------------------|----------------------------------------------------------------------------|------------------------------------------------------------| 
| 1      | LEVEL                 | Taxonomic rank                                                             |                                                            | 
| 2      | NAME                  | Taxonomic name                                                             |                                                            | 
| 3      | TAXID                 | Taxonomic ID                                                               |                                                            | 
| 4      | READ_COUNT            | Number of mapped reads                                                     |                                                            | 
| 5      | TOTAL_BP_MAPPED       | Total bases of mapped reads                                                |                                                            | 
| 6      | TOTAL_BP_MISMATCH     | Total mismatch bases of mapped reads                                       |                                                            | 
| 7      | LINEAR_LENGTH         | Number of non-overlapping bases covering the signatures                    |                                                            | 
| 8      | LINEAR_DOC            | Linear depth-of-coverage                                                   | = TOTAL_BP_MAPPED / LINEAR_LENGTH                          | 
| 9      | ROLLUP_DOC            | Rollup depth-of-coverage                                                   | = ∑ DOC of sub-level                                       | 
| 10     | REL_ABUNDANCE         | Relative abundance (normalized abundance)                                  | = ABUNDANCE / ∑ ABUNDANCE of given level                   | 
| 11     | LINEAR_COV            | Proportion of covered signatures to total signatures of mapped organism(s) | = LINEAR_LENGTH / SIG_LENGTH_TOL                           | 
| 12     | LINEAR_COV_MAPPED_SIG | Proportion of covered signatures to mapped signatures                      | = LINEAR_LENGTH / SIG_LENGTH_MAPPED                        | 
| 13     | BEST_LINEAR_COV       | Best linear coverage of corresponding taxons                               |                                                            | 
| 14     | DOC                   | Average depth-of-coverage                                                  | = TOTAL_BP_MAPPED / SIG_LENGTH_TOL                         | 
| 15     | BEST_DOC              | Best DOC of corresponding taxons                                           |                                                            | 
| 16     | SIG_LENGTH_TOL        | Length of all signatures in mapped organism(s)                             |                                                            | 
| 17     | SIG_LENGTH_MAPPED     | Length of signatures in mapped signature fragment(s)                       |                                                            | 
| 18     | ABUNDANCE             | abundance of the taxon                                                     | value of either ROLLUP_DOC, READ_COUNT or TOTAL_BP_MAPPED  |
| 19     | REL_ABU_ROLLUP_DOC    | Relative abundance calculated with ROLLUP_DOC                              | = ROLLUP_DOC / ∑ ROLLUP_DOC of given level                 | 
| 20     | REL_ABU_READ_COUNT    | Relative abundance calculated with READ_COUNT                              | = READ_COUNT / ∑ READ_COUNT of given level                 | 
| 21     | REL_ABU_TOL_BP_MAPPED | Relative abundance calculated with TOL_BP_MAPPED                           | = TOL_BP_MAPPED / ∑ TOL_BP_MAPPED of given level           | 
| 22     | NOTE                  | Information                                                                |                                                            | 

-------------------------------------------------------------------
## USAGE

```
usage: gottcha.py [-h] (-i [FASTQ] [[FASTQ] ...] | -s [SAMFILE])
                  [-d [BWA_INDEX]] [-l [LEVEL]] [-ti [FILE]] [-pm <INT>]
                  [-m {summary,full,tree,class,extract,lineage}] [-x [TAXID]]
                  [-r [FIELD]] [-t <INT>] [-o [DIR]] [-p <STR>] [-mc <FLOAT>]
                  [-mr <INT>] [-ml <INT>] [-nc] [-c] [--silent]

Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an annotation-
independent and signature-based metagenomic taxonomic profiling tool that has
significantly smaller FDR than other profiling tools. This program is a
wrapper to map input reads to pre-computed signature databases using BWA-MEM
and/or to profile mapped reads in SAM format. (VERSION: 2.2.2 BETA)

optional arguments:
  -h, --help            show this help message and exit
  -i [FASTQ] [[FASTQ] ...], --input [FASTQ] [[FASTQ] ...]
                        Input one or multiple FASTQ file(s). Use space to
                        separate multiple input files.
  -s [SAMFILE], --sam [SAMFILE]
                        Specify the input SAM file. Use '-' for standard
                        input.
  -d [BWA_INDEX], --database [BWA_INDEX]
                        The path of signature database. The database can be in
                        FASTA format or BWA index (5 files).
  -l [LEVEL], --dbLevel [LEVEL]
                        Specify the taxonomic level of the input database. You
                        can choose one rank from "superkingdom", "phylum",
                        "class", "order", "family", "genus", "species" and
                        "strain". The value will be auto-detected if the input
                        database ended with levels (e.g. GOTTCHA_db.species).
  -ti [FILE], --taxInfo [FILE]
                        Specify the path of taxonomy information file
                        (taxonomy.tsv). GOTTCHA2 will try to locate this file
                        when user doesn't specify a path. If '--database'
                        option is used, the program will try to find this file
                        in the directory of specified database. If not, the
                        'database' directory under the location of gottcha.py
                        will be used as default.
  -pm <INT>, --mismatch <INT>
                        Mismatch penalty for BWA-MEM (pass to option -B while
                        BWA-MEM is running). You can use 99 for not allowing
                        mismatch in alignments (except for extreme cases).
                        [default: 5]
  -m {summary,full,tree,class,extract,lineage}, --mode {summary,full,tree,class,extract,lineage}
                        You can specify one of the following output modes:
                        "summary" : report a summary of profiling result;
                        "full" : other than a summary result, this mode will
                        report unfiltered profiling results with more detail;
                        "tree" : report results with lineage of taxonomy;
                        "class" : output results of classified reads;
                        "extract" : extract mapped reads; "lineage" : output
                        abundance and lineage in a line; Note that only
                        results/reads belongs to descendants of TAXID will be
                        reported/extracted if option [--taxonomy TAXID] is
                        specified. [default: summary]
  -x [TAXID], --taxonomy [TAXID]
                        Specify a NCBI taxonomy ID. The program will only
                        report/extract the taxonomy you specified.
  -r [FIELD], --relAbu [FIELD]
                        The field will be used to calculate relative
                        abundance. You can specify one of the following
                        fields: "LINEAR_LENGTH", "TOTAL_BP_MAPPED",
                        "READ_COUNT" and "LINEAR_DOC". [default: ROLLUP_DOC]
  -t <INT>, --threads <INT>
                        Number of threads [default: 1]
  -o [DIR], --outdir [DIR]
                        Output directory [default: .]
  -p <STR>, --prefix <STR>
                        Prefix of the output file [default:
                        <INPUT_FILE_PREFIX>]
  -mc <FLOAT>, --minCov <FLOAT>
                        Minimum linear coverage to be considered valid in
                        abundance calculation [default: 0.005]
  -mr <INT>, --minReads <INT>
                        Minimum number of reads to be considered valid in
                        abundance calculation [default: 3]
  -ml <INT>, --minLen <INT>
                        Minimum unique length to be considered valid in
                        abundance calculation [default: 60]
  -nc, --noCutoff       Remove all cutoffs. This option is equivalent to use
                        [-mc 0 -mr 0 -ml 0].
  -c, --stdout          Write on standard output.
  --silent              Disable all messages.
  ```