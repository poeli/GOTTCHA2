# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)

<img src="https://gitlab.com/poeli/GOTTCHA2/badges/master/build.svg">
<img src="https://gitlab.com/poeli/GOTTCHA2/badges/master/coverage.svg">

GOTTCHA is an application of a novel, gene-independent and signature-based metagenomic taxonomic profiling 
method with significantly smaller false discovery rates (FDR) that is laptop deployable. Our algorithm was 
tested and validated on twenty synthetic and mock datasets ranging in community composition and complexity, 
was applied successfully to data generated from spiked environmental and clinical samples, and robustly 
demonstrates superior performance compared with other available tools.

GOTTCHAv2 is currently under development in BETA stage. Pre-built databases for v1 are incompatible with v2.

-------------------------------------------------------------------
## DEPENDENCIES

GOTTCHA2 profiler is written in Python3 and leverage minimap2 to map reads to signature sequences. In order to run GOTTCHA2 correctly, your system requires to have following dependencies installed correctly.

- Python 3.0+
- minimap2 2.1+
- GNU awk

-------------------------------------------------------------------
## QUICK START

1. Download GOTTCHA2 using git clone, cd to the cloned directory, and make a new directory for the database: 

        $ git clone https://gitlab.com/poeli/GOTTCHA2.git
        $ cd GOTTCHA2
        $ mkdir database

2. Download NCBI taxonomy dmp files, species signature index tar file and untar the file.

        $ wget https://edge-dl.lanl.gov/GOTTCHA2/RefSeq-Release90/taxdump.tar.gz
        $ wget https://edge-dl.lanl.gov/GOTTCHA2/RefSeq-Release90/RefSeq-r90.cg.BacteriaViruses.species.fna.tar
        $ tar -xf RefSeq-r90.cg.BacteriaViruses.species.fna.tar -C database
    
3. Run GOTTCHA2:

        $ ./gottcha2.py -d database/RefSeq-r90.cg.BacteriaViruses.species.fna -t 8 -i <FASTQ>

-------------------------------------------------------------------
## RESULT

GOTTCHA2 provides 5 different outputs that can be specified using `--mode` option (-m [summary|full|class|extract|lineage]):
- summary - report a summary of profiling results (10 columns) in taxonomic ranks breakdown 
- full - report full profiling results including unfiltered profiling results and 12 additional columns
- class - output read classifications
- extract - extract reads
- lineage - output lineage and abundance of a taxonomy per line

#### Summary and full reports

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
| 22     | NOTE                  | Only note the reason for being filtered out                                |                                                            | 

-------------------------------------------------------------------
## USAGE

```
usage: gottcha.py [-h] [-i [FASTQ] [[FASTQ] ...]] [-s [SAMFILE]]
                  [-d [MINIMAP2_INDEX]] [-l [LEVEL]] [-ti [FILE]] [-pm <INT>]
                  [-m {summary,full,class,extract,lineage}] [-x [TAXID]]
                  [-r [FIELD]] [-t <INT>] [-o [DIR]] [-p <STR>] [-xm <STR>]
                  [-mc <FLOAT>] [-mr <INT>] [-ml <INT>] [-mh <INT>] [-nc] [-c]
                  [-v] [--silent] [--debug]

Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an annotation-
independent and signature-based metagenomic taxonomic profiling tool that has
significantly smaller FDR than other profiling tools. This program is a
wrapper to map input reads to pre-computed signature databases using minimap2
and/or to profile mapped reads in SAM format. (VERSION: 2.1.3 BETA)

optional arguments:
  -h, --help            show this help message and exit
  -i [FASTQ] [[FASTQ] ...], --input [FASTQ] [[FASTQ] ...]
                        Input one or multiple FASTQ file(s). Use space to
                        separate multiple input files.
  -s [SAMFILE], --sam [SAMFILE]
                        Specify the input SAM file. Use '-' for standard
                        input.
  -d [MINIMAP2_INDEX], --database [MINIMAP2_INDEX]
                        The path of signature database. The database can be in
                        FASTA format or minimap2 index (5 files).
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
  -m {summary,full,class,extract,lineage}, --mode {summary,full,class,extract,lineage}
                        You can specify one of the following output modes:
                        "summary" : report a summary of profiling result;
                        "full" : other than a summary result, this mode will
                        report unfiltered profiling results with more detail;
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
  -xm <STR>, --presetx <STR>
                        The preset option (-x) for minimap2. Default value
                        'sr' for short reads. [default: sr]
  -mc <FLOAT>, --minCov <FLOAT>
                        Minimum linear coverage to be considered valid in
                        abundance calculation [default: 0.005]
  -mr <INT>, --minReads <INT>
                        Minimum number of reads to be considered valid in
                        abundance calculation [default: 3]
  -ml <INT>, --minLen <INT>
                        Minimum unique length to be considered valid in
                        abundance calculation [default: 60]
  -mh <INT>, --minMLRL <INT>
                        Minimum mean linear read length [default: 1]
  -nc, --noCutoff       Remove all cutoffs. This option is equivalent to use
                        [-mc 0 -mr 0 -ml 0].
  -c, --stdout          Write on standard output.
  -v, --version         Print version number.
  --silent              Disable all messages.
  --debug               Debug mode. Provide verbose running messages and keep
                        all temporary files.
```
