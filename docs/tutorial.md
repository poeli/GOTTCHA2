[![logo](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)

# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA2)

GOTTCHA2 is a gene-independent, signature-based metagenomic taxonomic profiler for sequencing reads. It is designed to reduce false discoveries while remaining practical to run on a workstation or laptop. Instead of relying on marker genes, it maps reads to precomputed unique signature fragments and estimates abundance from signature coverage and depth.

> GOTTCHA v1 databases are not compatible with GOTTCHA2.

---

## Table of contents

- [What's new](#whats-new)
- [Installation](#installation)
- [Dependencies](#dependencies)
- [Databases](#databases)
- [Quick start](#quick-start)
- [Command overview](#command-overview)
- [Profiling](#profiling)
- [Fast profile mode](#fast-profile-mode)
- [Read extraction](#read-extraction)
- [Output files](#output-files)
- [Thresholds and filtering](#thresholds-and-filtering)
- [Full report fields](#full-report-fields)
- [Troubleshooting](#troubleshooting)
- [License and citation](#license-and-citation)

---

## What's new

The current development version is v2.5.0. It includes several workflow changes that are worth knowing before you start:

- **Direct Oxford Nanopore profiling**: `-np/--nanopore` now maps intact ONT reads by default and resolves competing alignments at the species level. The earlier 150 bp chunk workflow remains available with `--ont-chunk`.
- **ONT mapping controls**: direct mode exposes the maximum number of secondary candidates, their minimum score ratio, and the minimum species support through `--ont-max-secondary`, `--ont-secondary-ratio`, and `--ont-min-species-support`.
- **Fast prefiltering mode**: `fast-profile` uses `sylph` to prefilter the reference set before read mapping while producing results comparable to the standard `profile` workflow. Depending on the sample and database, it often reduces runtime by about 5–10× and memory usage by roughly 2–10×.
- **Current CLI**: the supported entry points are `profile`, `fast-profile`, `extract`, `sam2bam`, `download`, and `version`.
- **Updated identity handling**: the reported `SNI_SCORE` is based on consensus identity rather than the legacy read-weighted identity metric.
- **BAM-based workflow**: runs use sorted and indexed BAM for downstream processing instead of keeping SAM as the main intermediate.
- **Legacy compatibility**: the older `gottcha2.py` workflow (SAM-based) is still available for compatibility, but it is frozen at v2.2.3.

---

## Installation

### Option A: Conda

```bash
conda install -c bioconda gottcha2
```

### Option B: Install from source

Install the external tools first, then install the Python package:

```bash
# required for profile
# add sylph as well if you plan to use fast-profile
git clone https://github.com/poeli/GOTTCHA2
cd GOTTCHA2
python -m pip install .

# development install
python -m pip install -e .
```

Confirm the installation:

```bash
gottcha2 version
gottcha2 profile --help
gottcha2 fast-profile --help
```

For containerized usage, see [DOCKER.md](../DOCKER.md).

---

## Dependencies

GOTTCHA2 requires Python 3.9+.

Runtime dependencies:

- `minimap2` 2.27 or newer for mapping
- `samtools` and `pysam` for BAM conversion and parsing
- `numpy` and `pandas`
- `requests`
- `tqdm`
- `biom-format` if you use `--format biom`
- `sylph` if you use `fast-profile`

A Conda environment file is provided as `environment.yml`.

---

## Databases

### Prebuilt databases

The database bundles can be downloaded using `gottcha2 download`. The default download target is:

```text
https://ref-db.edgebioinformatics.org/gottcha2/RefSeq-GTDB-multi-domain+/gottcha_db_fast.tar
```

This is the 32 GB prebuilt GOTTCHA2 database for the `gottcha2 fast-profile` mode. If you prefer, you can also download the database bundles manually from the same host. For the standard profiling mode, the 198 GB prebuilt database [`gottcha_db_standard.tar`](https://ref-db.edgebioinformatics.org/gottcha2/RefSeq-GTDB-multi-domain%2B/gottcha_db_standard.tar), is also available for download.

### Database bundle contents

A standard profiling database should include these files with the same prefix:

- `gottcha_db.<level>.fna.mmi` for `profile`
- `gottcha_db.<level>.fna.tax.tsv` taxonomy mapping
- `gottcha_db.<level>.fna.stats` signature and genome statistics

Fast mode uses the shared `.tax.tsv` and `.stats` files above together with:

- `gottcha_db.<level>.fna.syldb` `sylph` database for prefiltering
- `gottcha_db.<level>.fna.zip` archived signature sequences used to build the reduced reference

The `.mmi` index is required by `profile` but is not used by `fast-profile`.

Pass either the shared database prefix or the database directory to `-d/--database`. GOTTCHA2 locates the required sidecar files from that path. For example:

```text
/path/to/db/gottcha_db.species.fna
```

or

```text
/path/to/db
```

### Download helper

Use the built-in downloader to fetch and verify a database tarball, then extract it into a new `database/` directory. Fast mode is the default download:

```bash
gottcha2 download
```

Download the larger standard database instead with:

```bash
gottcha2 download -d standard
```

The downloader stops if a `database/` directory already exists in the current working directory. It verifies the archive against the published SHA-256 checksum before extraction.

See available options with:

```bash
gottcha2 download --help
```

---

## Quick start

These examples use the most common workflows. Replace the example paths and filenames with your own sample and database locations.

### Example conventions

- `gottcha2 profile` maps reads to the selected signature database and produces taxonomic reports.
- `gottcha2 fast-profile` first narrows the reference set, then runs the profiling workflow on the reduced reference.
- `gottcha2 extract` pulls reads assigned to selected taxa from an existing GOTTCHA2 BAM file.
- `-d/--database` points to either a database prefix, such as `/path/to/db/gottcha_db.species.fna`, or to a directory that contains the matching database files.
- `-i/--input` supplies one or more read files. Use two files for paired-end Illumina reads and one file for single-end or Nanopore reads.
- `-b/--bam` reuses an existing sorted and indexed BAM instead of remapping reads.
- `-np/--nanopore` enables ONT mode and maps intact long reads by default.
- `--ont-chunk`, used together with `-np`, selects the earlier 150 bp chunk workflow.
- `-t/--threads` controls the number of CPU threads used by mapping and related processing steps.
- `-o/--outdir` chooses the output directory. GOTTCHA2 creates it if needed.
- `-p/--prefix` sets the output filename prefix. If you omit it, GOTTCHA2 derives a prefix from the input filename or BAM name.
- The backslash (`\`) at the end of a line lets long shell commands continue on the next line. You can also write each example as a single line.

### 1) Profile Illumina paired-end reads

Use this when your sample has forward and reverse FASTQ files.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample_R1.fastq.gz sample_R2.fastq.gz \
  -t 8 \
  -o out \
  -p sample
```

What this command does:

- loads the species-level database specified with `-d`
- maps both paired-end read files supplied after `-i`
- uses 8 threads because of `-t 8`
- writes results into the `out/` directory
- names output files with the `sample` prefix, for example `sample.tsv`, `sample.full.tsv`, and `sample.gottcha_species.bam`

Use a sample-specific prefix whenever you process multiple samples into the same output directory.

### 2) Profile Illumina single-end reads

Use this when each sample has one FASTQ file.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out
```

Because `-p/--prefix` is omitted, GOTTCHA2 derives the output prefix from `sample.fastq.gz`. Add `-p sample_name` if you want a shorter or more explicit prefix.

### 3) Profile Oxford Nanopore reads

Nanopore mode expects exactly one input file. Add `-np` (short for `--nanopore`) so GOTTCHA2 maps the intact long reads with ONT-oriented settings.

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  -np \
  -t 8 \
  -o out \
  -p ont_sample
```

By default, `-np` uses direct ONT mode: it maps intact reads, retains competing candidate alignments, and keeps alignments supported by a consistent species assignment for each read. To use the earlier chunk-based workflow instead, add `--ont-chunk`:

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  -np --ont-chunk \
  -t 8 \
  -o out \
  -p ont_sample.chunk
```

See [Oxford Nanopore mode](#oxford-nanopore-mode) for the different defaults and tuning options.

### 4) Re-run profiling from an existing BAM

Use this when you already have a sorted and indexed GOTTCHA2 BAM and want to re-aggregate results with different thresholds. This avoids the slower read-mapping step.

```bash
gottcha2 profile \
  -b sample.gottcha_species.bam \
  -d /path/to/db/gottcha_db.species.fna \
  -Mc 0.01 \
  -Mr 10 \
  -mi 0.95 \
  -t 8 \
  -o out \
  -p sample.refiltered
```

What the non-default options mean:

- `-b sample.gottcha_species.bam` reads alignments from an existing BAM instead of using `-i` input reads.
- `-Mc 0.01` requires at least 1% signature coverage for abundance calculation.
- `-Mr 10` requires at least 10 mapped reads.
- `-mi 0.95` keeps only matches with at least 95% alignment identity.
- `-p sample.refiltered` keeps the re-filtered output separate from the original run.

The BAM must be coordinate-sorted and indexed. Keep the database path consistent with the database used for the original mapping.

### 5) Run the faster prefiltering workflow

Use `fast-profile` when you want to reduce runtime and memory usage while producing results comparable to the standard `profile` workflow. It does this by preselecting likely reference sequences before read mapping.

```bash
gottcha2 fast-profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i sample.fastq.gz \
  -t 8 \
  -o out \
  -p sample.fast
```

This command requires the fast-mode `.syldb` and `.zip` files plus the shared `.tax.tsv` and `.stats` files. It does not use the standard `.mmi` index. The outputs have the same general structure as `profile`, but mapping is performed against a reduced reference set selected by `sylph`.

### 6) Extract reads for a taxon from an existing BAM

Use `extract` after profiling when you want the reads assigned to one or more taxa. The example below extracts reads assigned to NCBI taxid `562`.

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 562 \
  -o out \
  -p sample.ecoli
```

Key options:

- `-b` points to the GOTTCHA2 BAM created by `profile` or `fast-profile`.
- `-e` selects the taxon or taxa to extract. You can use taxids, names, or `@file` syntax.
- `-o` and `-p` control where the extracted FASTA or FASTQ output is written.

### Check the results

After any profiling run, start with these files:

```bash
ls out
less out/sample.gottcha_species.log
column -t -s $'\t' out/sample.tsv | less -S
column -t -s $'\t' out/sample.full.tsv | less -S
```

The summary report (`*.tsv`, `*.csv`, or `*.biom`) contains taxa that passed the selected filters. The full report (`*.full.tsv`) includes both passing and filtered taxa and records filtering reasons in the `NOTE` column.

---

## Command overview

GOTTCHA2 uses a subcommand-style CLI:

```text
gottcha2 <command> [options]
```

| Command | Use it when you need to... | Typical starting point |
| ------- | -------------------------- | ---------------------- |
| `profile` | Map reads or reuse a BAM and generate taxonomic profiles. | `gottcha2 profile -d DB -i reads.fastq.gz -o out` |
| `fast-profile` | Prefilter the database with `sylph`, then run profiling on a reduced reference set. | `gottcha2 fast-profile -d DB -i reads.fastq.gz -o out` |
| `extract` | Extract reads assigned to one or more taxa from an existing BAM. | `gottcha2 extract -b sample.bam -e 562` |
| `sam2bam` | Convert legacy GOTTCHA2 SAM output into sorted, indexed BAM. | `gottcha2 sam2bam -i sample.sam -o sample.bam` |
| `download` | Download the default database bundle, when supported by your build. | `gottcha2 download` |
| `version` | Print the installed GOTTCHA2 version. | `gottcha2 version` |

Use `--help` after any command to see command-specific options, defaults, and examples:

```bash
gottcha2 profile --help
gottcha2 extract --help
```

---

## Profiling

### Key concepts

GOTTCHA2 profiles metagenomic samples by mapping sequencing reads directly to taxon-specific signature fragments. GOTTCHA2 consolidates alignments across each genome's signature space to compute coverage and depth statistics, then derives an ANI-like metric called the signature nucleotide identity score (`SNI_SCORE`). Genome-level results are subsequently aggregated to higher taxonomic ranks.

### Oxford Nanopore mode

Use `-np/--nanopore` for a single ONT FASTA or FASTQ file. In v2.5.0 development, this selects direct mode by default. Direct mode maps each intact read to the signature database, considers primary, secondary, and supplementary candidate alignments, and resolves competing species before calculating the profile.

For each read, GOTTCHA2 sums the minimap2 alignment scores for each candidate species. It retains the alignments for a species when that species contributes at least 60% of the read's total candidate score by default. This prevents several alignments from one long read from being counted as independent reads while preserving their aligned bases for coverage calculations.

The earlier chunk workflow remains available with `-np --ont-chunk`. It splits reads into non-overlapping 150 bp pieces, drops a trailing piece shorter than 150 bp, maps the pieces as short reads, and removes taxonomically inconsistent chunk assignments after mapping.

| Setting | Direct mode: `-np` | Chunk mode: `-np --ont-chunk` |
| ------- | ------------------ | ----------------------------- |
| Input used for mapping | Intact ONT reads | Non-overlapping 150 bp chunks |
| minimap2 preset | `lr:hq` | `sr` |
| `--matchIdentity` | `0.85` | `0.85` |
| `--matchFraction` | `0` | `0.85` |
| `--matchLength` | `100` bp | `100` bp |
| `--errorRate` | `0.01` | `0.03` |
| Candidate resolution | Alignment-score support by species | Most consistent taxid across chunks |

The lower direct-mode match fraction is intentional: a short signature alignment can cover only a small fraction of an intact long read. An alignment passes this threshold when the aligned span covers the required fraction of either the read or the signature fragment.

#### Direct-mode controls

Most users can keep the defaults. For samples that need different candidate sensitivity or species assignment stringency, direct mode provides:

- `--ont-max-secondary <INT>`: maximum secondary alignments requested per primary alignment. Default: `30`.
- `--ont-secondary-ratio <FLOAT>`: minimum secondary-to-primary minimap2 chaining-score ratio. Default: `0.5`; accepted range: `0` to `1`.
- `--ont-min-species-support <FLOAT>`: minimum fraction of a read's total candidate alignment score required to retain a species. Default: `0.6`; accepted range: `0` to `1`.
- `-xm/--presetx <STR>`: override the minimap2 preset. Direct mode defaults to `lr:hq`; other accepted values are `sr`, `map-pb`, and `map-ont`.
- `--m2options <STR>`: replace the automatically selected minimap2 tuning options. Use this only when you need explicit mapper control.

For example, this command considers up to 50 secondary candidates while requiring 70% species support:

```bash
gottcha2 profile \
  -d /path/to/db/gottcha_db.species.fna \
  -i ont_reads.fastq.gz \
  -np \
  --ont-max-secondary 50 \
  --ont-min-species-support 0.7 \
  -t 8 \
  -o out
```

These direct-mode options do not change the chunk workflow. `--ont-chunk` itself is only valid together with `-np/--nanopore`.

### Signature of interest

Use `--sigList` to provide a text file containing one accession or signature ID per line. This is useful for plasmids, spike-ins, or other targets you want to track during profiling.

Use `--sigListAction` to control how those reads are handled:

- `report_only` keeps all reads and reports the count in `SOI_READ_COUNT`
- `filter_out` removes reads matching listed accessions
- `filter_in` keeps only reads matching listed accessions

### Reporting level and database level

The database level is usually auto-detected from the database prefix or BAM name. For example, `gottcha_db.species.fna` implies `species`, and `sample.gottcha_species.bam` implies `species`.

If auto-detection is not possible, set it explicitly with `-l/--dbLevel`.

---

## Fast profile mode

`fast-profile` is a convenience wrapper for `profile --fast`. It adds a `sylph` prefiltering step before read mapping:

1. Query the `.syldb` database against the input sample.
2. Collect the subset of candidate signatures.
3. Extract those signatures from the `.zip` archive.
4. Map reads only against that reduced reference.

This mode is useful when you need faster execution with a smaller memory footprint. It still produces the standard GOTTCHA2 outputs, including the BAM and summary reports.

Nanopore selection works the same way in fast mode: `fast-profile -np` maps intact ONT reads by default, while `fast-profile -np --ont-chunk` uses the chunk workflow. Direct fast mode automatically uses mapping settings suited to the reduced reference extracted by the prefilter.

---

## Read extraction

GOTTCHA2 can extract reads for one or more taxa from an existing BAM file. Taxa can be provided as:

- comma-separated taxids, for example `-e "666,562"`
- comma-separated taxon names, for example `-e "Vibrio cholerae,Escherichia coli"`
- a file prefixed with `@`, for example `-e "@taxids.txt"`

The `extract` command is shorthand for running `profile` with `--extract` and `--extractOnly`.

### Example usages

Extract reads mapping to taxid `666`:

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 666
```

Extract with explicit match thresholds:

```bash
gottcha2 extract \
  -b sample.gottcha_species.bam \
  -e 666 \
  -mi 0.9 \
  -mf 0.9
```

Extract multiple taxa:

```bash
gottcha2 extract -b sample.gottcha_species.bam -e "1234,5678"
gottcha2 extract -b sample.gottcha_species.bam -e "@taxids.txt"
```

Limit the number of reads per taxon and choose the output format with `:N:FORMAT`:

```bash
# up to 1000 reads per taxon, FASTQ output
gottcha2 extract -b sample.gottcha_species.bam -e "@taxids.txt:1000:fastq"
```

Extract up to 20 representative sequences per profiled reference:

```bash
gottcha2 extract -b sample.gottcha_species.bam -ef
```

### Extracted record format

Each extracted FASTA or FASTQ header encodes the matched reference, interval, taxon, and match statistics:

```text
>{READ_NAME}{MATE}|{REFERENCE}:{START}..{END} LEVEL={LEVEL} NAME={NAME} TAXID={TAXID} AOI={AOI} MG={MG} MI={MI} MF={MF}
```

Field definitions:

- `READ_NAME`: read identifier
- `MATE`: paired-end suffix (`.1`, `.2`, or empty)
- `REFERENCE`: matched reference sequence name
- `START..END`: mapped reference interval (1-based)
- `LEVEL`: extracted taxonomic rank
- `NAME`: extracted taxon name
- `TAXID`: extracted taxonomy ID
- `AOI`: accession-of-interest flag
- `MG`: alignment length
- `MI`: mapping identity
- `MF`: mapping fraction

Example:

```text
>read123.1|chrA|1|300|GCF10000:10..120 LEVEL=species NAME=Escherichia_coli TAXID=562 AOI=False MG=148 MI=98.65 MF=0.99
ACGT...
```

---

## Output files

By default, outputs go to `--outdir` and use a prefix derived from `--prefix`, the first input filename, or the BAM name.

Typical outputs:

- `*.tsv`, `*.csv`, or `*.biom` - summary report at the requested reporting level
- `*.full.tsv` - full report including filtered taxa and notes
- `*.lineage.tsv` - lineage table for qualified taxa
- `*.mpa.tsv` - MetaPhlAn-style output when `--mpa` is enabled
- `*.extract.fasta` or `*.extract.fastq` - extracted reads when `--extract` or `extract` is used
- `*.gottcha_<level>.bam` and `*.bai` - sorted BAM and index for reuse
- `*.gottcha_<level>.log` - run log including thresholds and processing steps

---

## Thresholds and filtering

Most taxonomic cutoffs default to `0` and are disabled unless you set them explicitly. Alignment thresholds are applied by default unless you lower them yourself.

Use `--noCutoff` to disable taxonomic profiling cutoffs. This is equivalent to:

```text
-Mc 0 -Mr 0 -Ml 0 -Mz 0 -ss 0,0,0
```

### Alignment thresholds

- `-mi, --matchIdentity <FLOAT>`
  Minimum alignment identity for a valid match. Default: `0.95` for short reads and `0.85` for both Nanopore workflows.

- `-mf, --matchFraction <FLOAT>`
  Minimum aligned fraction of the read or signature fragment for a valid match. Default: `0.95` for short reads, `0.05` for direct Nanopore mode, and `0.85` for Nanopore chunk mode.

- `-mg, --matchLength <INT>`
  Minimum alignment length in bp. Default: `100`.

### Taxonomic profiling cutoffs

- `-er, --errorRate <FLOAT>`
  Estimated sequencing error rate. Default: `0.005` for short reads, `0.03` for Nanopore mode.

- `-ss, --sniScore <FLOAT>[,<FLOAT>,<FLOAT>]`
  SNI-score thresholds for `other,species,strain`. Default: `0.9,0.95,0.99`.

- `-Mc, --minCov <FLOAT>`
  Minimum signature coverage required for abundance calculation. Default: `0`.

- `-Mr, --minReads <INT>`
  Minimum number of mapped reads. Default: `0`.

- `-Ml, --minLen <INT>`
  Minimum covered signature length. Default: `0`.

- `-Mz, --maxZscore <FLOAT>`
  Maximum z-score for mapped-region depth distribution. Default: `0` (disabled).

Filtered taxa remain visible in `*.full.tsv`, with the reason recorded in `NOTE`.

---

## Full report fields

The full report (`<prefix>.full.tsv`) contains all computed metrics. The summary report contains the qualified rows shown at the requested reporting level.

| Field Name             | Description |
| ---------------------- | ----------- |
| LEVEL                  | Taxonomic rank (`superkingdom` through `strain`) |
| NAME                   | Taxon name |
| TAXID                  | NCBI taxonomy ID |
| READ_COUNT             | Reads mapped to this taxon |
| TOTAL_BP_MAPPED        | Total mapped bases across this taxon's signatures |
| SNI_SCORE              | Signature nucleotide identity used during filtering and aggregation |
| COVERED_SIG_LEN        | Total covered signature length |
| BEST_SIG_COV           | Highest signature coverage among rolled-up members |
| DEPTH                  | Depth of coverage (`TOTAL_BP_MAPPED / TOTAL_SIG_LEN`) |
| REL_ABUNDANCE_GC       | Relative abundance from genomic-content estimate |
| REL_ABUNDANCE          | Relative abundance from the field selected by `--relAbu` |
| PARENT_NAME            | Parent taxon name |
| PARENT_TAXID           | Parent taxonomy ID |
| AOI_READ_COUNT         | Reads matched to `--sigList` entries |
| TOTAL_READ_LEN         | Total aligned read length |
| TOTAL_BP_MISMATCH      | Total mismatched bases |
| TOTAL_BP_INDEL         | Total inserted and deleted bases |
| READ_WT_SNI            | Read-weighted identity estimate |
| CONSENSUS_SEQ_SNI      | Consensus-sequence identity estimate |
| SNI_CI95_LH            | Low and high 95% confidence bounds for identity |
| SIG_COV                | Signature coverage (`COVERED_SIG_LEN / TOTAL_SIG_LEN`) |
| MAPPED_SIG_LEN         | Signature length with at least one mapped read |
| TOTAL_SIG_LEN          | Total signature length for the taxon |
| COVERED_SIG_DEPTH      | Depth across covered signature only |
| COVERED_MAPPED_SIG_COV | Covered fraction of mapped signature |
| ZSCORE                 | Depth-distribution z-score |
| GENOMIC_CONTENT_EST    | Genomic-content estimate |
| ABUNDANCE              | Raw abundance value from `--relAbu` |
| REL_ABUNDANCE_DEPTH    | Relative abundance computed from depth |
| SIG_LEVEL              | Signature rank used for mapping |
| GENOME_COUNT           | Number of rolled-up genomes |
| GENOME_SIZE            | Combined genome size used for GC normalization |
| NOTE                   | Filtering or rollup note |

---

## Troubleshooting

### BAM input must be sorted and indexed

If you provide `-b/--bam`, the BAM must already be coordinate-sorted and indexed.

For legacy GOTTCHA2 SAM output, convert it with:

```bash
gottcha2 sam2bam -i sample.sam -o sample.bam -t 8
```

### Database sidecar files are required

For `profile`, keep the database sidecar files next to the database prefix. At minimum, GOTTCHA2 expects:

```text
<db>.mmi
<db>.tax.tsv
<db>.stats
```

For `fast-profile`, it expects:

```text
<db>.syldb
<db>.zip
<db>.tax.tsv
<db>.stats
```

### `-np/--nanopore` requires one input file

Nanopore mode only accepts a single FASTA or FASTQ input file. If you have multiple files, merge them first or process them separately.

### `--ont-chunk` requires Nanopore mode

The chunk workflow is selected with `-np --ont-chunk`. Using `--ont-chunk` without `-np/--nanopore` is rejected because chunk preprocessing and postprocessing are specific to ONT reads.

### Python and external dependency checks happen at runtime

GOTTCHA2 checks for:

- Python 3.9+
- `minimap2` 2.27 or newer
- `samtools`
- `sylph` when `fast-profile` is used

If one of these tools is missing from `PATH`, the run stops before mapping begins. Install the missing tool or activate the environment that contains it, then rerun the command.

### No taxa are reported

If the summary report is empty, check the full report and the log before rerunning:

- `*.full.tsv` shows filtered taxa and the reason in `NOTE`.
- `*.gottcha_<level>.log` records the thresholds and database files used.
- Lowering `-Mc`, `-Mr`, `-mi`, or `-mf` can increase sensitivity, but may also increase false positives.

### Fast profile cannot find `.syldb` or `.zip`

`fast-profile` requires `.syldb`, `.zip`, `.tax.tsv`, and `.stats` files with a shared prefix. It does not require the standard `.mmi` index. If the fast-mode files are missing, either download a fast-mode-compatible database bundle or use `profile` with the standard `.mmi` database.

### Identity and SNI changed from older releases

Modern GOTTCHA2 releases report `SNI_SCORE` from consensus identity. If you compare output against older `gottcha2.py` runs, expect differences in SNI-related columns and filtering behavior.

---

## License and citation

- License: BSD 3-Clause
- If you use GOTTCHA2 in publications, cite the GOTTCHA or GOTTCHA2 project, the database source, and the exact software version reported by `gottcha2 version`.
