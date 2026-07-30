[![logo](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)](https://gottcha2.readthedocs.io/en/latest/_images/gottcha_icon.png)

# Genomic Origin Through Taxonomic CHAllenge (GOTTCHA)

[![gottcha2](https://github.com/poeli/GOTTCHA2/actions/workflows/gottcha.yml/badge.svg?branch=master)](https://github.com/poeli/GOTTCHA2/actions/workflows/gottcha.yml)
[![bioconda](https://anaconda.org/bioconda/gottcha2/badges/version.svg)](https://anaconda.org/bioconda/gottcha2)


GOTTCHA is an application of a novel, gene-independent and signature-based metagenomic taxonomic profiling
method with significantly smaller false discovery rates (FDR) that is laptop deployable. Our algorithm was
tested and validated on twenty synthetic and mock datasets ranging in community composition and complexity,
was applied successfully to data generated from spiked environmental and clinical samples, and robustly
demonstrates superior performance compared with other available tools.

-------------------------------------------------------------------
## DEPENDENCIES

GOTTCHA2 profiler is written in Python 3 and uses minimap2 to map reads to signature sequences. To run GOTTCHA2, install the following dependencies. A ready-to-use Conda environment file is provided at `environment.yml`.

- Python 3.9+
- minimap2 2.17+
- samtools
- numpy
- pandas
- requests
- biom-format
- pysam
- tqdm
- sylph

-------------------------------------------------------------------
## QUICK START

1. Install the package:

        via conda `conda install -c bioconda gottcha2`

        OR

        Download or git clone GOTTCHA2 from this repository and run:

        `python -m pip install .`

        (For development installs: `python -m pip install -e .`)

2. Download GOTTCHA2 database bundle:

        The default database is used for `gottcha2 fast-profile`:

        `gottcha2 download`

        OR

        The database for the standard mode `gottcha2 profile`:

        `gottcha2 download -d standard`

3. Run GOTTCHA2:

        Fast-profiling mode:

        `gottcha2 fast-profile -d /path/to/db/ -t 8 -i <FASTQ>`

        OR

        Regular mode:

        `gottcha2 profile -d /path/to/db/ -t 8 -i <FASTQ>`

-------------------------------------------------------------------
## RESULT

GOTTCHA2 can output the profiling results in either CSV, TSV or BIOM format.
- summary (.tsv or .csv) - A summary of profiling results (10 columns) in taxonomic ranks breakdown
- full (.tsv or .csv) - A full profiling results including unfiltered profiling results and additional columns
- lineage (.lineage.tsv or .lineage.tsv) - output lineage and abundance of the profiled taxon per line
- extract (.extract[TAXID].fastq) - Extracted reads for a specific taxon.

-------------------------------------------------------------------
## DOCUMENTATION

Please refer to https://github.com/poeli/GOTTCHA2/wiki for more details.

## Notice of Copyright Assertion (O5124)

This program is Open-Source under the BSD-3 License.
 
Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 
Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 
Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 
Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.