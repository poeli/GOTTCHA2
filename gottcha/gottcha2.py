#!/usr/bin/env python3
__version__   = "2.4.3"
__author__    = "Po-E (Paul) Li, B-GEN, Bioscience Division, Los Alamos National Laboratory"
__credits__   = ["Po-E Li", "Anna Chernikov"]
try:
    from .utils import profile
    from .utils import download
    from .utils import sam_to_bam
except ImportError:
    # If the above relative imports fail, try absolute imports (for direct execution)
    from utils import profile
    from utils import download
    from utils import sam_to_bam

import sys

def usage():
    """Display usage information for GOTTCHA2 command-line tool."""
    print(f"""
GOTTCHA2 - Genomic Origin Through Taxonomic CHAllenge v{__version__}

Usage:
    gottcha2 <command> [options]

Commands:
    profile       Use GOTTCHA2 to profile metagenomic reads against a signature database

    fast-profile  Faster version of profile that uses a more aggressive prefiltering strategy to speed up the read-mapping process

    extract       Extract reads of a specific taxon from profiled results

    sam2bam       Convert GOTTCHA2 SAM to sorted/indexed BAM

    download      Download the latest GOTTCHA2 database

    version       Display version information
    
Examples:
    gottcha2 profile -i reads.fastq -d database/db_prefix

    gottcha2 fast-profile -i reads.fastq -d database/db_prefix

    gottcha2 extract -d prefix.bam -d database/db_prefix -e 666

    gottcha2 download -d fast

    gottcha2 sam2bam -i prefix.sam -o prefix.bam

For detailed help on a specific command:
    gottcha2 <command> --help
""")
    sys.exit(1)

def cli():
    args = sys.argv[1:]
    if len(args) < 1:
        usage()
    elif args[0] == "profile":
        profile.main(args)
    elif args[0] == "fast-profile":
        args.append("--fast")
        profile.main(args)
    elif args[0] == "download":
        download.main(args[1:])
    elif args[0] == "sam2bam":
        sam_to_bam.main(args[1:])
    elif args[0] == "version":
        print(f"{__version__}")
    elif args[0] == "extract":
        args.append("-eo")  # Add --extract-only flag for extract command
        profile.main(args)

    else:
        print(f"Error: '{args[0]}' is not a valid command")
        usage()
