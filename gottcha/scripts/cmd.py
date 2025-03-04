#!/usr/bin/env python3
import pull_database
import gottcha2
import sys

def usage():
    """Display usage information for GOTTCHA2 command-line tool."""
    version = gottcha2.__version__
    print(f"""
GOTTCHA2 - Genomic Origin Through Taxonomic CHAllenge v{version}

Usage:
    gottcha2 <command> [options]

Commands:
    profile    Taxonomic profiling of metagenomic reads
              (Map reads to signature database and classify)
    
Examples:
    gottcha2 profile -i reads.fastq -d database/db_prefix

For detailed help on a specific command:
    gottcha2 <command> --help
""")
    sys.exit(1)

def gottcha2_command():
    args = sys.argv[1:]
    if len(args) < 1:
        usage()
    elif args[0] == "profile":
        gottcha2.main(args[1:])
    # elif args[0] == "pull":
    #     pull_database.main(args[1:])
    else:
        print(f"Error: '{args[0]}' is not a valid command")
        usage()
