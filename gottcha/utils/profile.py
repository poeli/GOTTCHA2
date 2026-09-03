#!/usr/bin/env python3

import argparse as ap
from asyncio import threads
import re
import sys, os, time, subprocess
import pandas as pd
from pathlib import Path
import gc
from re import search
from multiprocessing import set_start_method
import logging

try:
    # Try relative import first (for package usage)
    import taxonomy
    import report
    import sam_to_bam
    import process_bam
    import ont_utils
    import read_mapping
    import aggregate_results
    import extract_reads
    import gottcha.utils.prefilter as prefilter
    import sig_archive
    from gottcha2 import __version__
except ImportError:
    # Fall back to direct import (for script usage)
    import gottcha.utils.report as report
    import gottcha.utils.taxonomy as taxonomy
    import gottcha.utils.ont_utils as ont_utils
    import gottcha.utils.sam_to_bam as sam_to_bam
    import gottcha.utils.process_bam as process_bam
    import gottcha.utils.aggregate_results as aggregate_results
    import gottcha.utils.read_mapping as read_mapping
    import gottcha.utils.extract_reads as extract_reads
    import gottcha.utils.prefilter as prefilter
    import gottcha.utils.sig_archive as sig_archive
    from gottcha.gottcha2 import __version__

def parse_args(ver, args):
    """
    Parse and validate command line arguments for GOTTCHA2.

    This function sets up the argument parser, defines all possible command-line
    options, parses the provided arguments, and performs validation to ensure
    the configuration is valid and complete.

    Parameters:
        ver (str): Version string to display in help messages
        args (list): Command line arguments to parse

    Returns:
        argparse.Namespace: Object containing all validated arguments

    Raises:
        SystemExit: If validation fails or --version is specified
    """
    command = args[0] if args and not args[0].startswith('-') else None
    parser_args = args[1:] if command else args
    taxonomic_levels = [
        'superkingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'strain',
    ]

    p = ap.ArgumentParser(
        prog=f'gottcha2 {command}' if command else 'gottcha2 profile',
        formatter_class=ap.RawTextHelpFormatter,
        description=(
            "Genomic Origin Through Taxonomic CHAllenge (GOTTCHA) is an annotation-independent,\n"
            "signature-based metagenomic taxonomic profiling tool with substantially low false\n"
            "discovery rates. This command maps input reads to precomputed signature databases\n"
            f"using minimap2 and profiles the organisms present in a sample. (Version: {ver})"
        ),
    )

    input_group = p.add_argument_group('Input source')
    eg = input_group.add_mutually_exclusive_group(required=True)
    eg.add_argument(
        '-i', '--input',
        metavar='FASTQ',
        nargs='+',
        type=str,
        help='Input FASTQ/FASTA file(s). Separate multiple files with spaces.',
    )
    eg.add_argument(
        '-b', '--bam',
        metavar='BAMFILE',
        type=str,
        help='Input sorted and indexed BAM file.',
    )

    database_group = p.add_argument_group('Database')
    database_group.add_argument(
        '-d', '--database',
        metavar='GOTTCHA2_DB',
        type=str,
        default=None,
        help='Path to the GOTTCHA2 database prefix, index file, or database directory.',
    )
    database_group.add_argument(
        '-l', '--dbLevel',
        metavar='LEVEL',
        type=str,
        default='',
        choices=taxonomic_levels,
        help=(
            'Taxonomic level of the input database.\n'
            'Choices: superkingdom, phylum, class, order, family, genus, species, strain.\n'
            'Auto-detected when the database prefix contains a rank, e.g. GOTTCHA_db.species.'
        ),
    )

    platform_group = p.add_argument_group('Sequencing and mapping')
    platform_group.add_argument(
        '-np', '--nanopore',
        action='store_true',
        help=(
            'Treat input reads as Oxford Nanopore (ONT) reads.\n'
            'Enables read preprocessing and ONT defaults: -er 0.03 -mi 0.85 -mf 0.85 -mg 100.'
        ),
    )
    platform_group.add_argument(
        '--ont-direct',
        action='store_true',
        help=(
            'Map intact ONT reads directly to GOTTCHA2 signature fragments and '
            'resolve multiple alignments at species level. Without this option, '
            'the legacy 150-bp chunk workflow is used.'
        ),
    )
    platform_group.add_argument(
        '--ont-max-secondary', type=int, default=50,
        help='Maximum minimap2 secondary candidates per primary in --ont-direct mode. [default: 50]',
    )
    platform_group.add_argument(
        '--ont-secondary-ratio', type=float, default=0.5,
        help='Minimum minimap2 secondary/primary chaining-score ratio (-p) in --ont-direct mode. [default: 0.5]',
    )
    platform_group.add_argument(
        '--ont-min-species-support',
        dest='ont_min_species_support', type=float, default=0.65,
        help='Minimum fraction of competing union-bp support required by the winning species. [default: 0.65]',
    )
    platform_group.add_argument(
        '-xm', '--presetx',
        metavar='STR',
        type=str,
        default=None,
        choices=['sr', 'map-pb', 'map-ont', 'lr:hq'],
        help='minimap2 preset passed with -x. [default: lr:hq for --nanopore --ont-direct; sr otherwise]',
    )
    platform_group.add_argument(
        '--m2options',
        metavar='STR',
        type=str,
        default='auto',
        help='Additional minimap2 options. Use with care. [default: platform-specific auto settings]',
    )
    platform_group.add_argument(
        '-er', '--errorRate',
        metavar='FLOAT',
        type=float,
        help='Estimated sequencing error rate. [default: 0.005, or 0.03 with --nanopore]',
    )
    platform_group.add_argument(
        '-t', '--threads',
        metavar='INT',
        type=int,
        default=1,
        help='Number of threads. [default: 1]',
    )

    alignment_group = p.add_argument_group('Alignment thresholds')
    alignment_group.add_argument(
        '-mi', '--matchIdentity',
        metavar='FLOAT',
        type=float,
        help='Minimum identity (0.0-1.0) required for a valid match. [default: 0.95, or 0.85 with --nanopore]',
    )
    alignment_group.add_argument(
        '-mf', '--matchFraction',
        metavar='FLOAT',
        type=float,
        help='Minimum aligned fraction (0.0-1.0) of the read or signature fragment. [default: 0.95, or 0.85 with --nanopore]',
    )
    alignment_group.add_argument(
        '-mg', '--matchLength',
        metavar='INT',
        type=int,
        help='Minimum alignment length in bp required for a valid match. [default: 100]',
    )

    profiling_group = p.add_argument_group('Profiling thresholds')
    profiling_group.add_argument(
        '-ss', '--sniScore',
        metavar='FLOAT[,FLOAT,FLOAT]',
        type=str,
        help=(
            'Signature nucleotide identity (SNI) thresholds for taxonomic aggregation.\n'
            'One value applies to all ranks; two values append strain default 0.99;\n'
            'three values mean other ranks, species, and strain. [default: 0.9,0.95,0.99]'
        ),
    )
    profiling_group.add_argument(
        '-Mc', '-mc', '--minCov',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Minimum signature coverage used in abundance calculation. [default: 0]',
    )
    profiling_group.add_argument(
        '-Mr', '-mr', '--minReads',
        metavar='INT',
        type=int,
        default=0,
        help='Minimum read count used in abundance calculation. [default: 0]',
    )
    profiling_group.add_argument(
        '-Ml', '-ml', '--minLen',
        metavar='INT',
        type=int,
        default=0,
        help='Minimum signature length used in abundance calculation. [default: 0]',
    )
    profiling_group.add_argument(
        '-Mz', '-mz', '--maxZscore',
        metavar='FLOAT',
        type=float,
        default=0,
        help='Maximum estimated z-score for mapped-region depths; 0 disables the filter. [default: 0]',
    )
    profiling_group.add_argument(
        '-nc', '--noCutoff',
        action='store_true',
        help='Disable profiling-stage cutoffs. Equivalent to -Mc 0 -Mr 0 -Ml 0 -Mz 0 -ss 0,0,0.',
    )
    profiling_group.add_argument(
        '-r', '--relAbu',
        metavar='FIELD',
        type=str,
        default='DEPTH',
        choices=['DEPTH', 'READ_COUNT', 'GENOMIC_CONTENT_EST'],
        help='Field used to calculate relative abundance. [default: DEPTH]',
    )

    signature_group = p.add_argument_group('Signature-of-interest filtering')
    signature_group.add_argument(
        '-sl', '--sigList',
        metavar='FILE',
        type=str,
        help='File containing accessions/signatures of interest, one per line.',
    )
    signature_group.add_argument(
        '-sa', '--sigListAction',
        choices=['filter_out', 'filter_in', 'report_only'],
        default='report_only',
        type=str,
        help=(
            'Action for aligned reads mapped to signatures of interest:\n'
            '  filter_out  discard matching reads\n'
            '  filter_in   keep only matching reads\n'
            '  report_only report matching reads as SOI_READ_COUNT without filtering [default]'
        ),
    )

    extraction_group = p.add_argument_group('Read extraction')
    extraction_group.add_argument(
        '-e', '--extract',
        metavar='TAXON[,TAXON2,...]',
        type=str,
        default=None,
        help=(
            'Extract mapped reads for specific taxa to FASTA or FASTQ.\n'
            'Accepted forms:\n'
            "  1234,5678              comma-separated taxon IDs\n"
            "  @taxids.txt            one taxon ID per line\n"
            "  @taxids.txt:1000:fasta limit reads per taxon and choose fasta/fastq\n"
            "  all                    extract all matching taxa/reads"
        ),
    )
    extraction_group.add_argument(
        '-ef', '--extractFullRef',
        action='store_true',
        help="Extract up to 20 sequences per reference as FASTA. Equivalent to -e 'all:20:fasta'.",
    )
    extraction_group.add_argument(
        '-eo', '--extractOnly',
        action='store_true',
        help='Only extract reads from the alignment file; skip profiling.',
    )

    output_group = p.add_argument_group('Output')
    output_group.add_argument(
        '-fm', '--format',
        metavar='STR',
        type=str,
        default='tsv',
        choices=['tsv', 'csv', 'biom'],
        help='Output format. [default: tsv]',
    )
    output_group.add_argument(
        '-o', '--outdir',
        metavar='DIR',
        type=str,
        default='.',
        help='Output directory. [default: .]',
    )
    output_group.add_argument(
        '-p', '--prefix',
        metavar='STR',
        type=str,
        help='Output file prefix. [default: input file prefix]',
    )
    output_group.add_argument(
        '--mpa',
        action='store_true',
        help='Generate output in MetaPhlAn format (*.mpa).',
    )

    fast_group = p.add_argument_group('Fast profile')
    fast_group.add_argument(
        '--fast',
        action='store_true',
        help='Enable fast-profile mode.',
    )
    fast_group.add_argument(
        '--fast-min-kmer',
        metavar='INT',
        type=int,
        default=5,
        help='Minimum k-mer size for fast-profile prefiltering. [default: 5]',
    )
    fast_group.add_argument(
        '--fast-min-ani',
        metavar='INT',
        type=int,
        default=80,
        help='Minimum ANI for fast-profile prefiltering. [default: 80]',
    )
    
    logging_group = p.add_argument_group('Logging')
    logging_group.add_argument(
        '--silent',
        action='store_true',
        help='Disable status messages.',
    )

    logging_group.add_argument(
        '--verbose',
        action='store_true',
        help='Show verbose status messages.',
    )
    logging_group.add_argument(
        '--debug',
        action='store_true',
        help='Enable debug logging and keep temporary files.',
    )
    
    eg.add_argument(
        '-v', '--version',
        action='store_true',
        help='Print version number and exit.',
    )
    args_parsed = p.parse_args(parser_args)

    """
    Checking options
    """
    if args_parsed.version:
        print(ver)
        sys.exit(0)

    if args_parsed.extract and args_parsed.extractFullRef:
        p.error('--extract and --extractFullRef are incompatible options.')

    if args_parsed.input and args_parsed.bam:
        p.error('--input / --bam are incompatible options.')

    if not args_parsed.extractOnly:
        if not args_parsed.database:
            p.error('--database option is missing.')
        if args_parsed.sniScore is None:
            args_parsed.sniScore='0.9,0.95,0.99'

    # Auto-detect database path and prefix, and check the existence of input files and database index
    if args_parsed.database:
        #assign default path for database name
        db_extfn = "syldb" if args_parsed.fast else "mmi"

        # find the database index file if a directory is provided, and set the database prefix accordingly;
        if Path(args_parsed.database).is_dir():
            dbs = list(Path(args_parsed.database).glob(f"*.{db_extfn}"))
            if len(dbs) > 1:
                p.error(f'Multiple .{db_extfn} files found in {args_parsed.database}. Please specify one with database prefix.')
            elif len(dbs) == 0:
                p.error(f'No .{db_extfn} file found in {args_parsed.database}. Please specify the database prefix or the path to the .{db_extfn} file.')
            else:
                args_parsed.database = str(dbs[0])

        # check if the database file is provided with or without the extension, and remove the extension if provided
        if args_parsed.database.endswith(f'.{db_extfn}'):
            args_parsed.database = args_parsed.database.replace(f'.{db_extfn}', '')

        # Only check the existence of the database index file if input reads are provided
        if args_parsed.input:
            if not Path(f'{args_parsed.database}.{db_extfn}').is_file():
                p.error(f'Database index {args_parsed.database}.{db_extfn} not found.')

        if args_parsed.fast:
            if not Path(f'{args_parsed.database}.{db_extfn}').is_file():
                p.error(f'Database index {args_parsed.database}.{db_extfn} not found.')
            if not Path(f'{args_parsed.database}.zip').is_file():
                p.error(f'Signature sequences file {args_parsed.database}.zip not found.')

        # check the existence of the taxonomic information file for the specified database
        if not Path(f'{args_parsed.database}.tax.tsv').is_file():
            p.error(f'Taxonomic file {args_parsed.database}.tax.tsv not found.')
        if not Path(f'{args_parsed.database}.stats').is_file():
            p.error(f'Database stats file {args_parsed.database}.stats not found.')

    if args_parsed.input:
        for path in args_parsed.input:
            if path == '-':
                p.error('--input does not support reading from stdin ("-"). Please provide a file path.')
            if not Path(path).is_file():
                p.error(f'Input file {path} not found.')

    if args_parsed.bam:
        if not Path(args_parsed.bam).is_file():
            p.error(f'BAM file {args_parsed.bam} not found.')

    if args_parsed.sigList:
        if not Path(args_parsed.sigList).is_file():
            p.error(f'Signature-of-interest list {args_parsed.sigList} not found.')
        args_parsed.sigList = os.path.abspath(args_parsed.sigList)

    if args_parsed.nanopore and args_parsed.input and len(args_parsed.input) != 1:
        p.error('--nanopore option requires a single input read file.')

    if not args_parsed.prefix:
        if args_parsed.input:
            name = search(r'([^\/\.]+)\..*$', args_parsed.input[0])
            args_parsed.prefix = name.group(1)
        elif args_parsed.bam:
            name = search(r'([^\/]+)\.\w+\.bam$', args_parsed.bam)
            args_parsed.prefix = name.group(1)
        else:
            args_parsed.prefix = "GOTTCHA_"

    if not args_parsed.dbLevel:
        if args_parsed.database:
            major_ranks = {"superkingdom":1,"phylum":2,"class":3,"order":4,"family":5,"genus":6,"species":7, "strain":8}
            parts = args_parsed.database.split('.')
            for part in parts:
                if part in major_ranks:
                    args_parsed.dbLevel = part
                    break
        elif args_parsed.bam:
            name = search(r'\.gottcha_(\w+).bam$', args_parsed.bam)
            try:
                args_parsed.dbLevel = name.group(1)
            except:
                pass

        if not args_parsed.dbLevel:
            p.error('--dbLevel is missing and cannot be auto-detected.')

    # Set SNI-SCORE default to 0.9, species 0.95, strain 0.99
    if args_parsed.sniScore:
        if args_parsed.sniScore.count(',') == 0:
            args_parsed.sniScore = ','.join([args_parsed.sniScore]*3)
        elif args_parsed.sniScore.count(',') == 1:
            args_parsed.sniScore = args_parsed.sniScore + ',0.99'
    else:
        args_parsed.sniScore = '0.9,0.95,0.99'

    # If mi/mf/mg are not specified, set them to default values based on whether --nanopore is specified
    # But if --extractOnly is specified, do not set default values for matchIdentity and matchFraction, the value will be load from the log file if not provided
    if args_parsed.matchIdentity is None:
        if not args_parsed.extractOnly:
            if args_parsed.nanopore:
                args_parsed.matchIdentity = 0.85
            else:
                args_parsed.matchIdentity = 0.95
    else:
        if args_parsed.matchIdentity < 0 or args_parsed.matchIdentity > 1:
            p.error('--matchIdentity must be between 0 and 1.')

    if args_parsed.matchFraction is None:
        if not args_parsed.extractOnly:
            if args_parsed.nanopore:
                args_parsed.matchFraction = 0.85
            else:
                args_parsed.matchFraction = 0.95
    else:
        if args_parsed.matchFraction < 0 or args_parsed.matchFraction > 1:
            p.error('--matchFraction must be between 0 and 1.')

    if args_parsed.matchLength is None:
        if not args_parsed.extractOnly:
            args_parsed.matchLength = 100
    else:
        if args_parsed.matchLength < 0:
            p.error('--matchLength must be a non-negative integer.')

    if args_parsed.extractFullRef:
        args_parsed.extract = 'all:20:fasta'

    if args_parsed.extractOnly:
        error_message = ""
        if not args_parsed.extract:
            error_message += "--extract must be specified. "

        if not args_parsed.bam:
            error_message += "--bam must be specified. "

        if error_message:
            p.error(error_message)

    if args_parsed.ont_direct and not args_parsed.nanopore:
        p.error('--ont-direct requires --nanopore.')
    if args_parsed.ont_direct:
        if args_parsed.ont_max_secondary < 0:
            p.error('--ont-max-secondary must be >= 0.')
        if not 0 <= args_parsed.ont_secondary_ratio <= 1:
            p.error('--ont-secondary-ratio must be between 0 and 1.')
        if not 0 <= args_parsed.ont_min_species_support <= 1:
            p.error('--ont-min-species-support must be between 0 and 1.')
        if not args_parsed.matchFraction:
            args_parsed.matchFraction = 0

    if args_parsed.presetx is None:
        args_parsed.presetx = 'lr:hq' if (args_parsed.nanopore and args_parsed.ont_direct) else 'sr'

    if args_parsed.m2options == 'auto':
        if args_parsed.nanopore and args_parsed.ont_direct:
            # Fast mode builds the signature index on the fly and adds k24/w12
            # below, so two minimizers are a useful short-fragment safeguard.
            # Prebuilt GOTTCHA2 .mmi indexes retain k28/w24; allow one seed to
            # initiate DP for ~100-bp signatures.
            seed_chain = '-n2' if args_parsed.fast else '-n1'
            args_parsed.m2options = f'{seed_chain} -m25 -s40 --no-long-join'
        else:
            args_parsed.m2options = '-s120'

    if args_parsed.noCutoff:
        args_parsed.sniScore = '0,0,0'

    if not args_parsed.errorRate:
        if args_parsed.nanopore:
            args_parsed.errorRate = 0.03
        else:
            args_parsed.errorRate = 0.005

    return args_parsed


def dependency_check(cmd: str) -> None:
    """
    Verify that external dependencies are available in the system.

    Attempts to execute the specified command with --help and checks if it runs
    successfully. Exits the program if the command is not found or fails.

    Parameters:
        cmd (str): Command to check

    Returns:
        None

    Raises:
        SystemExit: If the command is not found or fails
    """
    try:
        subprocess.check_call([cmd, "--help"], stdout=subprocess.DEVNULL)
    except Exception as e:
        sys.stderr.write(f"[ERROR] {cmd}: {e}\n")
        sys.exit(1)


def time_spend(start: float) -> str:
    """
    Calculate and format elapsed time since a given start time.

    Parameters:
        start (float): Starting time in seconds (as returned by time.time())

    Returns:
        str: Formatted time string in HH:MM:SS format
    """
    done = time.time()
    elapsed = done - start
    return time.strftime("%H:%M:%S", time.gmtime(elapsed))


def load_acc_list(filepath: str) -> set[str]:
    """
    Load a list of accession numbers to exclude from processing.

    Reads a file containing accession numbers (one per line) and returns them
    as a set for efficient lookup during processing. Empty lines are ignored.

    Parameters:
        filepath (str): Path to the file containing accession numbers to exclude

    Returns:
        set: Set of accession numbers to exclude. Returns empty set if input file is empty.

    Example:
        exclude_list = load_acc_list('exclude.txt')
    """
    with open(filepath) as f:
        acc_list = f.read().splitlines()

    if len(acc_list) == 0:
        logging.warning(f"Exclude accession list is empty.")
        return set()
    else:
        return set(acc_list)


def load_database_stats(db_stats_file: str) -> pd.DataFrame:
    """
    Load database signature statistics from a stats file.

    Reads a tab-delimited stats file containing information about
    taxonomic signatures and their lengths.

    Parameters:
        db_stats_file (str): Path to the database stats file

    Returns:
        pd.DataFrame: df indexed with taxid, contains signature lengths and genome sizes

    Note:
        The input stats file is an 9-column tab-delimited file with:
        1. Rank
        2. Name
        3. Taxid
        4. Superkingdom
        5. NumOfSeq
        6. Max
        7. Min
        8. TotalLength
        9. GenomeSize
       10. Note (optional)
    """

    # Determine the number of columns in the stats file
    header = pd.read_csv(db_stats_file, nrows=0, sep='\t', header=None)
    valid_col_count = len(header.columns)

    usecols = [0, 2, 7, 8]
    names = ['DB_level', 'Taxid', 'TotalLength', 'Note']

    if valid_col_count == 10:
        usecols=[0, 2, 7, 8, 9]
        names=['DB_level', 'Taxid', 'TotalLength', 'GenomeSize', 'Note']

    # Set header to None to support files without headers
    df_stats = pd.read_csv(db_stats_file,
                           low_memory=False,
                           sep='\t',
                           header=None,
                           usecols=usecols,
                           names=names,
                           dtype={'DB_level': str, 'Taxid': str},
                           index_col='Taxid')

    # If 'Note' column is not present, create it with empty strings
    if not 'Note' in df_stats:
        df_stats['Note'] = ''
    if not 'GenomeSize' in df_stats:
        df_stats['GenomeSize'] = 0

    # Remove the row with index 'Taxid' if it exists
    # This is to handle the case when the stats file having the header
    if 'Taxid' in df_stats.index:
        df_stats = df_stats.drop('Taxid')

    # Make sure the format is consistent
    try:
        df_stats['TotalLength'] = pd.to_numeric(df_stats['TotalLength'], errors='raise')
    except ValueError:
        logging.error(f"Error processing stats file. Please check the format of the file.")
        sys.exit(1)

    # This is to handle the case when the stats file does not have GenomeSize column
    # In that case, the 'Note' column will be loaded as 'GenomeSize' and filled with 0
    df_stats['GenomeSize'] = pd.to_numeric(df_stats['GenomeSize'], errors='coerce')
    df_stats['GenomeSize'] = df_stats['GenomeSize'].fillna(0).astype(int)

    return df_stats


def print_message(msg: str, silent: bool, start: float, logfile: Path, errorout: int = 0):
    """
    Print and log a timestamped message.

    Writes a message to the log file and optionally to stderr. Can also
    terminate the program with an error message.

    Parameters:
        msg (str): Message to print
        silent (bool): If True, suppress output to stderr
        start (float): Start time for timestamp calculation
        logfile (Path): Path to the log file
        errorout (int): If non-zero, exit with error after printing

    Returns:
        None

    Raises:
        SystemExit: If errorout is non-zero
    """
    message = "[%s] %s\n" % (time_spend(start), msg)

    with logfile.open("a", encoding="utf-8") as f:
        f.write(message)

    if errorout:
        sys.exit(message)
    elif not silent:
        sys.stderr.write(message)


def main(args):
    """
    Main execution function for GOTTCHA2.
    """
    global argvs
    global logfile
    global begin_t
    global df_stats
    global acc_list

    argvs = parse_args(__version__, args)
    direct_ont_flag = bool(argvs.nanopore and argvs.ont_direct)
    begin_t  = time.time()
    bamfile  = Path(argvs.bam) if argvs.bam else Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.bam"
    samfile  = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.sam"
    logfile  = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.log"
    set_start_method("fork") # for default multiprocessing method
    acc_list = set()
    split_read_flag = False
    multi_part_index_flag = False
    res_df = pd.DataFrame() # aggregated restuls
    logfile_prev = ""

    logging_level = logging.WARNING

    if argvs.debug:
        logging_level = logging.DEBUG
    elif argvs.silent:
        logging_level = logging.FATAL
    elif argvs.verbose:
        logging_level = logging.INFO

    logging.basicConfig(
        level=logging_level,
        format='[%(asctime)s] [%(levelname)s] [%(module)s] %(message)s',
        datefmt='%Y%m%d %H:%M:%S',
   )

    #dependency check
    if sys.version_info < (3,9):
        sys.exit("[ERROR] Python 3.9 or above is required.")

    dependency_check("minimap2")
    dependency_check("samtools")
    if argvs.fast:
        dependency_check("sylph")

    #prepare output object
    argvs.relAbu = argvs.relAbu.upper()
    outfile_full = Path(argvs.outdir) / f"{argvs.prefix}.full.tsv"
    outfile_lineage = Path(argvs.outdir) / f"{argvs.prefix}.lineage.tsv"
    outfile_mpa = Path(argvs.outdir) / f"{argvs.prefix}.mpa.tsv"

    #create output directory if not exists
    Path(argvs.outdir).mkdir(parents=True, exist_ok=True)

    outfile = Path(argvs.outdir) / f"{argvs.prefix}.tsv"
    if argvs.format == "csv":
        outfile = Path(argvs.outdir) / f"{argvs.prefix}.csv"
    elif argvs.format == "biom":
        outfile = Path(argvs.outdir) / f"{argvs.prefix}.biom"

    out_fp = outfile.open("w", encoding="utf-8")

    if argvs.extractOnly:
        # repalce bamfile name from ".gottcha_\w+.bam" to ".log"
        logfile_prev = bamfile.with_suffix(".log")
        (mi, mf, mg, sni_argv) = (None, None, None, None)

        # if match criteria (mi/mf/mg) are not provided, load them from the log file
        if logfile_prev.is_file():
            (mi, mf, mg, sni_argv) = extract_reads.load_criteria_from_log(logfile_prev)
        else:
            logfile_prev = None

        if (argvs.sniScore is None) and sni_argv is not None:
            argvs.sniScore = sni_argv

        if (argvs.matchIdentity is None) and mi and mi >= 0:
            argvs.matchIdentity = mi
        else:
            argvs.matchIdentity = 0

        if (argvs.matchFraction is None) and mf and mf >= 0:
            argvs.matchFraction = mf
        else:
            argvs.matchFraction = 0

        if (argvs.matchLength is None) and mg and mg >= 0:
            argvs.matchLength = mg
        else:
            argvs.matchLength = 0

    (sni_score_cutoff, sni_score_species, sni_score_strain) = [float(x) for x in argvs.sniScore.split(',')]

    # display the command line
    logging.info(' '.join(sys.argv))

    print_message(f"GOTTCHA (v{__version__})", argvs.silent, begin_t, logfile)
    print_message(f"Arguments and dependencies checked:", argvs.silent, begin_t, logfile)
    print_message(f" - Database           : {argvs.database}",    argvs.silent, begin_t, logfile)
    print_message(f" - Database level     : {argvs.dbLevel}",     argvs.silent, begin_t, logfile)
    print_message(f" - Abundance          : {argvs.relAbu}",      argvs.silent, begin_t, logfile)
    print_message(f" - Output directory   : {argvs.outdir}",      argvs.silent, begin_t, logfile)
    print_message(f" - Output prefix      : {argvs.prefix}",      argvs.silent, begin_t, logfile)
    print_message(f" - Threads            : {argvs.threads}",     argvs.silent, begin_t, logfile)
    print_message(f" - Fast mode          : {argvs.fast}",        argvs.silent, begin_t, logfile)
    if argvs.input:
        print_message(f" - Input Reads        : {argvs.input}",     argvs.silent, begin_t, logfile)
    if argvs.bam:
        print_message(f" - Input BAM File     : {bamfile}",           argvs.silent, begin_t, logfile)
    if argvs.nanopore:
        print_message(f" - Nanopore Mode      : Enabled",              argvs.silent, begin_t, logfile)
    if argvs.errorRate:
        print_message(f" - Read Error Rate    : {argvs.errorRate}", argvs.silent, begin_t, logfile)
    if argvs.sigList:
        print_message(f" - Sig-of-Int List    : {argvs.sigList}", argvs.silent, begin_t, logfile)
    if argvs.sigList:
        print_message(f" - Sig-of-Int Action  : {argvs.sigListAction}", argvs.silent, begin_t, logfile)
    if argvs.minCov > 0:
        print_message(f" - Minimal SIG Cov    : {argvs.minCov}",      argvs.silent, begin_t, logfile)
    if argvs.minLen > 0:
        print_message(f" - Minimal SIG Length : {argvs.minLen}",      argvs.silent, begin_t, logfile)
    if argvs.minReads > 0:
        print_message(f" - Minimal Reads      : {argvs.minReads}",    argvs.silent, begin_t, logfile)
    if argvs.extract:
        print_message(f" - Extract Taxa       : {argvs.extract}",     argvs.silent, begin_t, logfile)
    if argvs.extractOnly:
        print_message(f" - Extract Only       : {argvs.extractOnly}", argvs.silent, begin_t, logfile)
    if argvs.maxZscore > 0:
        print_message(f" - Maximal zScore     : {argvs.maxZscore}",   argvs.silent, begin_t, logfile)
    if logfile_prev:
        print_message(f" - Load criteria from : {logfile_prev}",      argvs.silent, begin_t, logfile)
    if argvs.matchIdentity != None:
        print_message(f" - Min Match Identity : {argvs.matchIdentity}", argvs.silent, begin_t, logfile)
    if argvs.matchFraction != None:
        print_message(f" - Min Match Fraction : {argvs.matchFraction}", argvs.silent, begin_t, logfile)
    if argvs.matchLength != None:
        print_message(f" - Min Match Length   : {argvs.matchLength}", argvs.silent, begin_t, logfile)
    if argvs.sniScore != None:
        print_message(f" - SNI-score (g,s,n)  : {argvs.sniScore}",    argvs.silent, begin_t, logfile)

    #load taxonomy for taxonomic aggregation and annotation
    if not argvs.extractOnly:
        print_message("Loading taxonomy information...", argvs.silent, begin_t, logfile)

        if Path(argvs.database + ".tax.tsv").exists():
            custom_taxa_tsv = Path(argvs.database + ".tax.tsv")
        elif Path(argvs.database + ".taxa").exists():
            custom_taxa_tsv = Path(argvs.database + ".taxa")

        taxonomy.loadTaxonomy(cus_taxonomy_file=custom_taxa_tsv, auto_download=False)
        print_message(f" - {len(taxonomy.taxNames):,} taxa loaded.", argvs.silent, begin_t, logfile)

        #load database stats
        print_message("Loading database stats...", argvs.silent, begin_t, logfile)
        if Path(argvs.database + ".stats").is_file():
            df_stats = load_database_stats(argvs.database+".stats")
        else:
            print_message(f"ERROR: {argvs.database+'.stats'} not found.", argvs.silent, begin_t, logfile, errorout=1)

        print_message(f" - {df_stats.shape[0]:,} entries loaded.", argvs.silent, begin_t, logfile)
        print_message(f" - signatures at {df_stats['DB_level'].unique().tolist()} levels loaded.", argvs.silent, begin_t, logfile)

    if argvs.sigList:
        print_message("Loading accession#s of interest list...", argvs.silent, begin_t, logfile)
        acc_list = load_acc_list(argvs.sigList)
        print_message(f" - {len(acc_list):,} accession/signature of interest loaded.", argvs.silent, begin_t, logfile)

    # Summary of the Main Process:
    #
    # Input Reads
    #     ↓
    # [Nanopore Preprocessing] (optional)
    #     ↓
    # [Run fast query] (sylph; optionall; if fast mode is on)
    #     ↓
    # [Extract queried signatures] (optionall; if fast mode is on)
    #     ↓
    # Read Mapping (minimap2)
    #     ↓
    # Alignments (SAM File)
    #     ↓
    # [Remove Multiple Hits] (if multi-part index)
    #     ↓
    # [Remove Inconsistent Chunks] (if nanopore)
    #     ↓
    # BAM Conversion + Indexing
    #     ↓
    # Parse & Filter Alignments
    #     ↓
    # Group to Strains
    #     ↓
    # Aggregate Taxonomy
    #     ↓
    # Generate Reports
    #     ↓
    # [Extract Reads] (optional)

    if argvs.input:
        # if fast mode is on, run Sylph query to prefilter the reference genomes and create a smaller reference for read mapping; 
        # otherwise, use the full database index for read mapping
        minimap2_index = "" if argvs.fast else f"{argvs.database}.mmi"
        
        # The original input reads for Sylph sketch and query; Not the split reads for minimap2 mapping if nanopore option is on
        sylph_input = argvs.input

        # if nanopore option is on, preprocessing reads
        if argvs.nanopore:
            print_message("Checking nanopore read files...", argvs.silent, begin_t, logfile)
            if direct_ont_flag:
                print_message(" - Direct ONT mode: mapping intact reads", argvs.silent, begin_t, logfile)
            else:
                argvs.input = ont_utils.preprocess_nanopore_reads(argvs.input, argvs.outdir, argvs.prefix, argvs.silent)
                split_read_flag = True

        if argvs.fast:
            print_message("Prefiltering reference genomes...", argvs.silent, begin_t, logfile)
            sylph_db = f"{argvs.database}.syldb"
            g2_archive = f"{argvs.database}.zip"
            sylph_query_tsv = Path(argvs.outdir) / f"{argvs.prefix}.sylph_query.tsv"
            queried_signatures_file = Path(argvs.outdir) / f"{argvs.prefix}.sylph_queried_signatures.txt"
            extracted_reference = Path(argvs.outdir) / f"{argvs.prefix}.sylph_extracted.fa.gz"
            
            # extract subsample (cXXX) rate from sylph_db string, default set to 100
            subsampling_rate = 100
            match = re.search(r'c(\d+)\.', sylph_db)
            if match:
                subsampling_rate = int(match.group(1))

            fast_min_kmer = argvs.fast_min_kmer
            fast_min_ani = argvs.fast_min_ani

            # Run Sylph sketch if the input file is in FASTA format
            if Path(sylph_input[0]).name.endswith(('.fa', '.fasta', '.fa.gz', '.fna', '.fna.gz', '.fasta.gz')):
                print_message("Generating sketchs for FASTA input reads...", argvs.silent, begin_t, logfile)
                try:
                    sylph_result = prefilter.run_sylph_sketch(
                        read_file=sylph_input[0],
                        outdir=str(argvs.outdir),
                        threads=argvs.threads,
                        subsampling_rate=subsampling_rate,
                    )
                except (FileNotFoundError, subprocess.CalledProcessError) as e:
                    print_message(f"ERROR: prefiltering failed: {e}", argvs.silent, begin_t, logfile, errorout=1)

                with logfile.open("a", encoding="utf-8") as f:
                    if sylph_result.stdout:
                        f.write(sylph_result.stdout)
                    if sylph_result.stderr:
                        f.write(sylph_result.stderr)
                
                sylph_input = [str(Path(argvs.outdir) / f"{Path(sylph_input[0]).name}.sylsp")]

            # Run Sylph query to get the list of signatures that are likely present in the input reads
            try:
                sylph_result = prefilter.run_sylph_query(
                    database=sylph_db,
                    reads=sylph_input,
                    output=str(sylph_query_tsv),
                    threads=argvs.threads,
                    subsampling_rate=subsampling_rate,
                    minimum_kmer=fast_min_kmer,
                    minimum_ani=fast_min_ani,
                    read_seq_id=float(100-argvs.errorRate*100)
                )
            except (FileNotFoundError, subprocess.CalledProcessError) as e:
                print_message(f"ERROR: prefiltering failed: {e}", argvs.silent, begin_t, logfile, errorout=1)

            with logfile.open("a", encoding="utf-8") as f:
                if sylph_result.stdout:
                    f.write(sylph_result.stdout)
                if sylph_result.stderr:
                    f.write(sylph_result.stderr)

            try:
                pd.read_csv(sylph_query_tsv, 
                            sep='\t', 
                            usecols=['Genome_file'], 
                            dtype={'Genome_file': str})['Genome_file'].dropna().str.strip().to_csv(queried_signatures_file, index=False, header=False)
            except pd.errors.EmptyDataError:
                queried_signatures = pd.Series(dtype=str)
            except (FileNotFoundError, ValueError) as e:
                print_message(f"ERROR: unable to parse Sylph query output {sylph_query_tsv}: {e}", argvs.silent, begin_t, logfile, errorout=1)

            filenames = sig_archive.read_file_list(queried_signatures_file, filename_only=True)
            print_message(f" - Identified {len(filenames):,} reference genomes.", argvs.silent, begin_t, logfile)
            
            if len(filenames) == 0:
                print_message("No references identified. GOTTCHA2 stopped.", argvs.silent, begin_t, logfile)
                sys.exit(0)

            # Extract those signatures from the archive to create a smaller reference for read mapping
            extracted_content, processed_files, skipped_files = sig_archive.quick_concat(g2_archive,
                                                                                         separator=str('\n').encode('utf-8'),
                                                                                         skip_missing=False, 
                                                                                         filenames=filenames)

            extracted_reference.write_bytes(extracted_content)
            if extracted_reference.stat().st_size == 0:
                print_message(
                    f"ERROR: no queried signatures could be extracted from {sylph_db}; {extracted_reference} is empty.",
                    argvs.silent,
                    begin_t,
                    logfile,
                    errorout=1
                )

            if skipped_files:
                print_message(f" - {len(skipped_files):,} queried signatures were not found in the archive.", argvs.silent, begin_t, logfile)
            print_message(f" - {len(processed_files):,} reference genomes extracted.", argvs.silent, begin_t, logfile)
            minimap2_index = str(extracted_reference)

        print_message("Running read-mapping...", argvs.silent, begin_t, logfile)
        exitcode, cmd, input_read_count, multi_part_index_flag = read_mapping.minimap2(
            argvs.input,
            minimap2_index,
            argvs.threads,
            argvs.m2options,
            argvs.presetx,
            samfile,
            logfile,
            allow_secondary=direct_ont_flag,
            max_secondary=argvs.ont_max_secondary,
            secondary_ratio=argvs.ont_secondary_ratio,
        )
        logging.info(f"COMMAND: {cmd}")

        if exitcode != 0:
            # if size of the samfile is zero
            print_message(f"Logfile saved to {logfile}.", argvs.silent, begin_t, logfile)
            sys.exit("[%s] ERROR: error occurred while running read mapping (exit: %s).\n" % (time_spend(begin_t), exitcode))
        else:
            print_message(f" - {input_read_count:,} input reads processed.", argvs.silent, begin_t, logfile)
            print_message(f"Mapped SAM file saved to {samfile}.", argvs.silent, begin_t, logfile)
        gc.collect()

    # remove multiple hits
    if multi_part_index_flag and not direct_ont_flag:
        # remove multiple hits from the SAM file
        print_message("Removing multiple hits from SAM file...", argvs.silent, begin_t, logfile)
        samfile_temp = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.sam.temp"
        flag, aln_count, top_hits_count = read_mapping.post_processing_sam(samfile, samfile_temp)
        if flag:
            samfile_temp.rename(samfile)
            # Note:
            # When input of the gottcha2 is a SAM file and new outdir/prefix is provided, the output will be saved to that location.
            # If not, the output will overwrite the original SAM file.
            print_message(f" - {aln_count:,} total alignments", argvs.silent, begin_t, logfile)
            print_message(f" - {top_hits_count:,} best hits among index-partitions", argvs.silent, begin_t, logfile)

        gc.collect()

    # preprocess SAM file for nanopore reads
    if direct_ont_flag and Path(samfile).is_file():
        print_message("Resolving direct ONT alignments...", argvs.silent, begin_t, logfile)
        samfile_temp = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.sam.temp"
        tol_alignment_cnt, tol_q_alignment_cnt = ont_utils.direct_ont_reads_samfile_postprocessing(samfile, samfile_temp, argvs.ont_min_species_support)
        samfile_temp.replace(samfile)
        print_message(f" - {tol_alignment_cnt:,} total alignments", argvs.silent, begin_t, logfile)
        print_message(f" - {tol_q_alignment_cnt:,} qualified-species alignments retained", argvs.silent, begin_t, logfile)
        if tol_q_alignment_cnt == 0:
            print_message("No direct ONT alignments remained after species resolution. Stopping.", argvs.silent, begin_t, logfile)
            sys.exit(0)
        gc.collect()
    elif argvs.nanopore and Path(samfile).is_file():
        print_message("Removing inconsistent read chunks from SAM file...", argvs.silent, begin_t, logfile)
        samfile_temp = Path(argvs.outdir) / f"{argvs.prefix}.gottcha_{argvs.dbLevel}.sam.temp"
        tol_chunks_count, tol_chunks_qualified = ont_utils.split_reads_samfile_postprocessing(samfile, samfile_temp)
        if tol_chunks_count > 0:
            samfile_temp.rename(samfile)
        print_message(f" - {tol_chunks_count:,} mapped read chunks processed", argvs.silent, begin_t, logfile)
        print_message(f" - {tol_chunks_count-tol_chunks_qualified:,} inconsistent hits removed", argvs.silent, begin_t, logfile)
        gc.collect()

    # processing alignments and generate results
    if not argvs.extractOnly:
        if os.path.isfile(os.path.abspath(samfile)):
            print_message("Converting to BAM file...", argvs.silent, begin_t, logfile)
            sam_to_bam.convert_sam_to_bam(input_sam=os.path.abspath(samfile),
                                          output_bam=os.path.abspath(bamfile),
                                          threads=argvs.threads,
                                          quiet=argvs.silent)
            print_message(f"BAM file saved to {bamfile}...", argvs.silent, begin_t, logfile)
            
            file_path = Path(samfile)
            if file_path.exists():
                file_path.unlink()
            gc.collect()

        if Path(bamfile).exists() and Path(f"{bamfile}.bai").exists():
            print_message("Processing alignments...", argvs.silent, begin_t, logfile)
            ref_chunk_results = process_bam.parse_aln_from_bam(
                bam_path=bamfile,
                processes=argvs.threads,
                min_frac=argvs.matchFraction,
                min_idt=argvs.matchIdentity,
                min_alen=argvs.matchLength,
                include_secondary=direct_ont_flag,
                include_supplementary=direct_ont_flag,
                split_read_flag=split_read_flag
            )

            str_df, soi_read_count = aggregate_results.group_refs_to_strains(ref_chunk_results, acc_list, argvs.sigListAction, df_stats)

            tol_read_count = str_df['READ_COUNT'].sum()
            tol_invalid_match_count = str_df['INVALID_ALNS'].sum()

            print_message(f" - {tol_invalid_match_count:,} alignments did not meet matching criteria", argvs.silent, begin_t, logfile)
            print_message(f" - {tol_read_count:,} qualified reads processed", argvs.silent, begin_t, logfile)

            if not tol_read_count:
                print_message("No qualified alignments found. Stopping.", argvs.silent, begin_t, logfile)
                sys.exit(0)

            # aggregate the results
            _args = (str_df,
                     argvs.relAbu,
                     argvs.dbLevel,
                     argvs.minCov,
                     argvs.minReads,
                     argvs.minLen,
                     argvs.maxZscore,
                     sni_score_species,
                     sni_score_strain,
                     sni_score_cutoff,
                     argvs.errorRate)
            res_df, soi_read_count = aggregate_results.aggregate_taxonomy(*_args)

            if acc_list:
                print_message(f" - {soi_read_count:,} reads mapped to accession-of-interest", argvs.silent, begin_t, logfile)
                read_count_after_soi = tol_read_count
                if argvs.sigListAction == 'filter_out':
                    read_count_after_soi = tol_read_count - soi_read_count
                elif argvs.sigListAction == 'filter_in':
                    read_count_after_soi = soi_read_count
                print_message(f" - {read_count_after_soi:,} reads after applying accession-of-interest action ({argvs.sigListAction})", argvs.silent, begin_t, logfile)

            print_message("Done taxonomy aggregation.", argvs.silent, begin_t, logfile)

            if not len(res_df):
                print_message("No qualified taxonomy profiled.", argvs.silent, begin_t, logfile)
            else:
                # generate output results
                if argvs.format == "biom":
                    report.generate_biom_file(res_df, out_fp, argvs.dbLevel, argvs.prefix)
                else:
                    report.generate_report_file(res_df, out_fp, outfile_full, argvs.format)
                # generate lineage file
                target_idx = (res_df['LEVEL']==argvs.dbLevel) & \
                                (res_df['NOTE'].str.contains('Filtered out', na=False) == False) & \
                                (res_df['NOTE'].str.contains('Not shown', na=False) == False)
                target_df = res_df.loc[target_idx, ['ABUNDANCE','TAXID']]
                tax_num = len(target_df)

                print_message(f"{tax_num} qualified {argvs.dbLevel} profiled.", argvs.silent, begin_t, logfile)

                if tax_num:
                    report.generate_lineage_file(target_df, outfile_lineage)

                    if argvs.mpa:
                        target_df = res_df.loc[target_idx, ['TAXID', 'REL_ABUNDANCE', 'REL_ABUNDANCE_GC','READ_COUNT', 'SIG_COV']]
                        report.generate_mpa_file(target_df, outfile_mpa)
                        print_message(f"MPA format file saved to {outfile_mpa}.", argvs.silent, begin_t, logfile)

                print_message(f"Results saved to {outfile}.", argvs.silent, begin_t, logfile)
        else:
            print_message(f"ERROR: BAM file {bamfile} or its index not found.", argvs.silent, begin_t, logfile, errorout=1)
            print_message("GOTTCHA2 stopped.", argvs.silent, begin_t, logfile)
            sys.exit(0)

    # extracting reads
    if argvs.extract:
        (taxa_arg, max_per_taxon, out_format) = (argvs.extract.split(':', maxsplit=2) + ['all', 'fasta'])[:3]

        print_message(f"Extracting {max_per_taxon} sequences per taxa in {out_format} format...", argvs.silent, begin_t, logfile)

        full_report_file = ""

        if max_per_taxon.isdigit() or max_per_taxon == 'all':
            max_per_taxon = int(max_per_taxon) if max_per_taxon != 'all' else 0

        if argvs.extractOnly:
            # repalce bamfile name to replace from ".gottcha_\w+.bam" to ".full.tsv"
            full_report_file = re.sub(r"\.gottcha_\w+\.bam$", ".full.tsv", str(bamfile))

        taxa_dict, ref_to_extract_taxid = extract_reads.parse_taxids(taxa_arg, 
                                                                     res_df, 
                                                                     full_report_file, 
                                                                     sni_score_cutoff,
                                                                     sni_score_species, 
                                                                     sni_score_strain
                                                                     )

        if not len(ref_to_extract_taxid):
            print_message("No qualified taxonomy profiled.", argvs.silent, begin_t, logfile)

        outfile = Path(argvs.outdir) / f"{argvs.prefix}.extract.{out_format.lower()}"

        _args = (os.path.abspath(bamfile),
                    taxa_dict,
                    ref_to_extract_taxid,
                    outfile.open("w", encoding="utf-8"),
                    argvs.threads,
                    argvs.matchFraction,
                    argvs.matchIdentity,
                    argvs.matchLength,
                    max_per_taxon,
                    acc_list,
                    argvs.sigListAction,
                    out_format)
        taxon_count, seq_count = extract_reads.extract_sequences_by_taxonomy(*_args)
        print_message(f"Done extracting {seq_count} sequences from {taxon_count} taxa to '{outfile}'.",
                        argvs.silent, begin_t, logfile)

if __name__ == '__main__':
    main(sys.argv[1:])
