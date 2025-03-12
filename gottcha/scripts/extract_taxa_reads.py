#!/usr/bin/env python3

import argparse
import gzip
import os
import sys
import subprocess
import tempfile
import logging
from pathlib import Path

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description="Extract reads from SAM/BAM files based on reference patterns.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('-i', '--input', required=True,
                        help="Input SAM or BAM file")
    
    parser.add_argument('-o', '--output', 
                        help="Output file prefix (default: extracted_reads)")
    
    # Reference options
    parser.add_argument('-t', '--taxa', 
                        help="Comma-separated list of taxa to extract")
    
    parser.add_argument('-f', '--taxa-file',
                        help="File containing taxa (one per line)")
    
    # Shared options
    parser.add_argument('--threads', type=int, default=1,
                        help="Number of threads for parallel processing")
    
    parser.add_argument('--verbose', action='store_true',
                        help="Enable verbose output")
    
    return parser.parse_args()

def setup_logging(verbose=False):
    """Configure logging based on verbosity level."""
    level = logging.INFO if verbose else logging.WARNING
    
    logging.basicConfig(
        level=level,
        format='%(asctime)s [%(levelname)s] %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def parse_taxa(taxa_arg=None, taxa_file=None):
    """Parse taxa from command line argument or file."""
    if taxa_arg:
        return [x.strip() for x in taxa_arg.split(',') if x.strip()]
    elif taxa_file:
        try:
            with open(taxa_file) as f:
                return [x.strip() for x in f.readlines() if x.strip()]
        except IOError as e:
            sys.stderr.write(f"Error reading taxa file {taxa_file}: {e}\n")
            sys.exit(1)
    return []

def extract_matching_reads(input_file, taxa_list, output_prefix, threads=1):
    """Extract reads containing specified taxa and output to gzipped FASTA and tab-delimited file."""
    # Check if input is BAM or SAM
    is_bam = input_file.lower().endswith('.bam')
    
    # Output filenames
    fasta_output = f"{output_prefix}.fasta.gz"  # Changed to .fasta.gz
    tab_output = f"{output_prefix}.tsv"
    
    # Prepare grep pattern for extraction
    taxa_pattern = '|'.join(taxa_list)
    grep_pattern = f'({taxa_pattern})'
    
    with tempfile.TemporaryDirectory() as temp_dir:
        # Create temp file for extracted SAM rows
        extracted_sam = os.path.join(temp_dir, "extracted_rows.sam")
        
        # Extract matching rows using grep
        if is_bam:
            logging.info(f"Extracting rows from BAM file: {input_file}")
            cmd = f"samtools view -h {input_file} | grep -E '{grep_pattern}' > {extracted_sam}"
        else:
            logging.info(f"Extracting rows from SAM file: {input_file}")
            cmd = f"grep -E '{grep_pattern}' {input_file} > {extracted_sam}"
        
        # Run extraction command
        subprocess.run(cmd, shell=True, check=True)
        
        # Process the extracted rows
        read_count = 0
        with open(extracted_sam, 'r') as sam, \
             gzip.open(fasta_output, 'wt') as fasta, \
             open(tab_output, 'w') as tab:
            
            for line in sam:
                if line.startswith('@'):  # Skip header lines
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 11:  # Not enough fields for a valid SAM record
                    continue
                
                # Extract required fields
                read_id = fields[0]
                ref_name = fields[2]
                flag = int(fields[1])
                cigar = fields[5]
                seq = fields[9]
                
                # Create FASTA header
                fasta_header = f"@{read_id} {ref_name} {cigar}"
                
                # Write to gzipped FASTA file
                fasta.write(f"{fasta_header}\n{seq}\n")
                
                # Write to tab-delimited file (uncompressed)
                tab.write(f"{fasta_header}\t{seq}\n")
                
                read_count += 1
    
    logging.info(f"Extracted {read_count} reads to {fasta_output} (gzipped) and {tab_output}")
    return read_count

def main():
    """Main function to run the extraction process."""
    args = parse_args()
    setup_logging(args.verbose)
    
    # Parse taxa
    taxa_list = parse_taxa(args.taxa, args.taxa_file)
    
    if not taxa_list:
        sys.exit("ERROR: You must specify taxa to extract with --taxa or --taxa-file")
    
    # Check if input file exists
    if not os.path.exists(args.input):
        sys.exit(f"ERROR: Input file {args.input} does not exist")
    
    # Generate output filename prefix if not specified
    output_prefix = args.output if args.output else "extracted_reads"
    
    # Extract reads
    logging.info(f"Starting extraction for {len(taxa_list)} taxa")
    extract_matching_reads(args.input, taxa_list, output_prefix, args.threads)
    logging.info("Extraction complete")

if __name__ == "__main__":
    main()