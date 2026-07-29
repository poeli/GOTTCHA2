#!/usr/bin/env python3

import pandas as pd
import numpy as np
import pysam
import argparse
import sys
import os
import subprocess
import tempfile
import logging
from typing import Dict, Set
import gc

def convert_sam_to_bam(input_sam: str, output_bam: str, threads=4, quiet=False) -> None:
    """Convert SAM to sorted BAM with proper headers - optimized for large files."""
    
    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Step 1: Create headers file from SAM (with streaming)
        header_file = os.path.join(temp_dir, "header.sam")
        
        logging.info(f"Extracting references from {input_sam}...")

        df = pd.read_csv(input_sam,
                         sep='\t',
                         header=None,
                         usecols=[2, 6],
                         names=['REF', 'MATE_REF'],
                         dtype='str')

        logging.debug(f"Sample of REF/MATE_REF columns:\n{df.head()}")

        # Extract unique references from both REF and MATE_REF columns, excluding '*' and '='
        idx = (df['MATE_REF'] != '=') & (df['MATE_REF'] != '*')
        refs = pd.concat([df['REF'], df.loc[idx, 'MATE_REF']]).unique()

        logging.debug(f"Unique references extracted: {len(refs)}")

        logging.info("Creating header file...")
        df = pd.DataFrame({'REF': refs})

        def convert_to_header(row):
            ref = row['REF']
            start, end = ref.split('|')[1:3]
            return f"@SQ\tSN:{ref}\tLN:{int(end)-int(start)+1}"
        
        with open(header_file, 'w') as out:
            out.write("@HD\tVN:1.0\tSO:coordinate\n")
            df.apply(convert_to_header, axis=1).to_csv(out, index=False, header=False)

        # open header_file and print the first 2s lines for debugging
        with open(header_file) as f:
            header_lines = [next(f) for _ in range(2)]
            logging.debug(f"Sample of header file:\n{''.join(header_lines)}")

        del df  # Free memory
        gc.collect()

        # Step 2: Create temp BAM with header
        temp_bam = os.path.join(temp_dir, "temp.bam")
        logging.info("Converting to BAM...")
        
        # Use samtools directly for better performance
        cmd = f"cat {header_file} {input_sam} | samtools view -b -@ {threads} -o {temp_bam}"
        if quiet:
            cmd += " 2>/dev/null"
        logging.debug(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        
        # Step 3: Sort BAM (with threads)
        logging.info("Sorting BAM file...")
        cmd = f"samtools sort -@ {threads} -o {output_bam} {temp_bam}"
        if quiet:
            cmd += " 2>/dev/null"
        logging.debug(f"Running command: {cmd}")
        subprocess.run(cmd, shell=True, check=True)
        
        # Step 4: Index BAM
        logging.info("Creating BAM index...")
        try:
            if quiet:
                subprocess.run(f"samtools index {output_bam} 2>/dev/null", shell=True, check=True)
            else:
                pysam.index(output_bam)
        except Exception as e:
            logging.warning(f"Warning: Could not create index: {e}")
        
    logging.info(f"Conversion complete: {output_bam}")

def main(argv):
    parser = argparse.ArgumentParser(
        prog="gottcha2 sam2bam",
        description="Convert SAM to sorted BAM - optimized for large files"
    )
    parser.add_argument('-i', '--input', required=True,
                       help="Input SAM file")
    parser.add_argument('-o', '--output', required=True,
                       help="Output sorted BAM file")
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help="Number of threads for processing (default: 4)")
    parser.add_argument('-q', '--quiet', action='store_true',
                       help="Suppress warning messages")
    
    args = parser.parse_args(argv)
    
    try:
        convert_sam_to_bam(args.input, args.output, args.threads, args.quiet)
        logging.info("SAM to BAM conversion successful.")
    except Exception as e:
        logging.error(f"Error converting file: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])
