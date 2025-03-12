#!/usr/bin/env python3

import pysam
import argparse
import sys
import os
import subprocess
import tempfile
from typing import Dict, Set

def convert_sam_to_bam(input_sam: str, output_bam: str, threads=4, quiet=False) -> None:
    """Convert SAM to sorted BAM with proper headers - optimized for large files."""
    
    # Create a temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as temp_dir:
        # Step 1: Create headers file from SAM (with streaming)
        refs_file = os.path.join(temp_dir, "refs.txt")
        header_file = os.path.join(temp_dir, "header.sam")
        
        print(f"Extracting references from {input_sam}...")
        with open(refs_file, 'w') as out:
            # Stream through file with grep/awk for better performance
            # Extract both reference (column 3) and mate reference (column 7)
            cmd = f"awk '!(/^@/) {{ print $3; if ($7 != \"=\") print $7 }}' {input_sam} | sort -u"
            subprocess.run(cmd, shell=True, stdout=out, check=True)
            
        print("Creating header file...")
        with open(header_file, 'w') as out:
            out.write("@HD\tVN:1.0\tSO:coordinate\n")
            with open(refs_file) as f:
                for line in f:
                    ref = line.strip()
                    if ref and ref != '*':
                        start, end = ref.split('|')[1:3]
                        out.write(f"@SQ\tSN:{ref}\tLN:{end-start+1}\n")

        # Step 2: Create temp BAM with header
        temp_bam = os.path.join(temp_dir, "temp.bam")
        print("Converting to BAM...")
        
        # Use samtools directly for better performance
        cmd = f"cat {header_file} {input_sam} | samtools view -b -@ {threads} -o {temp_bam}"
        if quiet:
            cmd += " 2>/dev/null"
        subprocess.run(cmd, shell=True, check=True)
        
        # Step 3: Sort BAM (with threads)
        print("Sorting BAM file...")
        cmd = f"samtools sort -@ {threads} -o {output_bam} {temp_bam}"
        if quiet:
            cmd += " 2>/dev/null"
        subprocess.run(cmd, shell=True, check=True)
        
        # Step 4: Index BAM
        print("Creating BAM index...")
        try:
            if quiet:
                subprocess.run(f"samtools index {output_bam} 2>/dev/null", shell=True, check=True)
            else:
                pysam.index(output_bam)
        except Exception as e:
            print(f"Warning: Could not create index: {e}", file=sys.stderr)
        
    print(f"Conversion complete: {output_bam}")

def main():
    parser = argparse.ArgumentParser(description="Convert SAM to sorted BAM - optimized for large files")
    parser.add_argument('-i', '--input', required=True,
                       help="Input SAM file")
    parser.add_argument('-o', '--output', required=True,
                       help="Output sorted BAM file")
    parser.add_argument('-t', '--threads', type=int, default=4,
                       help="Number of threads for processing (default: 4)")
    parser.add_argument('-q', '--quiet', action='store_true',
                       help="Suppress warning messages")
    
    args = parser.parse_args()
    
    try:
        convert_sam_to_bam(args.input, args.output, args.threads, args.quiet)
    except Exception as e:
        print(f"Error converting file: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
