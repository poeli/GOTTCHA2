#!/usr/bin/env python3
import argparse
import gzip
import itertools
import numpy as np
import sys
import pandas as pd
import os
import logging

try:
    from . import taxonomy as t
except ImportError:
    import gottcha.utils.taxonomy as t

def open_in(path: str):
    if path == "-":
        return sys.stdin.buffer
    if path.endswith(".gz"):
        return gzip.open(path, "rb")
    return open(path, "rb", buffering=5 * 1024 * 1024)

def open_out(path: str, gz_level: int):
    if path == "-":
        return sys.stdout.buffer
    if path.endswith(".gz"):
        # compresslevel=1 is much faster than default(9) for big outputs
        return gzip.open(path, "wb", compresslevel=gz_level)
    return open(path, "wb", buffering=5 * 1024 * 1024)

def detect_format(first_line: bytes) -> str:
    if not first_line:
        raise ValueError("Empty input")
    b0 = first_line[:1]
    if b0 == b">":
        return "fasta"
    if b0 == b"@":
        return "fastq"
    raise ValueError("Unknown format: expected FASTA (>) or FASTQ (@)")

def write_chunks_fasta(out, rid: bytes, seq: bytes, L: int, step: int, drop_tail: bool, min_tail: int, prefix: bytes):
    n = len(seq)
    if n == 0:
        return

    out_lines = []
    chunk_id = 0
    chunks_written = 0
    # full chunks
    last_full_start = n - L
    if last_full_start >= 0:
        for s in range(0, last_full_start + 1, step):
            e = s + L
            # >id|chunk=k
            header = b">" + prefix + rid + b"|chunk=" + str(chunk_id).encode() + b"\n"
            out_lines.append(header)
            out_lines.append(seq[s:e] + b"\n")
            chunk_id += 1
            chunks_written += 1

        # tail
        tail_start = ((last_full_start // step) + 1) * step
        if not drop_tail and tail_start < n:
            tail_len = n - tail_start
            if tail_len >= min_tail:
                header = b">" + prefix + rid + b"|chunk=" + str(chunk_id).encode() + b"\n"
                out_lines.append(header)
                out_lines.append(seq[tail_start:] + b"\n")
                chunk_id += 1
                chunks_written += 1
    else:
        # read shorter than L
        if not drop_tail and n >= min_tail:
            header = b">" + prefix + rid + b"|chunk=" + str(chunk_id).encode() + b"\n"
            out_lines.append(header)
            out_lines.append(seq + b"\n")
            chunk_id += 1
            chunks_written += 1

    out.writelines(out_lines)
    return chunks_written

def split_to_fasta(input_path: str, 
                   output_path: str, 
                   split_length: int, 
                   step_length: int,
                   drop_tail: bool, 
                   min_tail: int = 1,
                   prefix: str = "", 
                   gzip_level: int = 1) -> int:
    """
    Split long reads from input_path into fixed-length chunks and write to output_path.

    Returns the total number of chunks written.
    """
    prefix_b = prefix.encode()
    chunk_count = 0

    with open_in(input_path) as hin, open_out(output_path, gzip_level) as hout:
        first = hin.readline()
        if not first:
            return 0
        fmt = detect_format(first)

        if fmt == "fastq":
            # FASTQ assumed standard 4-line records (seq on one line).
            header = first
            while header:
                if not header.startswith(b"@"):
                    raise ValueError("Malformed FASTQ: expected '@' header")
                seq = hin.readline()
                plus = hin.readline()
                qual = hin.readline()
                if not seq or not plus or not qual:
                    break
                # read id up to first whitespace
                rid = header[1:].strip().split(None, 1)[0]
                seq = seq.strip()
                chunk_count += write_chunks_fasta(hout, rid, seq, split_length, step_length, drop_tail, min_tail, prefix_b)
                header = hin.readline()

        else:
            # FASTA: collect per record (still fast; avoids per-base shifting).
            header = first
            rid = None
            seq_parts = []
            while True:
                if header.startswith(b">"):
                    if rid is not None:
                        seq = b"".join(seq_parts)
                        chunk_count += write_chunks_fasta(hout, rid, seq, split_length, step_length, drop_tail, min_tail, prefix_b)
                    rid = header[1:].strip().split(None, 1)[0]
                    seq_parts = []
                else:
                    seq_parts.append(header.strip())

                header = hin.readline()
                if not header:
                    break

            if rid is not None:
                seq = b"".join(seq_parts)
                chunk_count += write_chunks_fasta(hout, rid, seq, split_length, step_length, drop_tail, min_tail, prefix_b)

    return chunk_count


def preprocess_nanopore_reads(reads, outdir, prefix, silent):
    """
    Split nanopore reads into shorter chunks before mapping.

    Parameters:
        reads (list): List of input read file objects (expects exactly one)
        outdir (str): Output directory for temporary files
        prefix (str): Prefix for the generated chunk file
        silent (bool): Silence flag for logging

    Returns:
        list: Updated list of read objects pointing to the chunked read file
    """
    os.makedirs(outdir, exist_ok=True)

    if len(reads) != 1:
        logging.fatal("ERROR: Nanopore read processing expects a single input file.")
        sys.exit(1)

    input_path = reads[0]
    output_path = os.path.join(outdir, f"{prefix}.split_reads.fasta.gz")

    try:
        chunk_count = split_to_fasta(input_path, output_path, split_length=150, step_length=150, drop_tail=True)
    except Exception as e:
        logging.fatal(f"ERROR: Failed to split nanopore reads: {e}")
        sys.exit(1)
    else:
        if chunk_count == 0:
            logging.fatal("ERROR: No reads were produced after splitting nanopore reads.")
            sys.exit(1)
        logging.info(f" - {chunk_count} chunks written to {output_path}")

    return [output_path]


def direct_ont_reads_samfile_postprocessing(
        samfile,
        samfile_temp,
        ont_min_species_ratio):

    df = pd.read_csv(
        samfile,
        sep="\t",
        header=None,
        usecols=[0, 2, 13],
        names=["QNAME", "REF", "AS"],
        converters={
            'AS': lambda x: x.replace('AS:i:', '')
        }
    )

    df[['AS']] = df[['AS']].astype('int16')

    logging.info(f"{len(df):,} alignments loaded from the SAM file.")

    taxids = df["REF"].str.rsplit("|", n=2).str[-2]
    taxid_to_species = {
        taxid: t.taxid2taxidOnRank(
            taxid,
            target_rank="species"
        )
        for taxid in taxids.unique()
    }

    df["SPECIES"] = taxids.map(taxid_to_species)

    logging.info(f"{df['SPECIES'].nunique():,} unique species taxids found in the SAM file.")

    # Total aligned length for each species within each read
    species_as = (
        df.groupby(["QNAME", "SPECIES"], sort=False, dropna=False)["AS"].transform("sum").to_numpy()
    )

    # Total aligned length across all species for each read
    read_as = (
        df.groupby("QNAME", sort=False)["AS"].transform("sum").to_numpy()
    )

    # Keep every alignment belonging to a species satisfying:
    #        species_AS / total_read_AS > ont_min_species_ratio
    qualified_mask = (
        species_as >= read_as * ont_min_species_ratio
    )

    n_qualified = int(qualified_mask.sum())
    n_total = len(qualified_mask)

    logging.info(
        f"{n_qualified:,} / {n_total:,} qualified alignments ({n_qualified / n_total:.2%}) processed."
        if n_total else
        "No alignments found."
    )

    logging.info("Writing qualified hits...")

    with open(samfile, "r") as fin, open(samfile_temp, "w") as fout:
        fout.writelines(itertools.compress(fin, qualified_mask))

    logging.info(f"Done writing {n_qualified:,} qualified alignments.")

    return n_total, n_qualified

def split_reads_samfile_postprocessing(samfile, samfile_temp):
    """
    Clean up SAM file by removing inconsistent split-read alignments.

    Parameters:
        samfile (str): Path to the input SAM file
        samfile_temp (str): Path to the output cleaned SAM file

    Returns:
        tuple: Total chunks and qualified hits counts
    """
    total_chunks = 0

    logging.info(f'Reading split-read alignments from the sam file...')

    df = pd.read_csv(samfile,
                sep='\t',
                header=None,
                usecols=[0, 2],
                names=['QNAME', 'REF'],
                dtype={'QNAME': 'str', 'REF': 'str'}
    )

    df['QNAME_MAIN'] = df['QNAME'].str.split('|').str[0]
    df['REF_TAXID'] = df['REF'].str.split('|').str[-2]

    logging.info(f'Filtering out inconsistent chunks of reads...')
    # get the index with the most frequent taxid for each read
    idxmax = df.assign(_taxid_cnt=df.groupby(["QNAME_MAIN", "REF_TAXID"])["REF_TAXID"].transform("size")) \
                .loc[lambda x: x["_taxid_cnt"] == x.groupby("QNAME_MAIN")["_taxid_cnt"].transform("max")] \
                .index
    # Create a set of indices for faster lookup
    idxmax_set = set(idxmax.values)

    # in the idxmax, get the index of the first occurrence of the chunk (QNAME_MAIN) for each read
    idx1st = df.loc[idxmax].drop_duplicates(subset="QNAME_MAIN", keep="first").index
    idx1st_set = set(idx1st.values)

    total_chunks = len(df)
    del idxmax
    del df

    logging.info(f'Writing qualified hits...')
    with open(samfile_temp, 'w') as fout, open(samfile, 'r') as fin:
        lines_to_write = []            
        for idx, line in enumerate(fin):
            if idx in idxmax_set:
                if idx in idx1st_set:
                    lines_to_write.append(f"{line.rstrip()}\tZC:i:1\n")
                else:
                    lines_to_write.append(line)

                if len(lines_to_write) >= 1000:
                    fout.writelines(lines_to_write)
                    lines_to_write.clear()
                    logging.debug(f'Written {idx} alignments...')

        if lines_to_write:
            fout.writelines(lines_to_write)

    logging.info(f'Done writing {len(idxmax_set)} hits.')

    return total_chunks, len(idxmax_set)


def main():
    ap = argparse.ArgumentParser(description="Parse and filter SAM file for ONT reads based on species support.")
    ap.add_argument("-s", "--samfile", required=True, help="Input SAM file")
    ap.add_argument("-o", "--out-temp", required=True, help="Output temporary SAM file for qualified hits")
    ap.add_argument("-c", "--custom-tsv", required=True, help="Custom taxonomy TSV file")
    ap.add_argument("--ont-min-species-support", type=float, default=0.6, help="Minimum species support for ONT reads, default 0.6")
    args = ap.parse_args()

    ont_min_species_support = args.ont_min_species_support
    samfile = args.samfile
    samfile_temp = args.out_temp

    logging.basicConfig(
        level=logging.DEBUG,
        format='[%(asctime)s] [%(levelname)s] [%(module)s] %(message)s',
        datefmt='%Y%m%d %H:%M:%S',
    )

    t.loadTaxonomy(cus_taxonomy_file=args.custom_tsv)

    tol_alignment_cnt, tol_q_alignment_cnt = direct_ont_reads_samfile_postprocessing(samfile, samfile_temp, ont_min_species_support)

    print(f'Total alignments: {tol_alignment_cnt}, Total qualified alignments: {tol_q_alignment_cnt}')


if __name__ == "__main__":
    main()
