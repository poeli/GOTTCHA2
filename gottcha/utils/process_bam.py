#!/usr/bin/env python3
"""
process_bam.py

Compute per-region coverage and consensus-mismatch metrics from a BAM file
(with .bai index present), **without** a reference FASTA.

Assumptions / notes:
- Alignments are from minimap2 with `--eqx`, so mismatches are encoded as CIGAR op `X`
  and matches as `=`.
- Depth is computed from aligned query bases (CIGAR ops M/= /X). Deletions (D) and
  refskips (N) do not contribute to depth in this implementation.
- "mismatches" counts total mismatched aligned bases across all reads (sum of `X`).
- "pileup_mismatch" counts reference positions where >50% of aligned reads at that
  position are mismatches (i.e., #positions where X_depth / depth > 0.5).

Parallelization:
- References are split into fixed-size chunks along their length. Each chunk is processed
  independently in a worker process.
- Each worker opens the BAM once (via Pool initializer) for performance.

Output columns (TSV):
- rname
- startpos
- endpos
- numreads
- covbases
- coverage
- mismatches
- pileup_mismatch
- meandepth
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import os
import sys
import logging
from typing import Iterable, List, Optional, Tuple

import numpy as np
import pysam

# Global BAM handle and config for worker processes
_BAM: Optional[pysam.AlignmentFile] = None
_CFG = {}

def _init_worker(
    bam_path: str,
    htslib_threads: int,
    min_mapq: int,
    min_frac: float,
    min_idt: float,
    min_alen: int,
    include_secondary: bool,
    include_supplementary: bool,
    include_duplicates: bool,
    include_qcfail: bool,
    split_read_flag: Optional[bool] = False,
) -> None:
    """Initializer for each worker process: open BAM once and stash filters."""
    global _BAM, _CFG
    _BAM = pysam.AlignmentFile(bam_path, "rb", threads=htslib_threads)
    _CFG = {
        "min_mapq": min_mapq,
        "min_frac": min_frac,
        "min_idt": min_idt,
        "min_alen": min_alen,
        "include_secondary": include_secondary,
        "include_supplementary": include_supplementary,
        "include_duplicates": include_duplicates,
        "include_qcfail": include_qcfail,
        "split_read_flag": split_read_flag,
    }


def _process_chunk(task: Tuple[str, int, int]) -> List:
    """
    Process one (rname, start0, end0) chunk.

    Returns:
      (rname, start0, end0, numreads, covbases, mismatches_total,
       consensus_diff, mean_depth)
    """
    global _BAM, _CFG
    assert _BAM is not None, "Worker BAM handle not initialized"

    rname, start0, end0 = task
    L = end0 - start0
    if L <= 0:
        return [rname, start0, end0, 0,0,0,0,0,0,0,0]

    # Difference arrays (signed) so we can do O(segments) updates and O(L) cumsums.
    # depth[pos] = #reads with an aligned base at that position (from CIGAR ops M/= /X)
    # mm[pos] = #reads with CIGAR X at that position
    depth_diff = np.zeros(L + 1, dtype=np.int32)
    mm_diff = np.zeros(L + 1, dtype=np.int32)

    min_mapq = _CFG["min_mapq"]
    min_frac = _CFG["min_frac"]
    min_idt = _CFG["min_idt"]
    inc_sec = _CFG["include_secondary"]
    inc_sup = _CFG["include_supplementary"]
    inc_dup = _CFG["include_duplicates"]
    inc_qcf = _CFG["include_qcfail"]
    min_alen = _CFG["min_alen"]
    split_read_flag = _CFG["split_read_flag"]

    numreads = 0
    readlength = 0
    indels = 0
    invalid_alns = 0
    bam = _BAM

    # Iterate reads overlapping this region.
    for aln in bam.fetch(rname, start0, end0):
        # Basic filters
        if aln.is_unmapped:
            continue
        if (not inc_sec) and aln.is_secondary:
            continue
        if (not inc_sup) and aln.is_supplementary:
            continue
        if (not inc_dup) and aln.is_duplicate:
            continue
        if (not inc_qcf) and aln.is_qcfail:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        # Note: aln.reference_start is 0-based leftmost coordinate of the alignment on the reference.
        # Only count reads that have their aligned portion starting within the chunk towards numreads, to avoid double-counting reads that span multiple chunks.
        if aln.reference_start >= start0:

            if min_idt > 0.0 and aln.has_tag('NM'):
                mm_idt = (1-aln.get_tag('NM')/aln.alen)
                if min_idt > mm_idt:
                    invalid_alns += 1
                    continue

            if min_frac > 0.0:
                if aln.query_length <= 0:
                    #for hard clips, query_length can be 0, recover it from CIGAR
                    query_length = sum(length for op, length in aln.cigartuples if op in (0, 1, 4, 5, 7, 8)) # M/I/S/H/=/X
                else:
                    query_length = aln.query_length

                if (aln.alen / query_length) < min_frac and (aln.alen / bam.get_reference_length(rname)) < min_frac:
                    invalid_alns += 1
                    continue

            if min_alen > 0 and aln.alen < min_alen:
                invalid_alns += 1
                continue

            # If split_read_flag is set, only count reads with ZC tag (the first chunked reads) towards numreads.
            if split_read_flag:
                if aln.has_tag('ZC'):
                    numreads += 1
            else:
                if aln.is_secondary or aln.is_supplementary:
                    pass
                else:
                    numreads += 1
            
            # count total read length (including softclips) for mean depth calculation
            readlength += aln.alen

        cig = aln.cigartuples
        if not cig:
            continue

        ref_pos = aln.reference_start
        block_start: Optional[int] = None

        # CIGAR operation codes in pysam:
        # 0=M, 1=I, 2=D, 3=N, 4=S, 5=H, 6=P, 7==, 8=X
        for op, length in cig:
            if length <= 0:
                continue

            if op in (0, 7, 8):  # aligned query bases consuming reference
                if block_start is None:
                    block_start = ref_pos

                if op == 8:  # X mismatches
                    seg_s = ref_pos
                    seg_e = ref_pos + length  # exclusive
                    if seg_e > start0 and seg_s < end0:
                        if seg_s < start0:
                            seg_s = start0
                        if seg_e > end0:
                            seg_e = end0
                        mm_diff[seg_s - start0] += 1
                        mm_diff[seg_e - start0] -= 1

                ref_pos += length

            elif op in (2, 3):  # D or N: consumes reference but not query -> breaks aligned block
                if block_start is not None:
                    seg_s = block_start
                    seg_e = ref_pos  # exclusive
                    if seg_e > start0 and seg_s < end0:
                        if seg_s < start0:
                            seg_s = start0
                        if seg_e > end0:
                            seg_e = end0
                        depth_diff[seg_s - start0] += 1
                        depth_diff[seg_e - start0] -= 1
                        indels += length
                    block_start = None

                ref_pos += length

            else:
                # I/S/H/P: does not consume reference; does not affect ref_pos.
                # We do not break the block on insertions/softclips because reference
                # positions remain contiguous.
                continue

        # Close any remaining aligned block
        if block_start is not None:
            seg_s = block_start
            seg_e = ref_pos
            if seg_e > start0 and seg_s < end0:
                if seg_s < start0:
                    seg_s = start0
                if seg_e > end0:
                    seg_e = end0
                depth_diff[seg_s - start0] += 1
                depth_diff[seg_e - start0] -= 1

    # Build per-base depth and mismatch arrays
    depth = np.cumsum(depth_diff[:-1], dtype=np.int32)  # length L
    mm = np.cumsum(mm_diff[:-1], dtype=np.int32)        # length L

    covbases = int(np.count_nonzero(depth))
    mismatches_total = int(mm.sum())
    mapped_bases = int(depth.sum())+indels # total aligned bases (including matches and mismatches)

    # Positions where mismatch fraction > 0.5 among reads with aligned bases
    # i.e. mm / depth > 0.5  ->  2*mm > depth
    consensus_diff = int(np.count_nonzero((depth > 0) & (mm * 2 > depth)))

    logging.debug(f"Processed {rname}: {numreads} reads, {covbases} covbases, {mismatches_total} mismatches, {consensus_diff} consensus_diff, {indels} indels, {mapped_bases} mapped_bases, {invalid_alns} invalid_alns")

    return [rname,
            start0,
            end0,
            numreads,
            covbases,
            mismatches_total,
            indels,
            consensus_diff,
            mapped_bases,
            invalid_alns,
            readlength]


def _iter_tasks(references: List[str], lengths: List[int], chunk_size: int) -> Iterable[Tuple[str, int, int]]:
    for rname, rlen in zip(references, lengths):
        if rlen <= 0:
            continue
        cs = rlen if chunk_size <= 0 else chunk_size
        for start0 in range(0, rlen, cs):
            end0 = min(start0 + cs, rlen)
            yield (rname, start0, end0)

def parse_aln_from_bam(bam_path: str,
                       processes: int, 
                       min_frac: float,
                       min_idt: float,
                       min_alen: int,
                       min_mapq: Optional[int] = 0,
                       htslib_threads: Optional[int] = 1,
                       chunk_size: Optional[int] = 10_000_000,
                       imap_chunksize: Optional[int] = 1,
                       include_secondary: Optional[bool] = False,
                       include_supplementary: Optional[bool] = False,
                       include_duplicates: Optional[bool] = False,
                       include_qcfail: Optional[bool] = False,
                       split_read_flag: Optional[bool] = False,
                       ) -> int:
    if not os.path.exists(bam_path):
        print(f"ERROR: BAM not found: {bam_path}", file=sys.stderr)
        return 2

    # Open BAM in main process to validate index and obtain reference lengths
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if not bam.has_index():
                print("ERROR: BAM index (.bai) not found or not readable. Pysam requires an index.", file=sys.stderr)
                return 2
            references = list(bam.references)
            lengths = list(bam.lengths)
    except Exception as e:
        print(f"ERROR: Failed to open BAM: {e}", file=sys.stderr)
        return 2

    logging.debug(f"Parsing {len(references)} references with {processes} processes...")
    logging.debug(f"Parameters: min_mapq={min_mapq}, min_frac={min_frac}, min_idt={min_idt}, min_alen={min_alen}, include_secondary={include_secondary}, include_supplementary={include_supplementary}, include_duplicates={include_duplicates}, include_qcfail={include_qcfail}, split_read_flag={split_read_flag}")

    tasks = _iter_tasks(references, lengths, chunk_size)

    # (default is one-based if neither flag used; argparse sets one_based True by default)
    # endpos will be end0 in both conventions; interpretation differs.

    pool = mp.Pool(
        processes=processes,
        initializer=_init_worker,
        initargs=(
            bam_path,
            htslib_threads,
            min_mapq,
            min_frac,
            min_idt,
            min_alen,
            include_secondary,
            include_supplementary,
            include_duplicates,
            include_qcfail,
            split_read_flag,
        ),
    )

    try:
        ref_chunk_results = []
        header = [
            "RNAME",             # rname,
            "STARTPOS",          # start0,
            "ENDPOS",            # end0,
            "NUMREADS",          # numreads,
            "COVBASES",          # covbases,
            "MISMATCHES",        # mismatches_total,
            "INDELS",            # indels,
            "CONSENSUS_DIFF",    # consensus_diff,
            "MAPPED_BASES",      # mapped_bases,
            "INVALID_ALNS",      # invalid_alns,
            "READLENGTH"         # readlength
        ]

        ref_chunk_results.append(header)
        mapper = pool.imap_unordered
        for result in mapper(_process_chunk, tasks, chunksize=imap_chunksize):
            result[1] +=  1 # start0 to 1-based startpos for output (endpos remains end0, which is exclusive in both conventions)
            ref_chunk_results.append(result)

    finally:
        pool.close()
        pool.join()

    logging.debug(f"Total signature fragments processed: {len(ref_chunk_results)-1}")

    return ref_chunk_results


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Compute coverage and consensus mismatch metrics from a BAM"
    )
    p.add_argument("bam", help="Input BAM path (requires .bai index).")
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    p.add_argument(
        "-c",
        "--chunk-size",
        type=int,
        default=1_000_000,
        help="Chunk size in reference bases for parallel tasks (default: 1,000,000). Use 0 for whole-contig.",
    )
    p.add_argument(
        "-p",
        "--processes",
        type=int,
        default=max(1, mp.cpu_count() - 1),
        help="Worker processes (default: cpu_count-1).",
    )
    p.add_argument(
        "-t",
        "--htslib-threads",
        type=int,
        default=1,
        help="HTSlib threads per worker for BAM decompression (default: 1).",
    )

    # Filters
    p.add_argument("--min-mapq", type=int, default=0, help="Minimum MAPQ to keep an alignment (default: 0).")
    p.add_argument("--min-frac", type=float, default=0.0, help="Minimum fraction to keep an alignment (default: 0.0).")
    p.add_argument("--min-idt", type=float, default=0.0, help="Minimum identity to keep an alignment (default: 0.0).")
    p.add_argument("--min-alen", type=int, default=0, help="Minimum alignment length to keep an alignment (default: 0).")
    p.add_argument("--include-secondary", action="store_true", help="Include secondary alignments (default: off).")
    p.add_argument("--include-supplementary", action="store_true", help="Include supplementary alignments (default: off).")
    p.add_argument("--include-duplicates", action="store_true", help="Include duplicate-marked reads (default: off).")
    p.add_argument("--include-qcfail", action="store_true", help="Include QC-failed reads (default: off).")

    # Coordinate output style
    p.add_argument(
        "--imap-chunksize",
        type=int,
        default=1,
        help="chunksize passed to multiprocessing imap/imap_unordered (default: 1).",
    )

    args = p.parse_args(argv)

    ref_results = parse_aln_from_bam(
        bam_path=args.bam,
        processes=args.processes,
        min_frac=args.min_frac,
        min_idt=args.min_idt,
        min_alen=args.min_alen,
        min_mapq=args.min_mapq,
        htslib_threads=args.htslib_threads,
        chunk_size=args.chunk_size,
        imap_chunksize=args.imap_chunksize,
        include_secondary=args.include_secondary,
        include_supplementary=args.include_supplementary,
        include_duplicates=args.include_duplicates,
        include_qcfail=args.include_qcfail,
    )
    
    out_path=args.out
    with open(out_path, "w", encoding="utf-8") as out:
        for res in ref_results:
            out.write("\t".join(map(str, res)) + "\n")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())
