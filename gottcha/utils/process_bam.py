#!/usr/bin/env python3
"""
process_bam.py

Compute per-region coverage and consensus-mismatch metrics from a BAM file
(with .bai index present), without a reference FASTA.

Also identify connected groups of species linked by qualified read mappings.

Example
-------
read1 -> sp1
read2 -> sp1, sp2, sp3
read3 -> sp2, sp4
read4 -> sp5
read5 -> sp6, sp7

Produces:

GROUP 1: sp1, sp2, sp3, sp4
GROUP 2: sp6, sp7

sp5 is omitted because it is not linked to another species.

Species are resolved from reference names using:

    taxid = rname.split('|')[-2]
    spe_taxid = t.taxid2taxidOnRank(taxid, target_rank="species")


Coverage assumptions / notes
----------------------------
- Alignments are from minimap2 with --eqx, so mismatches are encoded as
  CIGAR op X and matches as =.
- Depth is computed from aligned query bases (M/= /X).
- Deletions and refskips do not contribute to depth.
- "mismatches" counts total mismatched aligned bases across all reads.
- "pileup_mismatch" / CONSENSUS_DIFF counts reference positions where
  >50% of aligned reads at that position are mismatches.

Species-linkage algorithm
-------------------------
A separate fast streaming BAM pass uses union-find.

For each read, only the first observed species needs to be remembered:

    read -> anchor species

Subsequent qualified alignments from the same read are unioned with that
anchor. Thus:

    read -> sp1, sp2, sp3

requires only:

    union(sp1, sp2)
    union(sp1, sp3)

rather than constructing every pairwise species edge.

The resulting connected components represent groups of species linked
transitively by reads.

Secondary alignments
--------------------
Secondary alignments are INCLUDED for species linkage by default because
minimap2 normally represents alternative mappings as secondary alignments.

This is independent of --include-secondary, which still controls whether
secondary alignments contribute to coverage/mismatch calculations.

Use --groups-primary-only to prevent secondary alignments from being used
for species linkage.

Parallelization
---------------
- References are split into fixed-size chunks.
- Each chunk is processed independently by a worker.
- Each worker opens the BAM once.

Global group lookups
--------------------
_GROUPS:
    species_taxid -> group_id

_GROUP_MEMBERS:
    group_id -> list of species_taxids
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import os
import sys
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pysam

try:
    from . import taxonomy as t
except ImportError:
    import taxonomy as t


# ----------------------------------------------------------------------
# Globals
# ----------------------------------------------------------------------

# BAM handle/config for worker processes
_BAM: Optional[pysam.AlignmentFile] = None
_CFG = {}

# Species linkage results in the main process:
#
#   _GROUPS["12345"] = 1
#
# means species taxid 12345 belongs to linked group 1.
_GROUPS: Dict[str, int] = {}

# Reverse lookup:
#
#   _GROUP_MEMBERS[1] = ["12345", "67890", ...]
_GROUP_MEMBERS: Dict[int, List[str]] = {}


# ----------------------------------------------------------------------
# Alignment filtering
# ----------------------------------------------------------------------

def _alignment_filter_reason(
    aln: pysam.AlignedSegment,
    bam: pysam.AlignmentFile,
    cfg: dict,
) -> Optional[str]:
    """
    Return None if the alignment passes all filters.

    Otherwise return a rejection reason.

    The filtering logic is shared by:
      1. coverage/mismatch calculation
      2. species-linkage discovery
    """

    if aln.is_unmapped:
        return "unmapped"

    if (not cfg["include_secondary"]) and aln.is_secondary:
        return "secondary"

    if (not cfg["include_supplementary"]) and aln.is_supplementary:
        return "supplementary"

    if (not cfg["include_duplicates"]) and aln.is_duplicate:
        return "duplicate"

    if (not cfg["include_qcfail"]) and aln.is_qcfail:
        return "qcfail"

    if aln.mapping_quality < cfg["min_mapq"]:
        return "mapq"

    alen = aln.alen or 0

    if alen <= 0:
        return "length"

    # --------------------------------------------------------------
    # Alignment identity
    # --------------------------------------------------------------
    min_idt = cfg["min_idt"]

    if min_idt > 0.0 and aln.has_tag("NM"):
        nm = aln.get_tag("NM")
        identity = 1.0 - (nm / alen)

        if identity < min_idt:
            return "identity"

    # --------------------------------------------------------------
    # Fraction of query/reference represented by alignment
    # --------------------------------------------------------------
    min_frac = cfg["min_frac"]

    if min_frac > 0.0:

        query_length = aln.query_length

        if not query_length:
            # Recover length for hard-clipped reads from CIGAR.
            #
            # CIGAR:
            # 0=M, 1=I, 4=S, 5=H, 7==, 8=X
            cig = aln.cigartuples or ()

            query_length = sum(
                length
                for op, length in cig
                if op in (0, 1, 4, 5, 7, 8)
            )

        if query_length <= 0:
            return "fraction"

        reference_length = bam.lengths[aln.reference_id]

        if (
            alen / query_length < min_frac
            and
            alen / reference_length < min_frac
        ):
            return "fraction"

    # --------------------------------------------------------------
    # Minimum aligned length
    # --------------------------------------------------------------
    min_alen = cfg["min_alen"]

    if min_alen > 0 and alen < min_alen:
        return "length"

    return None


# ----------------------------------------------------------------------
# Taxonomy helpers
# ----------------------------------------------------------------------

def _rname_to_species_taxid(rname: str) -> Optional[str]:
    """
    Convert a BAM reference name to a species-level taxid.

    Expected reference-name convention:

        ...|taxid|...

    where:

        taxid = rname.split('|')[-2]
    """

    try:
        parts = rname.split("|")

        if len(parts) < 2:
            return None

        taxid = parts[-2]

        spe_taxid = t.taxid2taxidOnRank(taxid, target_rank="species",)

        if spe_taxid is None:
            return None

        if not spe_taxid or spe_taxid == "0":
            return None

        return spe_taxid

    except Exception as exc:
        logging.debug(
            "Unable to resolve species taxid from reference %s: %s",
            rname,
            exc,
        )
        return None


def _taxid_sort_key(taxid: str):
    """Sort numeric taxids numerically when possible."""

    try:
        return 0, int(taxid)
    except (TypeError, ValueError):
        return 1, str(taxid)


# ----------------------------------------------------------------------
# Species linkage
# ----------------------------------------------------------------------

def find_linked_species_groups(
    bam_path: str,
    min_frac: float,
    min_idt: float,
    min_alen: int,
    min_mapq: int = 0,
    htslib_threads: int = 1,
    include_secondary: bool = True,
    include_supplementary: bool = False,
    include_duplicates: bool = False,
    include_qcfail: bool = False,
) -> Tuple[Dict[int, List[str]], Dict[str, int]]:
    """
    Find connected groups of species linked by qualified read alignments.

    Returns
    -------
    groups
        group_id -> list of species taxids

        Example:
            {
                1: ["100", "101", "102"],
                2: ["200", "201"],
            }

    species_to_group
        species taxid -> group_id

        Example:
            {
                "100": 1,
                "101": 1,
                "102": 1,
                "200": 2,
                "201": 2,
            }

    Notes
    -----
    Singleton species are deliberately omitted.

    Complexity is approximately O(number of qualified alignments), because
    union-find operations are effectively constant-time.
    """

    cfg = {
        "min_mapq": min_mapq,
        "min_frac": min_frac,
        "min_idt": min_idt,
        "min_alen": min_alen,
        "include_secondary": include_secondary,
        "include_supplementary": include_supplementary,
        "include_duplicates": include_duplicates,
        "include_qcfail": include_qcfail,
    }

    # --------------------------------------------------------------
    # Map species taxids to compact integer indices.
    #
    # This makes the union-find arrays and per-read anchors much smaller
    # and faster than repeatedly storing taxid strings.
    # --------------------------------------------------------------

    species_to_idx: Dict[str, int] = {}
    idx_to_species: List[str] = []

    parent: List[int] = []
    component_size: List[int] = []

    def get_species_idx(species_taxid: str) -> int:
        idx = species_to_idx.get(species_taxid)

        if idx is not None:
            return idx

        idx = len(idx_to_species)

        species_to_idx[species_taxid] = idx
        idx_to_species.append(species_taxid)

        parent.append(idx)
        component_size.append(1)

        return idx

    def find(x: int) -> int:
        """
        Iterative union-find root lookup with path compression.
        """

        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]

        return x

    def union(a: int, b: int) -> None:
        """
        Union by component size.
        """

        root_a = find(a)
        root_b = find(b)

        if root_a == root_b:
            return

        if component_size[root_a] < component_size[root_b]:
            root_a, root_b = root_b, root_a

        parent[root_b] = root_a
        component_size[root_a] += component_size[root_b]

    # --------------------------------------------------------------
    # Remember only one anchor species for each read.
    #
    # Example:
    #
    #   read1 -> sp1, sp2, sp3
    #
    # produces:
    #
    #   union(sp1, sp2)
    #   union(sp1, sp3)
    #
    # rather than storing {sp1, sp2, sp3}.
    # --------------------------------------------------------------

    read_anchor: Dict[str, int] = {}

    n_alignments = 0
    n_qualified = 0
    n_species_resolved = 0
    n_cross_species_links = 0

    with pysam.AlignmentFile(
        bam_path,
        "rb",
        threads=htslib_threads,
    ) as bam:

        # Cache:
        #
        #   reference_id -> species index
        #
        # Sentinel values:
        #   -2 = not evaluated
        #   -1 = could not resolve species
        #
        # Taxonomy lookup can be relatively expensive, so we perform it
        # only once per BAM reference.
        ref_species_idx = [-2] * bam.nreferences

        for aln in bam.fetch(until_eof=True):

            n_alignments += 1

            reason = _alignment_filter_reason(
                aln,
                bam,
                cfg,
            )

            if reason is not None:
                continue

            n_qualified += 1

            ref_id = aln.reference_id

            species_idx = ref_species_idx[ref_id]

            # ------------------------------------------------------
            # Resolve reference -> species only once.
            # ------------------------------------------------------
            if species_idx == -2:

                rname = bam.references[ref_id]

                species_taxid = _rname_to_species_taxid(rname)

                if species_taxid is None:
                    ref_species_idx[ref_id] = -1
                    continue

                species_idx = get_species_idx(species_taxid)

                ref_species_idx[ref_id] = species_idx

                n_species_resolved += 1

            elif species_idx == -1:
                continue

            # ------------------------------------------------------
            # Link species through query/read name.
            # ------------------------------------------------------
            qname = aln.query_name

            anchor_idx = read_anchor.get(qname)

            if anchor_idx is None:
                read_anchor[qname] = species_idx
                continue

            # Same read mapping to another reference within the same
            # species does not create a species-to-species link.
            if anchor_idx == species_idx:
                continue

            union(anchor_idx, species_idx)
            n_cross_species_links += 1

    n_reads = len(read_anchor)

    # No longer needed after the union-find graph is complete.
    del read_anchor

    # --------------------------------------------------------------
    # Convert connected components back to species taxids.
    # --------------------------------------------------------------

    components: Dict[int, List[str]] = defaultdict(list)

    for idx, species_taxid in enumerate(idx_to_species):

        root = find(idx)

        components[root].append(species_taxid)

    linked_components: List[List[str]] = []

    for members in components.values():

        # Exclude singleton species.
        if len(members) < 2:
            continue

        members.sort(key=_taxid_sort_key)

        linked_components.append(members)

    # Deterministic group ordering.
    linked_components.sort(
        key=lambda members: _taxid_sort_key(members[0])
    )

    groups: Dict[int, List[str]] = {}
    species_to_group: Dict[str, int] = {}

    for group_id, members in enumerate(
        linked_components,
        start=1,
    ):

        groups[group_id] = members

        for species_taxid in members:
            species_to_group[species_taxid] = group_id

    logging.info(
        "Species linkage: "
        "%d BAM alignments, "
        "%d qualified alignments, "
        "%d reads, "
        "%d resolved species, "
        "%d cross-species links, "
        "%d connected species groups",
        n_alignments,
        n_qualified,
        n_reads,
        len(idx_to_species),
        n_cross_species_links,
        len(groups),
    )

    return groups, species_to_group


def write_species_groups(
    groups: Dict[int, List[str]],
    out_path: str,
) -> None:
    """
    Write connected species groups as TSV.

    Format:

        GROUP_ID    N_SPECIES    SPECIES_TAXIDS
        1           4            11,22,33,44
        2           2            55,66
    """

    with open(out_path, "w", encoding="utf-8") as out:

        out.write(
            "GROUP_ID\t"
            "N_SPECIES\t"
            "SPECIES_TAXIDS\n"
        )

        for group_id, species in groups.items():

            out.write(
                f"{group_id}\t"
                f"{len(species)}\t"
                f"{','.join(species)}\n"
            )


# ----------------------------------------------------------------------
# Multiprocessing worker
# ----------------------------------------------------------------------

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
    """
    Initializer for each worker process.

    Open BAM once and stash filters.
    """

    global _BAM, _CFG

    _BAM = pysam.AlignmentFile(
        bam_path,
        "rb",
        threads=htslib_threads,
    )

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
        # IMPORTANT:
        #
        # Apply filters to ALL alignments contributing to this chunk,
        # including alignments whose reference_start lies in a previous
        # chunk.
        #
        # In the previous implementation, threshold filters were inside:
        #
        #     if aln.reference_start >= start0
        #
        # so a read spanning a chunk boundary could contribute depth in
        # the next chunk even if it failed min-idt/min-frac/min-alen.
        filter_reason = _alignment_filter_reason(aln, bam, _CFG,)

        if filter_reason is not None:
            # Preserve existing INVALID_ALNS semantics:
            # only identity/fraction/length failures are counted, and
            # only in the chunk where the alignment starts.
            if (
                aln.reference_start >= start0
                and
                filter_reason
                in {
                    "identity",
                    "fraction",
                    "length",
                }
            ):
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


    # find_linked_species_groups
    logging.info("Finding linked species groups...")
    groups, species_to_group = find_linked_species_groups(bam_path, 
                               min_frac=min_frac,
                               min_idt=min_idt,
                               min_alen=min_alen,
                               min_mapq=min_mapq,
                               htslib_threads=htslib_threads,
                               include_secondary=include_secondary,
                               include_supplementary=include_supplementary,
                               include_duplicates=include_duplicates,
                               include_qcfail=include_qcfail)
    logging.info(f"Found {len(groups)} linked species groups.")
    logging.debug(f"Species to group mapping: {species_to_group}")
    logging.debug(f"Group members: {groups}")

    write_species_groups(groups, f"{bam_path}.species_groups.tsv")

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
