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
Species-linkage observations are collected inside the same parallel chunk
pass used for coverage, so the BAM is not scanned a second time. Each worker
returns the unique qualified read names whose alignments start in its chunk,
together with that chunk's species index. The parent merges these observations
into a unique set of species for each read.

After the BAM pass, linkage is treated as a weighted species graph. For every
read observed on multiple species, each species pair receives one unit of
shared-read support. An edge is retained only when both:

    shared_reads >= min_link_reads

and
    shared_reads / min(reads_species1, reads_species2) >= min_link_fraction

Optionally, --same-genus-only restricts retained edges to species assigned to
the same genus. --max-link-species-per-read can exclude highly ambiguous reads
from pair support while still counting them in per-species read totals.

Only retained edges are unioned. Connected components of retained edges form
the final linked-species groups. Diagnostic TSVs report every observed edge
and the distribution of species-per-read values.

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
from collections import Counter, defaultdict
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

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


def _rname_to_species_taxid(rname: str) -> Optional[str]:
    """Convert a BAM reference name to a species-level taxid."""

    try:
        taxid = rname.split("|")[-2]
        spe_taxid = t.taxid2taxidOnRank(taxid, target_rank="species")
        if spe_taxid is None:
            return None
        spe_taxid = str(spe_taxid)
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


def _taxid_to_rank_taxid(taxid: str, target_rank: str) -> Optional[str]:
    """Resolve a taxid to a requested rank, returning a normalized string."""

    try:
        rank_taxid = t.taxid2taxidOnRank(taxid, target_rank=target_rank)
        if rank_taxid is None:
            return None
        rank_taxid = str(rank_taxid)
        if not rank_taxid or rank_taxid == "0":
            return None
        return rank_taxid
    except Exception as exc:
        logging.debug(
            "Unable to resolve taxid %s to rank %s: %s",
            taxid,
            target_rank,
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
# Species linkage helpers
# ----------------------------------------------------------------------

ReadSpeciesState = Union[int, Set[int]]


def _add_read_species_observation(
    read_species: Dict[str, ReadSpeciesState],
    qname: str,
    species_idx: int,
) -> bool:
    """
    Add one read->species observation, deduplicating globally across chunks.

    To reduce memory for the common case, a read mapping to exactly one species
    is stored as a single integer. It is promoted to a set only after a second
    distinct species is observed.

    Returns True only when this is a new species observation for the read.
    """

    state = read_species.get(qname)

    if state is None:
        read_species[qname] = species_idx
        return True

    if isinstance(state, int):
        if state == species_idx:
            return False
        read_species[qname] = {state, species_idx}
        return True

    if species_idx in state:
        return False

    state.add(species_idx)
    return True


def _build_species_link_support(
    read_species: Dict[str, ReadSpeciesState],
    n_species: int,
    max_link_species_per_read: int = 0,
) -> Tuple[List[int], Dict[Tuple[int, int], int], Counter, int, int]:
    """
    Convert per-read species observations into weighted pair support.

    species_read_counts counts every distinct read observed on each species.
    Pair support is counted once per read per species pair. If
    max_link_species_per_read > 0, reads mapping to more than that many species
    are excluded from pair support but remain in species_read_counts and the
    degree histogram.

    Returns
    -------
    species_read_counts
        Unique qualified reads observed for each species.
    edge_read_counts
        (species_idx1, species_idx2) -> number of shared reads.
    read_degree_hist
        number_of_species_per_read -> number_of_reads.
    n_multi_species_reads
        Number of reads observed on >=2 species.
    n_excluded_hyperambiguous_reads
        Multi-species reads excluded from edge support by the maximum-degree
        filter.
    """

    species_read_counts = [0] * n_species
    edge_read_counts: Dict[Tuple[int, int], int] = defaultdict(int)
    read_degree_hist: Counter = Counter()
    n_multi_species_reads = 0
    n_excluded_hyperambiguous_reads = 0

    for state in read_species.values():
        if isinstance(state, int):
            species = (state,)
        else:
            species = tuple(sorted(state))

        degree = len(species)
        read_degree_hist[degree] += 1

        for species_idx in species:
            species_read_counts[species_idx] += 1

        if degree < 2:
            continue

        n_multi_species_reads += 1

        if max_link_species_per_read > 0 and degree > max_link_species_per_read:
            n_excluded_hyperambiguous_reads += 1
            continue

        for i in range(degree - 1):
            a = species[i]
            for j in range(i + 1, degree):
                b = species[j]
                edge_read_counts[(a, b)] += 1

    return (
        species_read_counts,
        edge_read_counts,
        read_degree_hist,
        n_multi_species_reads,
        n_excluded_hyperambiguous_reads,
    )


def _finalize_species_groups(
    idx_to_species: List[str],
    idx_to_genus: List[Optional[str]],
    species_read_counts: List[int],
    edge_read_counts: Dict[Tuple[int, int], int],
    min_link_reads: int,
    min_link_fraction: float,
    same_genus_only: bool,
) -> Tuple[
    Dict[int, List[str]],
    Dict[str, int],
    List[Tuple],
    int,
]:
    """
    Filter weighted species edges and construct deterministic components.

    link_fraction is a containment-style measure:

        shared_reads / min(reads_species1, reads_species2)

    A species pair is unioned only when it passes min_link_reads,
    min_link_fraction, and (optionally) the same-genus requirement.
    """

    n_species = len(idx_to_species)
    parent = list(range(n_species))
    component_size = [1] * n_species

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: int, b: int) -> bool:
        root_a = find(a)
        root_b = find(b)

        if root_a == root_b:
            return False

        if component_size[root_a] < component_size[root_b]:
            root_a, root_b = root_b, root_a

        parent[root_b] = root_a
        component_size[root_a] += component_size[root_b]
        return True

    # Store edge diagnostics before assigning final group IDs.
    edge_rows_raw = []
    n_component_unions = 0

    for (idx1, idx2), shared_reads in edge_read_counts.items():
        reads1 = species_read_counts[idx1]
        reads2 = species_read_counts[idx2]
        denominator = min(reads1, reads2)
        link_fraction = shared_reads / denominator if denominator > 0 else 0.0

        genus1 = idx_to_genus[idx1]
        genus2 = idx_to_genus[idx2]
        same_genus = (
            genus1 is not None
            and genus2 is not None
            and genus1 == genus2
        )

        pass_min_reads = shared_reads >= min_link_reads
        pass_min_fraction = link_fraction >= min_link_fraction
        pass_genus = (not same_genus_only) or same_genus
        retained = pass_min_reads and pass_min_fraction and pass_genus

        if retained and union(idx1, idx2):
            n_component_unions += 1

        edge_rows_raw.append(
            (
                idx1,
                idx2,
                shared_reads,
                reads1,
                reads2,
                link_fraction,
                genus1,
                genus2,
                same_genus,
                pass_min_reads,
                pass_min_fraction,
                pass_genus,
                retained,
            )
        )

    components: Dict[int, List[str]] = defaultdict(list)

    for idx, species_taxid in enumerate(idx_to_species):
        components[find(idx)].append(species_taxid)

    linked_components: List[List[str]] = []

    for members in components.values():
        if len(members) < 2:
            continue
        members.sort(key=_taxid_sort_key)
        linked_components.append(members)

    linked_components.sort(
        key=lambda members: _taxid_sort_key(members[0])
    )

    groups: Dict[int, List[str]] = {}
    species_to_group: Dict[str, int] = {}

    for group_id, members in enumerate(linked_components, start=1):
        groups[group_id] = members
        for species_taxid in members:
            species_to_group[species_taxid] = group_id

    edge_rows: List[Tuple] = []

    for row in edge_rows_raw:
        (
            idx1,
            idx2,
            shared_reads,
            reads1,
            reads2,
            link_fraction,
            genus1,
            genus2,
            same_genus,
            pass_min_reads,
            pass_min_fraction,
            pass_genus,
            retained,
        ) = row

        species1 = idx_to_species[idx1]
        species2 = idx_to_species[idx2]
        group1 = species_to_group.get(species1)
        group2 = species_to_group.get(species2)
        same_final_group = (
            group1 is not None
            and group2 is not None
            and group1 == group2
        )

        edge_rows.append(
            (
                species1,
                species2,
                shared_reads,
                reads1,
                reads2,
                link_fraction,
                genus1,
                genus2,
                same_genus,
                pass_min_reads,
                pass_min_fraction,
                pass_genus,
                retained,
                group1,
                group2,
                same_final_group,
            )
        )

    edge_rows.sort(
        key=lambda row: (
            not row[12],          # retained edges first
            -row[2],              # then strongest shared-read support
            -row[5],              # then strongest link fraction
            _taxid_sort_key(row[0]),
            _taxid_sort_key(row[1]),
        )
    )

    return groups, species_to_group, edge_rows, n_component_unions


def write_species_groups(
    groups: Dict[int, List[str]],
    out_path: str,
) -> None:
    """Write connected species groups as TSV."""

    with open(out_path, "w", encoding="utf-8") as out:
        out.write("GROUP_ID\tN_SPECIES\tSPECIES_TAXIDS\n")

        for group_id, species in groups.items():
            out.write(
                f"{group_id}\t"
                f"{len(species)}\t"
                f"{','.join(species)}\n"
            )


def write_species_links(
    edge_rows: List[Tuple],
    out_path: str,
) -> None:
    """Write weighted species-link diagnostics as TSV."""

    header = [
        "SPECIES1",
        "SPECIES2",
        "SHARED_READS",
        "READS_SPECIES1",
        "READS_SPECIES2",
        "LINK_FRACTION",
        "GENUS1",
        "GENUS2",
        "SAME_GENUS",
        "PASS_MIN_READS",
        "PASS_MIN_FRACTION",
        "PASS_GENUS",
        "RETAINED",
        "GROUP1",
        "GROUP2",
        "SAME_FINAL_GROUP",
    ]

    with open(out_path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")

        for row in edge_rows:
            values = list(row)
            values[5] = f"{values[5]:.8f}"
            values[6] = "" if values[6] is None else values[6]
            values[7] = "" if values[7] is None else values[7]
            values[13] = "" if values[13] is None else values[13]
            values[14] = "" if values[14] is None else values[14]
            out.write("\t".join(map(str, values)) + "\n")


def write_linkage_read_degree_histogram(
    read_degree_hist: Counter,
    out_path: str,
) -> None:
    """Write the number of distinct mapped species observed per read."""

    with open(out_path, "w", encoding="utf-8") as out:
        out.write("N_SPECIES_PER_READ\tN_READS\n")
        for degree in sorted(read_degree_hist):
            out.write(f"{degree}\t{read_degree_hist[degree]}\n")


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


def _process_chunk(task: Tuple[str, int, int, Optional[int]]) -> Tuple[List, Optional[int], List[str]]:
    """
    Process one (rname, start0, end0, species_idx) chunk.

    In addition to the coverage/mismatch metrics, return the qualified read
    names whose alignments START in this chunk. The parent process combines
    these observations across chunks/references to build the linked-species
    connected components.

    Returning only alignments that start in this chunk prevents an alignment
    spanning a chunk boundary from being emitted more than once.
    """
    global _BAM, _CFG
    assert _BAM is not None, "Worker BAM handle not initialized"

    rname, start0, end0, species_idx = task
    L = end0 - start0

    if L <= 0:
        metrics = [rname, start0, end0, 0, 0, 0, 0, 0, 0, 0, 0]
        return metrics, species_idx, []

    # Difference arrays (signed) so we can do O(segments) updates and O(L)
    # cumulative sums.
    depth_diff = np.zeros(L + 1, dtype=np.int32)
    mm_diff = np.zeros(L + 1, dtype=np.int32)

    split_read_flag = _CFG["split_read_flag"]

    numreads = 0
    readlength = 0
    indels = 0
    invalid_alns = 0
    bam = _BAM

    # A chunk corresponds to one reference, hence one species_idx. We only
    # need unique read names for this species within this chunk.
    linkage_reads = set()

    # Iterate alignments overlapping this region.
    for aln in bam.fetch(rname, start0, end0):
        # Apply filters to ALL alignments contributing coverage to this chunk,
        # including alignments that start in the preceding chunk.
        filter_reason = _alignment_filter_reason(aln, bam, _CFG)

        if filter_reason is not None:
            # Count threshold failures once, in the chunk where the alignment
            # starts.
            if (
                aln.reference_start >= start0
                and filter_reason in {"identity", "fraction", "length"}
            ):
                invalid_alns += 1
            continue

        # An alignment overlapping multiple chunks contributes depth to each
        # overlapping chunk, but NUMREADS/READLENGTH and linkage must be
        # recorded exactly once.
        owns_alignment = aln.reference_start >= start0

        if owns_alignment:
            if species_idx is not None and aln.query_name:
                linkage_reads.add(aln.query_name)

            if split_read_flag:
                if aln.has_tag("ZC"):
                    numreads += 1
            else:
                if not aln.is_secondary and not aln.is_supplementary:
                    numreads += 1

            # Preserve the original use of aligned length for READLENGTH.
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
                    seg_e = ref_pos + length
                    if seg_e > start0 and seg_s < end0:
                        seg_s = max(seg_s, start0)
                        seg_e = min(seg_e, end0)
                        mm_diff[seg_s - start0] += 1
                        mm_diff[seg_e - start0] -= 1

                ref_pos += length

            elif op in (2, 3):  # D or N: consumes reference, not query
                if block_start is not None:
                    seg_s = block_start
                    seg_e = ref_pos
                    if seg_e > start0 and seg_s < end0:
                        seg_s = max(seg_s, start0)
                        seg_e = min(seg_e, end0)
                        depth_diff[seg_s - start0] += 1
                        depth_diff[seg_e - start0] -= 1
                        indels += length
                    block_start = None

                ref_pos += length

            else:
                # I/S/H/P do not consume reference and do not break the
                # contiguous reference block.
                continue

        # Close any remaining aligned block.
        if block_start is not None:
            seg_s = block_start
            seg_e = ref_pos
            if seg_e > start0 and seg_s < end0:
                seg_s = max(seg_s, start0)
                seg_e = min(seg_e, end0)
                depth_diff[seg_s - start0] += 1
                depth_diff[seg_e - start0] -= 1

    # Build per-base depth and mismatch arrays.
    depth = np.cumsum(depth_diff[:-1], dtype=np.int32)
    mm = np.cumsum(mm_diff[:-1], dtype=np.int32)

    covbases = int(np.count_nonzero(depth))
    mismatches_total = int(mm.sum())
    mapped_bases = int(depth.sum()) + indels

    # mm / depth > 0.5  <=>  2*mm > depth
    consensus_diff = int(
        np.count_nonzero((depth > 0) & (mm * 2 > depth))
    )

    logging.debug(
        "Processed %s: %d reads, %d covbases, %d mismatches, "
        "%d consensus_diff, %d indels, %d mapped_bases, %d invalid_alns, "
        "%d linkage reads",
        rname,
        numreads,
        covbases,
        mismatches_total,
        consensus_diff,
        indels,
        mapped_bases,
        invalid_alns,
        len(linkage_reads),
    )

    metrics = [
        rname,
        start0,
        end0,
        numreads,
        covbases,
        mismatches_total,
        indels,
        consensus_diff,
        mapped_bases,
        invalid_alns,
        readlength,
    ]

    return metrics, species_idx, list(linkage_reads)


def _iter_tasks(
    references: List[str],
    lengths: List[int],
    reference_species_idx: List[Optional[int]],
    chunk_size: int,
) -> Iterable[Tuple[str, int, int, Optional[int]]]:
    for rname, rlen, species_idx in zip(
        references, lengths, reference_species_idx
    ):
        if rlen <= 0:
            continue
        cs = rlen if chunk_size <= 0 else chunk_size
        for start0 in range(0, rlen, cs):
            end0 = min(start0 + cs, rlen)
            yield (rname, start0, end0, species_idx)


def parse_aln_from_bam(
    bam_path: str,
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
    min_link_reads: int = 3,
    min_link_fraction: float = 0.01,
    same_genus_only: bool = False,
    max_link_species_per_read: int = 0,
) -> List[List]:
    """
    Process coverage and construct weighted linked-species groups in one BAM pass.

    Workers emit unique qualified read names for their reference chunks. The
    parent merges those observations into a unique species set per read. After
    the BAM pass, shared-read support is counted for every species pair and
    only sufficiently supported edges are unioned.
    """
    global _GROUPS, _GROUP_MEMBERS

    if min_link_reads < 1:
        raise ValueError("min_link_reads must be >= 1")
    if not 0.0 <= min_link_fraction <= 1.0:
        raise ValueError("min_link_fraction must be between 0 and 1")
    if max_link_species_per_read < 0:
        raise ValueError("max_link_species_per_read must be >= 0")

    if not os.path.exists(bam_path):
        print(f"ERROR: BAM not found: {bam_path}", file=sys.stderr)
        return []

    # Open BAM in the main process to validate the index and obtain references.
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if not bam.has_index():
                print(
                    "ERROR: BAM index (.bai) not found or not readable. "
                    "Pysam requires an index.",
                    file=sys.stderr,
                )
                return []
            references = list(bam.references)
            lengths = list(bam.lengths)
    except Exception as exc:
        print(f"ERROR: Failed to open BAM: {exc}", file=sys.stderr)
        return []

    logging.debug(
        "Parsing %d references with %d processes...",
        len(references),
        processes,
    )
    logging.debug(
        "Parameters: min_mapq=%s, min_frac=%s, min_idt=%s, min_alen=%s, "
        "include_secondary=%s, include_supplementary=%s, "
        "include_duplicates=%s, include_qcfail=%s, split_read_flag=%s, "
        "min_link_reads=%s, min_link_fraction=%s, same_genus_only=%s, "
        "max_link_species_per_read=%s",
        min_mapq,
        min_frac,
        min_idt,
        min_alen,
        include_secondary,
        include_supplementary,
        include_duplicates,
        include_qcfail,
        split_read_flag,
        min_link_reads,
        min_link_fraction,
        same_genus_only,
        max_link_species_per_read,
    )

    # Resolve each reference to species once in the parent. Tasks carry only a
    # compact integer species index, keeping taxonomy work out of the workers.
    species_to_idx: Dict[str, int] = {}
    idx_to_species: List[str] = []
    idx_to_genus: List[Optional[str]] = []
    reference_species_idx: List[Optional[int]] = []

    for rname in references:
        species_taxid = _rname_to_species_taxid(rname)

        if species_taxid is None:
            reference_species_idx.append(None)
            continue

        species_idx = species_to_idx.get(species_taxid)
        if species_idx is None:
            species_idx = len(idx_to_species)
            species_to_idx[species_taxid] = species_idx
            idx_to_species.append(species_taxid)
            idx_to_genus.append(
                _taxid_to_rank_taxid(species_taxid, target_rank="genus")
            )

        reference_species_idx.append(species_idx)

    logging.debug(
        "Resolved %d species across %d BAM references",
        len(idx_to_species),
        len(references),
    )

    # Parent-process read state. Most reads map to one species and are stored as
    # one integer; only multi-species reads are promoted to a set.
    read_species: Dict[str, ReadSpeciesState] = {}

    tasks = _iter_tasks(
        references,
        lengths,
        reference_species_idx,
        chunk_size,
    )

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

    n_chunk_read_species_observations = 0
    n_unique_read_species_observations = 0

    try:
        ref_chunk_results: List[List] = []
        header = [
            "RNAME",
            "STARTPOS",
            "ENDPOS",
            "NUMREADS",
            "COVBASES",
            "MISMATCHES",
            "INDELS",
            "CONSENSUS_DIFF",
            "MAPPED_BASES",
            "INVALID_ALNS",
            "READLENGTH",
        ]
        ref_chunk_results.append(header)

        mapper = pool.imap_unordered

        for metrics, species_idx, linkage_reads in mapper(
            _process_chunk, tasks, chunksize=imap_chunksize
        ):
            if species_idx is not None:
                n_chunk_read_species_observations += len(linkage_reads)

                for qname in linkage_reads:
                    if _add_read_species_observation(
                        read_species,
                        qname,
                        species_idx,
                    ):
                        n_unique_read_species_observations += 1

            # start0 -> one-based STARTPOS. END remains numerically correct.
            metrics[1] += 1
            ref_chunk_results.append(metrics)

    finally:
        pool.close()
        pool.join()

    logging.debug(
        "Total signature fragments processed: %d",
        len(ref_chunk_results) - 1,
    )

    (
        species_read_counts,
        edge_read_counts,
        read_degree_hist,
        n_multi_species_reads,
        n_excluded_hyperambiguous_reads,
    ) = _build_species_link_support(
        read_species,
        len(idx_to_species),
        max_link_species_per_read=max_link_species_per_read,
    )

    groups, species_to_group, edge_rows, n_component_unions = (
        _finalize_species_groups(
            idx_to_species=idx_to_species,
            idx_to_genus=idx_to_genus,
            species_read_counts=species_read_counts,
            edge_read_counts=edge_read_counts,
            min_link_reads=min_link_reads,
            min_link_fraction=min_link_fraction,
            same_genus_only=same_genus_only,
        )
    )

    _GROUPS = species_to_group
    _GROUP_MEMBERS = groups

    n_retained_edges = sum(1 for row in edge_rows if row[12])
    max_species_per_read = max(read_degree_hist, default=0)

    logging.info(
        "Species linkage from chunk pass: %d distinct reads, "
        "%d unique read-species observations, %d multi-species reads, "
        "%d observed species-pair edges, %d retained edges, "
        "%d component unions, %d linked species groups",
        len(read_species),
        n_unique_read_species_observations,
        n_multi_species_reads,
        len(edge_read_counts),
        n_retained_edges,
        n_component_unions,
        len(groups),
    )
    logging.info(
        "Species-link thresholds: min_link_reads=%d, min_link_fraction=%.6f, "
        "same_genus_only=%s, max_link_species_per_read=%d; "
        "max observed species/read=%d, excluded hyper-ambiguous reads=%d",
        min_link_reads,
        min_link_fraction,
        same_genus_only,
        max_link_species_per_read,
        max_species_per_read,
        n_excluded_hyperambiguous_reads,
    )
    logging.debug(
        "Chunk read-species observations before global deduplication: %d",
        n_chunk_read_species_observations,
    )
    logging.debug("Species to group mapping: %s", species_to_group)
    logging.debug("Group members: %s", groups)

    groups_out = f"{bam_path}.species_groups.tsv"
    links_out = f"{bam_path}.species_links.tsv"
    degree_out = f"{bam_path}.species_linkage_read_degree.tsv"

    write_species_groups(groups, groups_out)
    write_species_links(edge_rows, links_out)
    write_linkage_read_degree_histogram(read_degree_hist, degree_out)

    logging.info("Species groups written to %s", groups_out)
    logging.info("Species-link diagnostics written to %s", links_out)
    logging.info("Species/read degree histogram written to %s", degree_out)

    # Linkage state is no longer needed after diagnostics/groups are finalized.
    del read_species

    return ref_chunk_results


def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description="Compute coverage and consensus mismatch metrics from a BAM"
    )
    p.add_argument("bam", help="Input BAM path (requires .bai index).")
    p.add_argument("-o", "--out", required=True, help="Output TSV path.")
    p.add_argument("-c", "--chunk-size", type=int, default=1_000_000,
        help="Chunk size in reference bases for parallel tasks (default: 1,000,000). Use 0 for whole-contig.",
    )
    p.add_argument("-p", "--processes", type=int, default=max(1, mp.cpu_count() - 1),
        help="Worker processes (default: cpu_count-1).",
    )
    p.add_argument("-t", "--htslib-threads", type=int, default=1,
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


    # Species-link filtering
    p.add_argument("--min-link-reads", type=int, default=3,
        help=(
            "Minimum number of distinct shared reads required to retain a "
            "species-pair edge (default: 3)."
        ),
    )
    p.add_argument("--min-link-fraction", type=float, default=0.01,
        help=(
            "Minimum shared_reads / min(reads_species1, reads_species2) "
            "required to retain an edge (default: 0.01)."
        ),
    )
    p.add_argument("--same-genus-only", action="store_true",
        help=(
            "Retain species-link edges only when both species resolve to the "
            "same genus (default: off)."
        ),
    )
    p.add_argument(
        "--max-link-species-per-read", type=int, default=0,
        help=(
            "Exclude reads mapping to more than this many species from pair "
            "support. They still contribute to per-species read totals. "
            "Use 0 for no limit (default: 0)."
        ),
    )

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
        min_link_reads=args.min_link_reads,
        min_link_fraction=args.min_link_fraction,
        same_genus_only=args.same_genus_only,
        max_link_species_per_read=args.max_link_species_per_read,
    )
    
    out_path=args.out
    with open(out_path, "w", encoding="utf-8") as out:
        for res in ref_results:
            out.write("\t".join(map(str, res)) + "\n")

    return 0

if __name__ == "__main__":
    raise SystemExit(main())