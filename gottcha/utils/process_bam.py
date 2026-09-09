#!/usr/bin/env python3
"""
process_bam.py

Compute per-region coverage and consensus-mismatch metrics from a BAM file
(with .bai index present), **without** a reference FASTA.

Assumptions / notes:
- Alignments are from minimap2 with `--eqx`, so mismatches are encoded as CIGAR
  op `X` and matches as `=`.
- Depth is computed from aligned query bases (CIGAR ops M/= /X). Deletions (D)
  and refskips (N) do not contribute to depth in this implementation.
- "mismatches" counts total mismatched aligned bases across all reads (sum of
  `X`).
- "pileup_mismatch" counts reference positions where >50% of aligned reads at
  that position are mismatches (i.e. #positions where X_depth / depth > 0.5).

Species linkage:
- Each BAM reference is resolved to a species-level taxid before chunking.
- During the normal parallel BAM pass, each worker records qualified primary
  and secondary/supplementary read observations for the chunk's species.
- Because the BAM is coordinate-sorted, alignments for one read may occur on
  different references/chunks. Workers therefore do not construct final
  cross-species groups independently.
- The parent process merges the read observations into directed links from the
  primary species to each distinct secondary/supplementary species for that
  read, counts distinct supporting reads per directed species pair, and then
  identifies strongly connected components (SCCs).
- Only SCCs containing more than one species are reported as linked groups.

Parallelization:
- References are split into fixed-size chunks along their length. Each chunk is
  processed independently in a worker process.
- Each worker opens the BAM once (via Pool initializer) for performance.
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import os
import sys
from collections import Counter, defaultdict
from typing import Counter as CounterType
from typing import Dict, Iterable, List, Optional, Set, Tuple, Union

import numpy as np
import pysam

try:
    from . import taxonomy as t
except ImportError:
    import taxonomy as t


# Global BAM handle and config for worker processes
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

# A read key distinguishes paired-end mates that share the same QNAME.
ReadKey = Tuple[str, int]
DirectedLink = Tuple[str, str]
ChunkTask = Tuple[str, int, int, Optional[str]]
ChunkResult = Tuple[List, Optional[str], Set[ReadKey], Set[ReadKey]]


def _rname_to_taxid_at_rank(
    rname: str,
    target_rank: str,
) -> Optional[str]:
    """Convert a BAM reference name to a taxid at the requested rank."""

    try:
        taxid = rname.split("|")[-2]
        resolved = t.taxid2taxidOnRank(taxid, target_rank=target_rank)
        if resolved in (None, "", "0", "unknown"):
            return None
        return str(resolved)
    except Exception as exc:
        logging.debug(
            "Unable to resolve taxid at rank %s from reference %s: %s",
            target_rank,
            rname,
            exc,
        )
        return None


def _taxid_sort_key(taxid: str) -> Tuple[int, Union[int, str]]:
    """Sort numeric taxids numerically and other identifiers lexicographically."""

    try:
        return (0, int(taxid))
    except (TypeError, ValueError):
        return (1, str(taxid))


def _read_key(aln: pysam.AlignedSegment) -> Optional[ReadKey]:
    """Return a read identifier that keeps paired-end mates separate."""

    qname = aln.query_name
    if not qname:
        return None

    if aln.is_read1:
        mate = 1
    elif aln.is_read2:
        mate = 2
    else:
        mate = 0

    return (qname, mate)


def _match_filter_reason(
    aln: pysam.AlignedSegment,
    bam: pysam.AlignmentFile,
    rname: str,
    min_frac: float,
    min_idt: float,
    min_alen: int,
) -> Optional[str]:
    """
    Return the matching-threshold rejection reason, or None when qualified.

    This contains the same identity/fraction/aligned-length criteria used by
    the existing chunk parser. Role-specific inclusion of secondary and
    supplementary alignments is intentionally handled outside this helper.
    """

    alen = aln.alen or 0

    if min_idt > 0.0 and aln.has_tag("NM"):
        if alen <= 0:
            return "identity"
        mm_idt = 1.0 - (aln.get_tag("NM") / alen)
        if min_idt > mm_idt:
            return "identity"

    if min_frac > 0.0:
        query_length = aln.query_length
        if not query_length:
            # For hard clips, query_length can be 0; recover it from CIGAR.
            query_length = sum(
                length
                for op, length in (aln.cigartuples or ())
                if op in (0, 1, 4, 5, 7, 8)  # M/I/S/H/=/X
            )

        if query_length <= 0:
            return "fraction"

        reference_length = bam.get_reference_length(rname)
        if reference_length <= 0:
            return "fraction"

        if (
            (alen / query_length) < min_frac
            and (alen / reference_length) < min_frac
        ):
            return "fraction"

    if min_alen > 0 and alen < min_alen:
        return "length"

    return None


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


def _process_chunk(task: ChunkTask) -> ChunkResult:
    """
    Process one ``(rname, start0, end0, species_taxid)`` chunk.

    Returns
    -------
    metrics
        The existing per-chunk coverage/mismatch result row.
    species_taxid
        Species-level taxid for this reference, or None if unresolved.
    primary_reads
        Unique qualified primary read keys whose alignments START in this
        chunk.
    alternative_reads
        Unique qualified secondary/supplementary read keys whose alignments
        START in this chunk.

    Notes
    -----
    Linkage observations are emitted only by the chunk that owns an alignment
    (the alignment starts in that chunk), preventing duplicate observations for
    alignments that span chunk boundaries.

    Secondary/supplementary alignments are always considered for linkage. The
    existing ``include_secondary`` / ``include_supplementary`` settings retain
    their original meaning for coverage/mismatch accumulation only.
    """

    global _BAM, _CFG
    assert _BAM is not None, "Worker BAM handle not initialized"

    rname, start0, end0, species_taxid = task
    L = end0 - start0
    if L <= 0:
        metrics = [rname, start0, end0, 0, 0, 0, 0, 0, 0, 0, 0]
        return metrics, species_taxid, set(), set()

    # Difference arrays (signed) so we can do O(segments) updates and O(L)
    # cumsums.
    # depth[pos] = #reads with an aligned base at that position (M/= /X)
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

    # One chunk belongs to one BAM reference/species. Sets deduplicate repeated
    # alignments of the same read to that species within this chunk.
    primary_reads: Set[ReadKey] = set()
    alternative_reads: Set[ReadKey] = set()

    # Iterate reads overlapping this region.
    for aln in bam.fetch(rname, start0, end0):
        if aln.is_unmapped:
            continue
        if (not inc_dup) and aln.is_duplicate:
            continue
        if (not inc_qcf) and aln.is_qcfail:
            continue
        if aln.mapping_quality < min_mapq:
            continue

        # Only the chunk in which the alignment begins should emit a linkage
        # observation or count it as a new alignment/read.
        owns_alignment = aln.reference_start >= start0

        match_filter_reason: Optional[str] = None
        if owns_alignment:
            match_filter_reason = _match_filter_reason(
                aln,
                bam,
                rname,
                min_frac,
                min_idt,
                min_alen,
            )

            # Linkage qualification is independent of whether secondary or
            # supplementary alignments are included in coverage metrics.
            if match_filter_reason is None and species_taxid is not None:
                key = _read_key(aln)
                if key is not None:
                    if aln.is_secondary or aln.is_supplementary:
                        alternative_reads.add(key)
                    else:
                        primary_reads.add(key)

        # Preserve the existing role-specific coverage behavior.
        if (not inc_sec) and aln.is_secondary:
            continue
        if (not inc_sup) and aln.is_supplementary:
            continue

        # The original implementation evaluates min-idt/min-frac/min-alen only
        # in the chunk where the alignment starts. Preserve that behavior for
        # coverage/mismatch aggregation.
        if owns_alignment:
            if match_filter_reason is not None:
                invalid_alns += 1
                continue

            # If split_read_flag is set, only count reads with ZC tag (the first
            # chunked reads) towards numreads.
            if split_read_flag:
                if aln.has_tag("ZC"):
                    numreads += 1
            else:
                if aln.is_secondary or aln.is_supplementary:
                    pass
                else:
                    numreads += 1

            # Count aligned read length for mean-depth/SNI calculations.
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

            elif op in (2, 3):  # D or N: consumes reference, not query
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
                # We do not break the block on insertions/softclips because
                # reference positions remain contiguous.
                continue

        # Close any remaining aligned block.
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

    # Build per-base depth and mismatch arrays.
    depth = np.cumsum(depth_diff[:-1], dtype=np.int32)
    mm = np.cumsum(mm_diff[:-1], dtype=np.int32)

    covbases = int(np.count_nonzero(depth))
    mismatches_total = int(mm.sum())
    mapped_bases = int(depth.sum()) + indels
    consensus_diff = int(np.count_nonzero((depth > 0) & (mm * 2 > depth)))

    logging.debug(
        "Processed %s: %d reads, %d covbases, %d mismatches, %d consensus_diff, "
        "%d indels, %d mapped_bases, %d invalid_alns, %d primary linkage reads, "
        "%d alternative linkage reads",
        rname,
        numreads,
        covbases,
        mismatches_total,
        consensus_diff,
        indels,
        mapped_bases,
        invalid_alns,
        len(primary_reads),
        len(alternative_reads),
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
    return metrics, species_taxid, primary_reads, alternative_reads


def _iter_tasks(
    references: List[str],
    lengths: List[int],
    species_taxids: List[Optional[str]],
    chunk_size: int,
) -> Iterable[ChunkTask]:
    """Yield chunk tasks with the pre-resolved species taxid for each reference."""

    for rname, rlen, species_taxid in zip(references, lengths, species_taxids):
        if rlen <= 0:
            continue
        cs = rlen if chunk_size <= 0 else chunk_size
        for start0 in range(0, rlen, cs):
            end0 = min(start0 + cs, rlen)
            yield (rname, start0, end0, species_taxid)


def _add_directed_link_for_read(
    read_key: ReadKey,
    primary_species: str,
    alternative_species: str,
    link_counts: CounterType[DirectedLink],
    linked_targets_by_read: Dict[ReadKey, Set[str]],
) -> None:
    """Add one distinct-read support observation for a directed species link."""

    if primary_species == alternative_species:
        return

    linked_targets = linked_targets_by_read.setdefault(read_key, set())
    if alternative_species in linked_targets:
        return

    linked_targets.add(alternative_species)
    link_counts[(primary_species, alternative_species)] += 1


def _merge_linkage_observations(
    species_taxid: Optional[str],
    primary_reads: Set[ReadKey],
    alternative_reads: Set[ReadKey],
    primary_by_read: Dict[ReadKey, Optional[str]],
    pending_alternatives_by_read: Dict[ReadKey, Set[str]],
    linked_targets_by_read: Dict[ReadKey, Set[str]],
    link_counts: CounterType[DirectedLink],
) -> int:
    """
    Merge one worker's read-role observations into global directed link counts.

    Returns the number of read keys found to have conflicting primary species.
    Such reads are excluded from linkage rather than allowing an ambiguous
    primary assignment to create spurious directed edges.
    """

    if species_taxid is None:
        return 0

    conflicting_primary_reads = 0

    for read_key in primary_reads:
        if read_key not in primary_by_read:
            primary_by_read[read_key] = species_taxid

            pending = pending_alternatives_by_read.pop(read_key, None)
            if pending:
                for alternative_species in pending:
                    _add_directed_link_for_read(
                        read_key,
                        species_taxid,
                        alternative_species,
                        link_counts,
                        linked_targets_by_read,
                    )
            continue

        previous_primary = primary_by_read[read_key]

        # None marks a read already known to have conflicting primaries.
        if previous_primary is None:
            continue

        if previous_primary == species_taxid:
            continue

        # Conflicting primaries should be rare in a well-formed BAM. If one is
        # discovered after links were already emitted, retract those supports.
        for alternative_species in linked_targets_by_read.pop(read_key, set()):
            edge = (previous_primary, alternative_species)
            link_counts[edge] -= 1
            if link_counts[edge] <= 0:
                del link_counts[edge]

        pending_alternatives_by_read.pop(read_key, None)
        primary_by_read[read_key] = None
        conflicting_primary_reads += 1

    for read_key in alternative_reads:
        if read_key not in primary_by_read:
            pending_alternatives_by_read[read_key].add(species_taxid)
            continue

        primary_species = primary_by_read[read_key]
        if primary_species is None:
            continue

        _add_directed_link_for_read(
            read_key,
            primary_species,
            species_taxid,
            link_counts,
            linked_targets_by_read,
        )

    return conflicting_primary_reads


def _strongly_connected_components(
    link_counts: CounterType[DirectedLink],
) -> List[List[str]]:
    """Return deterministic SCCs from directed links using iterative Kosaraju."""

    graph: Dict[str, Set[str]] = defaultdict(set)
    reverse_graph: Dict[str, Set[str]] = defaultdict(set)
    nodes: Set[str] = set()

    for (source, target), count in link_counts.items():
        if count <= 0 or source == target:
            continue
        graph[source].add(target)
        reverse_graph[target].add(source)
        nodes.add(source)
        nodes.add(target)

    if not nodes:
        return []

    # First pass: finishing order on the original directed graph.
    visited: Set[str] = set()
    finish_order: List[str] = []

    for start in sorted(nodes, key=_taxid_sort_key):
        if start in visited:
            continue

        stack: List[Tuple[str, bool]] = [(start, False)]
        while stack:
            node, expanded = stack.pop()

            if expanded:
                finish_order.append(node)
                continue

            if node in visited:
                continue

            visited.add(node)
            stack.append((node, True))

            # Reverse sort because stack is LIFO; this keeps traversal stable.
            for neighbor in sorted(
                graph.get(node, ()),
                key=_taxid_sort_key,
                reverse=True,
            ):
                if neighbor not in visited:
                    stack.append((neighbor, False))

    # Second pass: traverse the transpose in reverse finishing order.
    visited.clear()
    components: List[List[str]] = []

    for start in reversed(finish_order):
        if start in visited:
            continue

        component: List[str] = []
        stack = [start]
        visited.add(start)

        while stack:
            node = stack.pop()
            component.append(node)

            for neighbor in reverse_graph.get(node, ()):
                if neighbor not in visited:
                    visited.add(neighbor)
                    stack.append(neighbor)

        component.sort(key=_taxid_sort_key)
        components.append(component)

    components.sort(
        key=lambda members: tuple(_taxid_sort_key(taxid) for taxid in members)
    )
    return components


def _build_species_groups(
    link_counts: CounterType[DirectedLink],
) -> Dict[int, List[str]]:
    """Build linked groups from SCCs and populate the module-level lookups."""

    global _GROUPS, _GROUP_MEMBERS

    _GROUPS.clear()
    _GROUP_MEMBERS.clear()

    linked_components = [
        component
        for component in _strongly_connected_components(link_counts)
        if len(component) > 1
    ]

    groups: Dict[int, List[str]] = {}
    for group_id, members in enumerate(linked_components, start=1):
        groups[group_id] = members
        _GROUP_MEMBERS[group_id] = list(members)
        for species_taxid in members:
            _GROUPS[species_taxid] = group_id

    return groups


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
) -> Union[int, Tuple[List[List], Dict[int, List[str]]]]:
    """
    Parse BAM chunks and identify strongly connected species groups.

    On success returns ``(ref_chunk_results, groups)``. The existing per-chunk
    result rows are unchanged; ``groups`` maps group id to species-level taxids.
    """

    _GROUPS.clear()
    _GROUP_MEMBERS.clear()

    if not os.path.exists(bam_path):
        print(f"ERROR: BAM not found: {bam_path}", file=sys.stderr)
        return 2

    # Open BAM in main process to validate index and obtain reference lengths.
    try:
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            if not bam.has_index():
                print(
                    "ERROR: BAM index (.bai) not found or not readable. "
                    "Pysam requires an index.",
                    file=sys.stderr,
                )
                return 2
            references = list(bam.references)
            lengths = list(bam.lengths)
    except Exception as exc:
        print(f"ERROR: Failed to open BAM: {exc}", file=sys.stderr)
        return 2

    # Resolve each reference once in the main process. This avoids repeated
    # taxonomy traversal in workers/chunks and works cleanly with spawn-based
    # multiprocessing as well as fork-based multiprocessing.
    species_taxids = [
        _rname_to_taxid_at_rank(rname, target_rank="species")
        for rname in references
    ]
    unresolved_references = sum(taxid is None for taxid in species_taxids)
    if unresolved_references:
        logging.debug(
            "Unable to resolve %d/%d BAM references to species-level taxids",
            unresolved_references,
            len(references),
        )

    logging.debug(
        "Parsing %d references with %d processes...",
        len(references),
        processes,
    )
    logging.debug(
        "Parameters: min_mapq=%s, min_frac=%s, min_idt=%s, min_alen=%s, "
        "include_secondary=%s, include_supplementary=%s, include_duplicates=%s, "
        "include_qcfail=%s, split_read_flag=%s",
        min_mapq,
        min_frac,
        min_idt,
        min_alen,
        include_secondary,
        include_supplementary,
        include_duplicates,
        include_qcfail,
        split_read_flag,
    )

    tasks = _iter_tasks(references, lengths, species_taxids, chunk_size)

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

    # Main-process linkage state. Because a coordinate-sorted BAM can place a
    # read's primary and alternative alignments on unrelated chunks/references,
    # the read key must be retained until those observations are reunited.
    primary_by_read: Dict[ReadKey, Optional[str]] = {}
    pending_alternatives_by_read: Dict[ReadKey, Set[str]] = defaultdict(set)
    linked_targets_by_read: Dict[ReadKey, Set[str]] = {}
    link_counts: CounterType[DirectedLink] = Counter()
    conflicting_primary_reads = 0

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
        for (
            result,
            species_taxid,
            primary_reads,
            alternative_reads,
        ) in mapper(_process_chunk, tasks, chunksize=imap_chunksize):
            result[1] += 1  # 0-based start0 -> 1-based STARTPOS
            ref_chunk_results.append(result)

            conflicting_primary_reads += _merge_linkage_observations(
                species_taxid,
                primary_reads,
                alternative_reads,
                primary_by_read,
                pending_alternatives_by_read,
                linked_targets_by_read,
                link_counts,
            )

    finally:
        pool.close()
        pool.join()

    groups = _build_species_groups(link_counts)

    logging.debug(
        "Total signature fragments processed: %d",
        len(ref_chunk_results) - 1,
    )
    logging.info(
        "Species linkage: %d directed species pairs supported by %d distinct "
        "read-pair observations; %d SCC groups; %d reads with conflicting "
        "primary species; %d alternative-only reads without an observed primary",
        len(link_counts),
        sum(link_counts.values()),
        len(groups),
        conflicting_primary_reads,
        len(pending_alternatives_by_read),
    )
    logging.debug("Species directed-link counts: %s", dict(link_counts))
    logging.debug("Species groups: %s", groups)

    return ref_chunk_results, groups


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
        help=(
            "Chunk size in reference bases for parallel tasks (default: "
            "1,000,000). Use 0 for whole-contig."
        ),
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
    p.add_argument(
        "--min-mapq",
        type=int,
        default=0,
        help="Minimum MAPQ to keep an alignment (default: 0).",
    )
    p.add_argument(
        "--min-frac",
        type=float,
        default=0.0,
        help="Minimum fraction to keep an alignment (default: 0.0).",
    )
    p.add_argument(
        "--min-idt",
        type=float,
        default=0.0,
        help="Minimum identity to keep an alignment (default: 0.0).",
    )
    p.add_argument(
        "--min-alen",
        type=int,
        default=0,
        help="Minimum alignment length to keep an alignment (default: 0).",
    )
    p.add_argument(
        "--include-secondary",
        action="store_true",
        help="Include secondary alignments in coverage metrics (default: off).",
    )
    p.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Include supplementary alignments in coverage metrics (default: off).",
    )
    p.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include duplicate-marked reads (default: off).",
    )
    p.add_argument(
        "--include-qcfail",
        action="store_true",
        help="Include QC-failed reads (default: off).",
    )
    p.add_argument(
        "--imap-chunksize",
        type=int,
        default=1,
        help="chunksize passed to multiprocessing imap/imap_unordered (default: 1).",
    )

    args = p.parse_args(argv)

    parsed = parse_aln_from_bam(
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

    if parsed == 2:
        return 2

    ref_results, groups = parsed

    with open(args.out, "w", encoding="utf-8") as out:
        for res in ref_results:
            out.write("\t".join(map(str, res)) + "\n")

    if groups:
        logging.info("Identified linked species groups: %s", groups)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
