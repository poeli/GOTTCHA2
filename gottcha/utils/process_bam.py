#!/usr/bin/env python3
"""
process_bam.py

Compute per-region coverage and consensus-mismatch metrics from an indexed BAM
without a reference FASTA. During the same chunk pass, optionally collect
query-segment-aware evidence linking alignments to different species.

Species linkage is delegated to ``species_linkage.py``. Primary and secondary
alignments that cover substantially overlapping intervals of the same query are
treated as competing explanations. Non-overlapping intervals of long reads do
not create species links, supplementary alignments are excluded by default,
and paired-end mates are kept separate even when they share a QNAME.

Alternative alignments are weighted by query overlap, shared aligned length,
and normalized AS-score difference. Evidence is aggregated over distinct reads,
query intervals, and signature pairs. Final groups are anchor-centered: every
member must link directly to the anchor, preventing weak transitive chains from
creating a giant connected component.

Coverage assumptions
--------------------
- minimap2 ``--eqx`` is expected; CIGAR X records mismatches and = matches.
- M/= /X contribute to depth; D and N do not.
- MISMATCHES is the sum of X bases.
- CONSENSUS_DIFF counts positions where X_depth / depth > 0.5.
"""

from __future__ import annotations

import argparse
import logging
import multiprocessing as mp
import os
import sys
from typing import Dict, Iterable, List, Optional, Tuple

import numpy as np
import pysam

try:
    from . import taxonomy as t
except ImportError:
    import taxonomy as t

try:
    from .species_linkage import (
        EvidenceConfig,
        GroupConfig,
        LinkageObservation,
        LinkageObservationStore,
        build_linkage_evidence,
        finalize_species_groups,
        make_read_key,
        observation_rank,
        query_interval_from_cigar,
        write_degree_histogram,
        write_species_evidence,
        write_species_groups,
        write_species_links,
    )
except ImportError:
    from species_linkage import (
        EvidenceConfig,
        GroupConfig,
        LinkageObservation,
        LinkageObservationStore,
        build_linkage_evidence,
        finalize_species_groups,
        make_read_key,
        observation_rank,
        query_interval_from_cigar,
        write_degree_histogram,
        write_species_evidence,
        write_species_groups,
        write_species_links,
    )


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

_GROUP_ANCHORS: Dict[int, str] = {}
_SHADOWS: Dict[str, str] = {}


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
    """Convert a reference name to a normalized species-level taxid."""

    try:
        taxid = rname.split("|")[-2]
        species_taxid = t.taxid2taxidOnRank(taxid, target_rank="species")
        if species_taxid is None:
            return None
        species_taxid = str(species_taxid)
        if not species_taxid or species_taxid == "0":
            return None
        return species_taxid
    except Exception as exc:
        logging.debug(
            "Unable to resolve species taxid from reference %s: %s",
            rname,
            exc,
        )
        return None


def _taxid_to_rank_taxid(taxid: str, target_rank: str) -> Optional[str]:
    """Resolve and normalize a taxid at the requested rank."""

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


# ----------------------------------------------------------------------
# Species linkage state
# ----------------------------------------------------------------------


def get_species_linkage_state() -> Dict[str, object]:
    """Return copies of linkage results from the latest call."""

    return {
        "species_to_group": dict(_GROUPS),
        "group_members": {
            group_id: list(members)
            for group_id, members in _GROUP_MEMBERS.items()
        },
        "group_anchors": dict(_GROUP_ANCHORS),
        "shadow_species": dict(_SHADOWS),
    }


# ----------------------------------------------------------------------
# Multiprocessing worker
# ----------------------------------------------------------------------

def _init_worker(
    bam_path: str,
    htslib_threads: int,
    coverage_cfg: Dict[str, object],
    linkage_cfg: Dict[str, object],
    split_read_flag: bool = False,
) -> None:
    """Open one BAM handle per worker and cache both filter sets."""

    global _BAM, _CFG
    _BAM = pysam.AlignmentFile(bam_path, "rb", threads=htslib_threads)
    _CFG = {
        "coverage": coverage_cfg,
        "linkage": linkage_cfg,
        "split_read_flag": split_read_flag,
    }


def _alignment_score_and_source(
    aln: pysam.AlignedSegment,
    aligned_length: int,
) -> Tuple[float, int]:
    """Return minimap2 AS or a conservative aligned-length-minus-NM proxy."""

    if aln.has_tag("AS"):
        return float(aln.get_tag("AS")), 1
    nm = float(aln.get_tag("NM")) if aln.has_tag("NM") else 0.0
    return float(aligned_length) - nm, 0


def _alignment_identity(
    aln: pysam.AlignedSegment,
    aligned_length: int,
) -> float:
    if aligned_length <= 0 or not aln.has_tag("NM"):
        return -1.0
    return max(0.0, 1.0 - float(aln.get_tag("NM")) / aligned_length)


def _linkage_read_key(aln: pysam.AlignedSegment) -> str:
    """Separate read groups and paired-end mates before interval clustering."""

    read_group = str(aln.get_tag("RG")) if aln.has_tag("RG") else ""
    if getattr(aln, "is_paired", False):
        if getattr(aln, "is_read1", False):
            mate_number = 1
        elif getattr(aln, "is_read2", False):
            mate_number = 2
        else:
            mate_number = 0
    else:
        mate_number = 0
    return make_read_key(aln.query_name, read_group, mate_number)


def _process_chunk(
    task: Tuple[str, int, int, Optional[int], int],
) -> Tuple[List, List[LinkageObservation]]:
    """Process one reference chunk and emit each linkage alignment once."""

    global _BAM, _CFG
    assert _BAM is not None, "Worker BAM handle not initialized"

    rname, start0, end0, species_idx, reference_idx = task
    length = end0 - start0
    if length <= 0:
        return [rname, start0, end0, 0, 0, 0, 0, 0, 0, 0, 0], []

    depth_diff = np.zeros(length + 1, dtype=np.int32)
    mismatch_diff = np.zeros(length + 1, dtype=np.int32)
    coverage_cfg = _CFG["coverage"]
    linkage_cfg = _CFG["linkage"]
    split_read_flag = bool(_CFG["split_read_flag"])

    numreads = 0
    readlength = 0
    indels = 0
    invalid_alns = 0
    bam = _BAM
    linkage_by_key: Dict[Tuple, LinkageObservation] = {}

    for aln in bam.fetch(rname, start0, end0):
        owns_alignment = start0 <= aln.reference_start < end0
        coverage_reason = _alignment_filter_reason(aln, bam, coverage_cfg)

        if (
            coverage_reason is not None
            and owns_alignment
            and coverage_reason in {"identity", "fraction", "length"}
        ):
            invalid_alns += 1

        # Linkage is intentionally independent of coverage inclusion flags.
        if owns_alignment and species_idx is not None and aln.query_name:
            linkage_reason = _alignment_filter_reason(aln, bam, linkage_cfg)
            if linkage_reason is None:
                query_interval = query_interval_from_cigar(
                    aln.cigartuples or (),
                    bool(aln.is_reverse),
                )
                if query_interval is not None:
                    qstart, qend = query_interval
                    aligned_length = max(
                        1,
                        int(
                            aln.query_alignment_length
                            or aln.alen
                            or (qend - qstart)
                        ),
                    )
                    score, has_as = _alignment_score_and_source(
                        aln,
                        aligned_length,
                    )
                    alignment_class = (
                        2
                        if aln.is_supplementary
                        else 1 if aln.is_secondary else 0
                    )
                    observation = LinkageObservation(
                        _linkage_read_key(aln),
                        qstart,
                        qend,
                        species_idx,
                        reference_idx,
                        score,
                        _alignment_identity(aln, aligned_length),
                        aligned_length,
                        alignment_class,
                        has_as,
                    )
                    key = (
                        observation.read_key,
                        observation.qstart,
                        observation.qend,
                        observation.species_idx,
                        observation.reference_idx,
                    )
                    current = linkage_by_key.get(key)
                    if (
                        current is None
                        or observation_rank(observation)
                        > observation_rank(current)
                    ):
                        linkage_by_key[key] = observation

        if coverage_reason is not None:
            continue

        if owns_alignment:
            if split_read_flag:
                if aln.has_tag("ZC"):
                    numreads += 1
            elif not aln.is_secondary and not aln.is_supplementary:
                numreads += 1
            readlength += int(aln.alen or 0)

        cigar = aln.cigartuples
        if not cigar:
            continue

        ref_pos = aln.reference_start
        block_start: Optional[int] = None
        for operation, operation_length in cigar:
            if operation_length <= 0:
                continue

            if operation in (0, 7, 8):  # M/= /X
                if block_start is None:
                    block_start = ref_pos

                if operation == 8:
                    segment_start = ref_pos
                    segment_end = ref_pos + operation_length
                    if segment_end > start0 and segment_start < end0:
                        segment_start = max(segment_start, start0)
                        segment_end = min(segment_end, end0)
                        mismatch_diff[segment_start - start0] += 1
                        mismatch_diff[segment_end - start0] -= 1
                ref_pos += operation_length

            elif operation in (2, 3):  # D/N
                if block_start is not None:
                    segment_start = block_start
                    segment_end = ref_pos
                    if segment_end > start0 and segment_start < end0:
                        segment_start = max(segment_start, start0)
                        segment_end = min(segment_end, end0)
                        depth_diff[segment_start - start0] += 1
                        depth_diff[segment_end - start0] -= 1
                        indels += operation_length
                    block_start = None
                ref_pos += operation_length

        if block_start is not None:
            segment_start = block_start
            segment_end = ref_pos
            if segment_end > start0 and segment_start < end0:
                segment_start = max(segment_start, start0)
                segment_end = min(segment_end, end0)
                depth_diff[segment_start - start0] += 1
                depth_diff[segment_end - start0] -= 1

    depth = np.cumsum(depth_diff[:-1], dtype=np.int32)
    mismatch_depth = np.cumsum(mismatch_diff[:-1], dtype=np.int32)
    covbases = int(np.count_nonzero(depth))
    mismatches_total = int(mismatch_depth.sum())
    mapped_bases = int(depth.sum()) + indels
    consensus_diff = int(
        np.count_nonzero((depth > 0) & (mismatch_depth * 2 > depth))
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
    return metrics, list(linkage_by_key.values())

def _iter_tasks(
    references: List[str],
    lengths: List[int],
    reference_species_idx: List[Optional[int]],
    chunk_size: int,
) -> Iterable[Tuple[str, int, int, Optional[int], int]]:
    for reference_idx, (rname, rlen, species_idx) in enumerate(
        zip(references, lengths, reference_species_idx)
    ):
        if rlen <= 0:
            continue
        current_chunk_size = rlen if chunk_size <= 0 else chunk_size
        for start0 in range(0, rlen, current_chunk_size):
            end0 = min(start0 + current_chunk_size, rlen)
            yield rname, start0, end0, species_idx, reference_idx

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
    link_min_mapq: Optional[int] = 0,
    link_min_frac: Optional[float] = None,
    link_min_idt: Optional[float] = None,
    link_min_alen: Optional[int] = None,
    link_include_secondary: bool = True,
    link_include_supplementary: bool = False,
    link_min_query_overlap: float = 0.80,
    link_max_species_per_segment: int = 5,
    link_score_tau: float = 0.05,
    link_length_saturation: int = 100,
    link_min_alignment_weight: float = 0.10,
    link_independent_score_margin: float = 0.05,
    min_link_reads: int = 3,
    min_link_segments: int = 3,
    min_link_loci: int = 2,
    min_link_weight: float = 1.5,
    same_genus_only: bool = True,
    min_shadow_containment: float = 0.80,
    min_mutual_ambiguity: float = 0.80,
    max_shadow_independent_fraction: float = 0.20,
    min_independent_loci_to_retain: int = 2,
    min_anchor_independent_loci_for_shadow: int = 2,
    linkage_storage: str = "sqlite",
    linkage_temp_dir: Optional[str] = None,
    keep_linkage_temp: bool = False,
    linkage_prefix: Optional[str] = None,
) -> List[List]:
    """Process coverage and query-segment-aware linkage in one chunk pass."""

    global _GROUPS, _GROUP_MEMBERS, _GROUP_ANCHORS, _SHADOWS

    if not os.path.exists(bam_path):
        print(f"ERROR: BAM not found: {bam_path}", file=sys.stderr)
        return []

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

    if link_min_frac is None:
        link_min_frac = min_frac
    if link_min_idt is None:
        link_min_idt = min_idt
    if link_min_alen is None:
        link_min_alen = max(int(min_alen), 60)
    if link_min_mapq is None:
        link_min_mapq = 0

    coverage_cfg: Dict[str, object] = {
        "min_mapq": int(min_mapq or 0),
        "min_frac": float(min_frac),
        "min_idt": float(min_idt),
        "min_alen": int(min_alen),
        "include_secondary": bool(include_secondary),
        "include_supplementary": bool(include_supplementary),
        "include_duplicates": bool(include_duplicates),
        "include_qcfail": bool(include_qcfail),
    }
    linkage_cfg: Dict[str, object] = {
        "min_mapq": int(link_min_mapq),
        "min_frac": float(link_min_frac),
        "min_idt": float(link_min_idt),
        "min_alen": int(link_min_alen),
        "include_secondary": bool(link_include_secondary),
        "include_supplementary": bool(link_include_supplementary),
        "include_duplicates": bool(include_duplicates),
        "include_qcfail": bool(include_qcfail),
    }

    evidence_config = EvidenceConfig(
        min_query_overlap=link_min_query_overlap,
        max_species_per_segment=link_max_species_per_segment,
        score_tau=link_score_tau,
        length_saturation=link_length_saturation,
        min_alignment_weight=link_min_alignment_weight,
        independent_score_margin=link_independent_score_margin,
    )
    group_config = GroupConfig(
        min_link_reads=min_link_reads,
        min_link_segments=min_link_segments,
        min_link_loci=min_link_loci,
        min_link_weight=min_link_weight,
        same_genus_only=same_genus_only,
        min_shadow_containment=min_shadow_containment,
        min_mutual_ambiguity=min_mutual_ambiguity,
        max_shadow_independent_fraction=max_shadow_independent_fraction,
        min_independent_loci_to_retain=min_independent_loci_to_retain,
        min_anchor_independent_loci_for_shadow=(
            min_anchor_independent_loci_for_shadow
        ),
    )
    evidence_config.validate()
    group_config.validate()

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

    tasks = _iter_tasks(
        references,
        lengths,
        reference_species_idx,
        int(chunk_size or 0),
    )
    store = LinkageObservationStore(
        mode=linkage_storage,
        temp_dir=linkage_temp_dir,
        keep=keep_linkage_temp,
    )

    ref_chunk_results: List[List] = [[
        "RNAME", "STARTPOS", "ENDPOS", "NUMREADS", "COVBASES",
        "MISMATCHES", "INDELS", "CONSENSUS_DIFF", "MAPPED_BASES",
        "INVALID_ALNS", "READLENGTH",
    ]]

    try:
        pool = mp.Pool(
            processes=processes,
            initializer=_init_worker,
            initargs=(
                bam_path,
                int(htslib_threads or 1),
                coverage_cfg,
                linkage_cfg,
                bool(split_read_flag),
            ),
        )
        try:
            for metrics, observations in pool.imap_unordered(
                _process_chunk,
                tasks,
                chunksize=int(imap_chunksize or 1),
            ):
                store.add_many(observations)
                metrics[1] += 1
                ref_chunk_results.append(metrics)
        finally:
            pool.close()
            pool.join()

        logging.info(
            "Coverage pass: %d chunks; %d linkage observations",
            len(ref_chunk_results) - 1,
            store.n_observations,
        )
        evidence = build_linkage_evidence(
            store.iter_reads(),
            n_species=len(idx_to_species),
            config=evidence_config,
        )
    finally:
        sqlite_path = store.path
        store.close()

    groups = finalize_species_groups(
        idx_to_species,
        idx_to_genus,
        evidence,
        config=group_config,
    )

    _GROUPS = groups.species_to_group
    _GROUP_MEMBERS = groups.groups
    _GROUP_ANCHORS = {
        group_id: str(meta["anchor_taxid"])
        for group_id, meta in groups.group_meta.items()
    }
    _SHADOWS = groups.shadows

    output_prefix = linkage_prefix or str(bam_path)
    groups_out = f"{output_prefix}.species_groups.tsv"
    links_out = f"{output_prefix}.species_links.tsv"
    evidence_out = f"{output_prefix}.species_evidence.tsv"
    read_degree_out = f"{output_prefix}.species_linkage_read_degree.tsv"
    segment_degree_out = f"{output_prefix}.species_linkage_segment_degree.tsv"

    write_species_groups(groups, groups_out)
    write_species_links(groups.edge_rows, links_out)
    write_species_evidence(
        idx_to_species,
        idx_to_genus,
        evidence,
        groups,
        evidence_out,
    )
    write_degree_histogram(
        evidence.read_degree_histogram,
        read_degree_out,
        "N_SPECIES_PER_READ",
    )
    write_degree_histogram(
        evidence.segment_degree_histogram,
        segment_degree_out,
        "N_SPECIES_PER_QUERY_SEGMENT",
    )

    logging.info(
        "Species linkage: %d reads, %d query segments, %d pairs, "
        "%d groups, %d conservative shadow candidates",
        evidence.stats["reads"],
        evidence.stats["query_segments"],
        len(evidence.pair_evidence),
        len(groups.groups),
        len(groups.shadows),
    )
    logging.info(
        "Linkage outputs: %s, %s, %s",
        groups_out,
        links_out,
        evidence_out,
    )
    if keep_linkage_temp and sqlite_path:
        logging.info("Linkage SQLite retained at %s", sqlite_path)
    if evidence.stats["fallback_score_observations"]:
        logging.warning(
            "%d linkage observations lacked AS and used aligned_length-NM",
            evidence.stats["fallback_score_observations"],
        )

    return ref_chunk_results

def main(argv: Optional[List[str]] = None) -> int:
    p = argparse.ArgumentParser(
        description=(
            "Compute coverage and query-segment-aware species linkage from BAM."
        )
    )
    p.add_argument("bam", help="Input BAM path (requires .bai index).")
    p.add_argument("-o", "--out", required=True, help="Output coverage TSV.")
    p.add_argument("-c", "--chunk-size", type=int, default=1_000_000)
    p.add_argument(
        "-p", "--processes", type=int,
        default=max(1, mp.cpu_count() - 1),
    )
    p.add_argument("-t", "--htslib-threads", type=int, default=1)

    # Coverage filters.
    p.add_argument("--min-mapq", type=int, default=0)
    p.add_argument("--min-frac", type=float, default=0.0)
    p.add_argument("--min-idt", type=float, default=0.0)
    p.add_argument("--min-alen", type=int, default=0)
    p.add_argument("--include-secondary", action="store_true")
    p.add_argument("--include-supplementary", action="store_true")
    p.add_argument("--include-duplicates", action="store_true")
    p.add_argument("--include-qcfail", action="store_true")

    # Linkage alignment selection, independent of coverage flags.
    p.set_defaults(link_include_secondary=True, same_genus_only=True)
    p.add_argument("--link-min-mapq", type=int, default=0)
    p.add_argument("--link-min-frac", type=float, default=None)
    p.add_argument("--link-min-idt", type=float, default=None)
    p.add_argument("--link-min-alen", type=int, default=None)
    p.add_argument(
        "--link-primary-only",
        dest="link_include_secondary",
        action="store_false",
    )
    p.add_argument(
        "--link-include-secondary",
        dest="link_include_secondary",
        action="store_true",
    )
    p.add_argument("--link-include-supplementary", action="store_true")
    p.add_argument("--link-min-query-overlap", type=float, default=0.80)
    p.add_argument("--link-max-species-per-segment", type=int, default=5)
    p.add_argument("--link-score-tau", type=float, default=0.05)
    p.add_argument("--link-length-saturation", type=int, default=100)
    p.add_argument("--link-min-alignment-weight", type=float, default=0.10)
    p.add_argument(
        "--link-independent-score-margin", type=float, default=0.05
    )

    # Aggregated edge/group filters.
    p.add_argument("--min-link-reads", type=int, default=3)
    p.add_argument("--min-link-segments", type=int, default=3)
    p.add_argument("--min-link-loci", type=int, default=2)
    p.add_argument("--min-link-weight", type=float, default=1.5)
    p.add_argument(
        "--same-genus-only",
        dest="same_genus_only",
        action="store_true",
    )
    p.add_argument(
        "--allow-cross-genus-links",
        dest="same_genus_only",
        action="store_false",
    )
    p.add_argument("--min-shadow-containment", type=float, default=0.80)
    p.add_argument("--min-mutual-ambiguity", type=float, default=0.80)
    p.add_argument(
        "--max-shadow-independent-fraction", type=float, default=0.20
    )
    p.add_argument("--min-independent-loci-to-retain", type=int, default=2)
    p.add_argument(
        "--min-anchor-independent-loci-for-shadow", type=int, default=2
    )

    # Storage and diagnostic outputs.
    p.add_argument(
        "--linkage-storage",
        choices=("memory", "sqlite"),
        default="sqlite",
    )
    p.add_argument("--linkage-temp-dir", default=None)
    p.add_argument("--keep-linkage-temp", action="store_true")
    p.add_argument("--linkage-prefix", default=None)
    p.add_argument("--imap-chunksize", type=int, default=1)

    args = p.parse_args(argv)
    results = parse_aln_from_bam(
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
        link_min_mapq=args.link_min_mapq,
        link_min_frac=args.link_min_frac,
        link_min_idt=args.link_min_idt,
        link_min_alen=args.link_min_alen,
        link_include_secondary=args.link_include_secondary,
        link_include_supplementary=args.link_include_supplementary,
        link_min_query_overlap=args.link_min_query_overlap,
        link_max_species_per_segment=args.link_max_species_per_segment,
        link_score_tau=args.link_score_tau,
        link_length_saturation=args.link_length_saturation,
        link_min_alignment_weight=args.link_min_alignment_weight,
        link_independent_score_margin=args.link_independent_score_margin,
        min_link_reads=args.min_link_reads,
        min_link_segments=args.min_link_segments,
        min_link_loci=args.min_link_loci,
        min_link_weight=args.min_link_weight,
        same_genus_only=args.same_genus_only,
        min_shadow_containment=args.min_shadow_containment,
        min_mutual_ambiguity=args.min_mutual_ambiguity,
        max_shadow_independent_fraction=(
            args.max_shadow_independent_fraction
        ),
        min_independent_loci_to_retain=args.min_independent_loci_to_retain,
        min_anchor_independent_loci_for_shadow=(
            args.min_anchor_independent_loci_for_shadow
        ),
        linkage_storage=args.linkage_storage,
        linkage_temp_dir=args.linkage_temp_dir,
        keep_linkage_temp=args.keep_linkage_temp,
        linkage_prefix=args.linkage_prefix,
    )
    if not results:
        return 2

    with open(args.out, "w", encoding="utf-8") as out:
        for result in results:
            out.write("\t".join(map(str, result)) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
