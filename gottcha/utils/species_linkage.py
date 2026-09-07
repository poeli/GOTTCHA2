#!/usr/bin/env python3
"""Query-segment-aware evidence model for cross-species read alignments.

This module is deliberately independent of pysam. ``process_bam.py`` converts
qualified BAM alignments into :class:`LinkageObservation` records, while this
module clusters homologous query intervals, accumulates weighted evidence, and
forms conservative anchor-centered species groups.
"""

from __future__ import annotations

import math
import os
import sqlite3
import tempfile
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from itertools import groupby
from typing import (
    Dict,
    Iterable,
    Iterator,
    List,
    NamedTuple,
    Optional,
    Sequence,
    Set,
    Tuple,
)


class LinkageObservation(NamedTuple):
    """Compact alignment observation transferred from a BAM worker."""

    read_key: str
    qstart: int
    qend: int
    species_idx: int
    reference_idx: int
    score: float
    identity: float
    aligned_length: int
    alignment_class: int  # 0=primary, 1=secondary, 2=supplementary
    has_as: int


@dataclass(frozen=True)
class EvidenceConfig:
    """Controls query-interval clustering and per-alignment weighting."""

    min_query_overlap: float = 0.80
    max_species_per_segment: int = 5
    score_tau: float = 0.05
    length_saturation: int = 100
    min_alignment_weight: float = 0.10
    independent_score_margin: float = 0.05

    def validate(self) -> None:
        if not 0.0 <= self.min_query_overlap <= 1.0:
            raise ValueError("min_query_overlap must be between 0 and 1")
        if self.max_species_per_segment < 0:
            raise ValueError("max_species_per_segment must be >= 0")
        if self.score_tau <= 0.0:
            raise ValueError("score_tau must be > 0")
        if self.length_saturation < 0:
            raise ValueError("length_saturation must be >= 0")
        if self.min_alignment_weight < 0.0:
            raise ValueError("min_alignment_weight must be >= 0")
        if self.independent_score_margin < 0.0:
            raise ValueError("independent_score_margin must be >= 0")


@dataclass(frozen=True)
class GroupConfig:
    """Controls edge retention, grouping, and conservative shadow labels."""

    min_link_reads: int = 3
    min_link_segments: int = 3
    min_link_loci: int = 2
    min_link_weight: float = 1.5
    same_genus_only: bool = True
    min_shadow_containment: float = 0.80
    min_mutual_ambiguity: float = 0.80
    max_shadow_independent_fraction: float = 0.20
    min_independent_loci_to_retain: int = 2
    min_anchor_independent_loci_for_shadow: int = 2

    def validate(self) -> None:
        if min(
            self.min_link_reads,
            self.min_link_segments,
            self.min_link_loci,
            self.min_independent_loci_to_retain,
            self.min_anchor_independent_loci_for_shadow,
        ) < 0:
            raise ValueError("count thresholds must be >= 0")
        if self.min_link_weight < 0.0:
            raise ValueError("min_link_weight must be >= 0")
        for name, value in (
            ("min_shadow_containment", self.min_shadow_containment),
            ("min_mutual_ambiguity", self.min_mutual_ambiguity),
            (
                "max_shadow_independent_fraction",
                self.max_shadow_independent_fraction,
            ),
        ):
            if not 0.0 <= value <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")


@dataclass
class SpeciesEvidence:
    """Aggregated query-segment evidence for one species."""

    total_weight: float = 0.0
    independent_weight: float = 0.0
    segment_count: int = 0
    independent_segment_count: int = 0
    read_count: int = 0
    signature_ids: Set[int] = field(default_factory=set)
    independent_signature_ids: Set[int] = field(default_factory=set)

    @property
    def independent_fraction(self) -> float:
        if self.total_weight <= 0.0:
            return 0.0
        return min(1.0, self.independent_weight / self.total_weight)


@dataclass
class PairEvidence:
    """Aggregated evidence for an unordered pair of species."""

    shared_weight: float = 0.0
    raw_weight: float = 0.0
    shared_segments: int = 0
    shared_reads: int = 0
    locus_pairs: Set[Tuple[int, int]] = field(default_factory=set)
    loci1: Set[int] = field(default_factory=set)
    loci2: Set[int] = field(default_factory=set)
    anchor1_segments: int = 0
    anchor2_segments: int = 0
    weighted_score_gap_sum: float = 0.0
    weighted_overlap_sum: float = 0.0
    mean_weight_denominator: float = 0.0


@dataclass
class EvidenceResult:
    species_evidence: List[SpeciesEvidence]
    pair_evidence: Dict[Tuple[int, int], PairEvidence]
    read_degree_histogram: Counter
    segment_degree_histogram: Counter
    stats: Counter


@dataclass
class GroupResult:
    groups: Dict[int, List[str]]
    species_to_group: Dict[str, int]
    group_meta: Dict[int, Dict[str, object]]
    shadows: Dict[str, str]
    edge_rows: List[Dict[str, object]]


class LinkageObservationStore:
    """Memory- or SQLite-backed storage for coordinate-sorted BAM input."""

    def __init__(
        self,
        mode: str = "sqlite",
        temp_dir: Optional[str] = None,
        keep: bool = False,
    ) -> None:
        if mode not in {"memory", "sqlite"}:
            raise ValueError("linkage storage must be 'memory' or 'sqlite'")

        self.mode = mode
        self.keep = keep
        self.n_observations = 0
        self.path: Optional[str] = None
        self._memory: Optional[Dict[str, List[LinkageObservation]]] = None
        self._conn: Optional[sqlite3.Connection] = None

        if mode == "memory":
            self._memory = defaultdict(list)
            return

        fd, path = tempfile.mkstemp(
            prefix="gottcha2_linkage_",
            suffix=".sqlite",
            dir=temp_dir,
        )
        os.close(fd)
        self.path = path
        self._conn = sqlite3.connect(path)
        self._conn.execute("PRAGMA journal_mode=OFF")
        self._conn.execute("PRAGMA synchronous=OFF")
        self._conn.execute("PRAGMA temp_store=FILE")
        self._conn.execute("PRAGMA cache_size=-65536")  # approximately 64 MiB
        self._conn.execute("PRAGMA locking_mode=EXCLUSIVE")
        self._conn.execute(
            """
            CREATE TABLE observations (
                read_key TEXT NOT NULL,
                qstart INTEGER NOT NULL,
                qend INTEGER NOT NULL,
                species_idx INTEGER NOT NULL,
                reference_idx INTEGER NOT NULL,
                score REAL NOT NULL,
                identity REAL NOT NULL,
                aligned_length INTEGER NOT NULL,
                alignment_class INTEGER NOT NULL,
                has_as INTEGER NOT NULL
            )
            """
        )

    def add_many(self, observations: Sequence[LinkageObservation]) -> None:
        if not observations:
            return

        self.n_observations += len(observations)
        if self.mode == "memory":
            assert self._memory is not None
            for observation in observations:
                self._memory[observation.read_key].append(observation)
            return

        assert self._conn is not None
        self._conn.executemany(
            "INSERT INTO observations VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)",
            observations,
        )

    def iter_reads(self) -> Iterator[Tuple[str, List[LinkageObservation]]]:
        """Yield all observations for one biological query at a time."""

        if self.mode == "memory":
            assert self._memory is not None
            yield from self._memory.items()
            return

        assert self._conn is not None
        self._conn.commit()
        self._conn.execute(
            "CREATE INDEX IF NOT EXISTS idx_observations_read "
            "ON observations(read_key, qstart, qend, score DESC)"
        )
        cursor = self._conn.execute(
            "SELECT read_key, qstart, qend, species_idx, reference_idx, "
            "score, identity, aligned_length, alignment_class, has_as "
            "FROM observations "
            "ORDER BY read_key, qstart, qend, score DESC"
        )
        for read_key, rows in groupby(cursor, key=lambda row: row[0]):
            yield read_key, [LinkageObservation(*row) for row in rows]

    def close(self) -> None:
        if self._conn is not None:
            self._conn.close()
            self._conn = None
        if self.path and not self.keep:
            try:
                os.remove(self.path)
            except FileNotFoundError:
                pass


def make_read_key(
    query_name: str,
    read_group: str = "",
    mate_number: int = 0,
) -> str:
    """Combine read group, QNAME, and mate number into an exact query key."""

    return f"{read_group}\x1f{query_name}\x1f{mate_number}"


def query_interval_from_cigar(
    cigartuples: Sequence[Tuple[int, int]],
    is_reverse: bool,
) -> Optional[Tuple[int, int]]:
    """Return the aligned interval in original-read coordinates.

    Hard and soft clipping are included in the reconstructed read coordinate
    system. CIGAR operation codes follow the SAM/pysam convention.
    """

    if not cigartuples:
        return None

    query_consuming = {0, 1, 4, 5, 7, 8}  # M/I/S/H/=/X
    aligned_query_consuming = {0, 1, 7, 8}  # M/I/=/X
    query_span = sum(
        length for operation, length in cigartuples
        if operation in query_consuming
    )
    aligned_query_length = sum(
        length for operation, length in cigartuples
        if operation in aligned_query_consuming
    )
    if query_span <= 0 or aligned_query_length <= 0:
        return None

    left_clip = 0
    for operation, length in cigartuples:
        if operation not in (4, 5):
            break
        left_clip += length

    right_clip = 0
    for operation, length in reversed(cigartuples):
        if operation not in (4, 5):
            break
        right_clip += length

    if is_reverse:
        qstart = right_clip
        qend = query_span - left_clip
    else:
        qstart = left_clip
        qend = query_span - right_clip

    if qstart < 0 or qend <= qstart:
        return None
    return qstart, qend


def observation_rank(observation: LinkageObservation) -> Tuple:
    """Rank observations within one homologous query interval."""

    return (
        observation.score,
        observation.identity,
        observation.aligned_length,
        -observation.alignment_class,
    )


def query_overlap_fraction(
    first: LinkageObservation,
    second: LinkageObservation,
) -> float:
    """Intersection divided by the shorter aligned query interval."""

    intersection = shared_query_length(first, second)
    shorter = min(first.qend - first.qstart, second.qend - second.qstart)
    if shorter <= 0:
        return 0.0
    return intersection / shorter


def shared_query_length(
    first: LinkageObservation,
    second: LinkageObservation,
) -> int:
    return max(
        0,
        min(first.qend, second.qend) - max(first.qstart, second.qstart),
    )


def cluster_read_observations(
    observations: Sequence[LinkageObservation],
    min_query_overlap: float,
) -> List[List[LinkageObservation]]:
    """Form non-transitive, anchor-centered homologous query clusters."""

    best_exact: Dict[Tuple[int, int, int, int], LinkageObservation] = {}
    for observation in observations:
        key = (
            observation.qstart,
            observation.qend,
            observation.species_idx,
            observation.reference_idx,
        )
        current = best_exact.get(key)
        if current is None or observation_rank(observation) > observation_rank(current):
            best_exact[key] = observation

    ordered = sorted(best_exact.values(), key=observation_rank, reverse=True)
    clusters: List[List[LinkageObservation]] = []

    for observation in ordered:
        best_cluster_idx: Optional[int] = None
        best_overlap = -1.0
        for cluster_idx, cluster in enumerate(clusters):
            overlap = query_overlap_fraction(observation, cluster[0])
            if overlap < min_query_overlap:
                continue
            if overlap > best_overlap:
                best_overlap = overlap
                best_cluster_idx = cluster_idx

        if best_cluster_idx is None:
            clusters.append([observation])
        else:
            clusters[best_cluster_idx].append(observation)

    return clusters


def _best_species_candidates(
    cluster: Sequence[LinkageObservation],
    max_species_per_segment: int,
) -> Tuple[List[LinkageObservation], int]:
    best_by_species: Dict[int, LinkageObservation] = {}
    for observation in cluster:
        current = best_by_species.get(observation.species_idx)
        if current is None or observation_rank(observation) > observation_rank(current):
            best_by_species[observation.species_idx] = observation

    candidates = sorted(
        best_by_species.values(),
        key=observation_rank,
        reverse=True,
    )
    original_degree = len(candidates)
    if max_species_per_segment > 0:
        candidates = candidates[:max_species_per_segment]
    return candidates, original_degree


def build_linkage_evidence(
    read_groups: Iterable[Tuple[str, List[LinkageObservation]]],
    n_species: int,
    config: Optional[EvidenceConfig] = None,
) -> EvidenceResult:
    """Aggregate score-weighted evidence across reads and query intervals."""

    config = config or EvidenceConfig()
    config.validate()

    species_evidence = [SpeciesEvidence() for _ in range(n_species)]
    pair_evidence: Dict[Tuple[int, int], PairEvidence] = {}
    read_degree_histogram: Counter = Counter()
    segment_degree_histogram: Counter = Counter()
    stats: Counter = Counter()

    for _read_key, observations in read_groups:
        stats["reads"] += 1
        stats["observations"] += len(observations)
        stats["fallback_score_observations"] += sum(
            1 for observation in observations if not observation.has_as
        )

        clusters = cluster_read_observations(
            observations,
            min_query_overlap=config.min_query_overlap,
        )
        read_species: Set[int] = set()
        pairs_seen_on_read: Set[Tuple[int, int]] = set()

        for cluster in clusters:
            candidates, original_degree = _best_species_candidates(
                cluster,
                max_species_per_segment=config.max_species_per_segment,
            )
            if not candidates:
                continue

            stats["query_segments"] += 1
            if original_degree > len(candidates):
                stats["truncated_query_segments"] += 1

            anchor = candidates[0]
            alternatives = []
            for alternative in candidates[1:]:
                overlap = query_overlap_fraction(anchor, alternative)
                shared_length = shared_query_length(anchor, alternative)
                if overlap < config.min_query_overlap or shared_length <= 0:
                    continue

                normalized_score_gap = max(
                    0.0,
                    (anchor.score - alternative.score) / shared_length,
                )
                score_factor = math.exp(
                    -normalized_score_gap / config.score_tau
                )
                length_factor = (
                    1.0
                    if config.length_saturation <= 0
                    else min(1.0, shared_length / config.length_saturation)
                )
                raw_weight = overlap * length_factor * score_factor
                if raw_weight < config.min_alignment_weight:
                    stats["weak_alternatives_removed"] += 1
                    continue
                alternatives.append(
                    (
                        alternative,
                        overlap,
                        normalized_score_gap,
                        raw_weight,
                    )
                )

            anchor_evidence = species_evidence[anchor.species_idx]
            anchor_evidence.total_weight += 1.0
            anchor_evidence.segment_count += 1
            anchor_evidence.signature_ids.add(anchor.reference_idx)
            read_species.add(anchor.species_idx)

            if not alternatives:
                anchor_evidence.independent_weight += 1.0
                anchor_evidence.independent_segment_count += 1
                anchor_evidence.independent_signature_ids.add(anchor.reference_idx)
                segment_degree_histogram[1] += 1
                continue

            raw_weight_sum = sum(item[3] for item in alternatives)
            scale = max(1.0, raw_weight_sum)
            segment_degree_histogram[1 + len(alternatives)] += 1

            best_alternative_gap = min(item[2] for item in alternatives)
            if best_alternative_gap >= config.independent_score_margin:
                anchor_evidence.independent_weight += 1.0
                anchor_evidence.independent_segment_count += 1
                anchor_evidence.independent_signature_ids.add(anchor.reference_idx)

            for alternative, overlap, normalized_score_gap, raw_weight in alternatives:
                capped_weight = raw_weight / scale
                alternative_evidence = species_evidence[alternative.species_idx]
                alternative_evidence.total_weight += capped_weight
                alternative_evidence.segment_count += 1
                alternative_evidence.signature_ids.add(alternative.reference_idx)
                read_species.add(alternative.species_idx)

                if anchor.species_idx < alternative.species_idx:
                    pair_key = (anchor.species_idx, alternative.species_idx)
                    locus_pair = (anchor.reference_idx, alternative.reference_idx)
                    anchor_is_species1 = True
                else:
                    pair_key = (alternative.species_idx, anchor.species_idx)
                    locus_pair = (alternative.reference_idx, anchor.reference_idx)
                    anchor_is_species1 = False

                pair = pair_evidence.setdefault(pair_key, PairEvidence())
                pair.shared_weight += capped_weight
                pair.raw_weight += raw_weight
                pair.shared_segments += 1
                pair.locus_pairs.add(locus_pair)
                pair.loci1.add(locus_pair[0])
                pair.loci2.add(locus_pair[1])
                pair.weighted_score_gap_sum += (
                    normalized_score_gap * capped_weight
                )
                pair.weighted_overlap_sum += overlap * capped_weight
                pair.mean_weight_denominator += capped_weight
                if anchor_is_species1:
                    pair.anchor1_segments += 1
                else:
                    pair.anchor2_segments += 1
                pairs_seen_on_read.add(pair_key)

        for species_idx in read_species:
            species_evidence[species_idx].read_count += 1
        for pair_key in pairs_seen_on_read:
            pair_evidence[pair_key].shared_reads += 1

        read_degree_histogram[len(read_species)] += 1
        if len(read_species) > 1:
            stats["multi_species_reads"] += 1

    return EvidenceResult(
        species_evidence=species_evidence,
        pair_evidence=pair_evidence,
        read_degree_histogram=read_degree_histogram,
        segment_degree_histogram=segment_degree_histogram,
        stats=stats,
    )


def _taxid_sort_key(taxid: str) -> Tuple[int, object]:
    try:
        return 0, int(taxid)
    except (TypeError, ValueError):
        return 1, str(taxid)


def finalize_species_groups(
    idx_to_species: List[str],
    idx_to_genus: List[Optional[str]],
    evidence_result: EvidenceResult,
    config: Optional[GroupConfig] = None,
) -> GroupResult:
    """Build direct anchor-centered groups and conservative shadow labels."""

    config = config or GroupConfig()
    config.validate()
    species_evidence = evidence_result.species_evidence
    pair_evidence = evidence_result.pair_evidence
    decisions: Dict[Tuple[int, int], Dict[str, object]] = {}

    for (idx1, idx2), pair in pair_evidence.items():
        evidence1 = species_evidence[idx1]
        evidence2 = species_evidence[idx2]
        containment1 = (
            min(1.0, pair.shared_weight / evidence1.total_weight)
            if evidence1.total_weight > 0.0
            else 0.0
        )
        containment2 = (
            min(1.0, pair.shared_weight / evidence2.total_weight)
            if evidence2.total_weight > 0.0
            else 0.0
        )
        mutual_ambiguity = min(containment1, containment2)
        genus1 = idx_to_genus[idx1]
        genus2 = idx_to_genus[idx2]
        same_genus = (
            genus1 is not None and genus2 is not None and genus1 == genus2
        )

        pass_reads = pair.shared_reads >= config.min_link_reads
        pass_segments = pair.shared_segments >= config.min_link_segments
        # Require independent signature support on both sides. Counting only
        # signature pairs can be inflated when one conserved signature maps to
        # many fragments in the other species.
        loci1 = pair.loci1 or {first for first, _ in pair.locus_pairs}
        loci2 = pair.loci2 or {second for _, second in pair.locus_pairs}
        shared_loci = min(len(loci1), len(loci2))
        pass_loci = shared_loci >= config.min_link_loci
        pass_weight = pair.shared_weight >= config.min_link_weight
        pass_genus = (not config.same_genus_only) or same_genus
        pass_containment = (
            max(containment1, containment2)
            >= config.min_shadow_containment
            or mutual_ambiguity >= config.min_mutual_ambiguity
        )
        base_qualified = (
            pass_reads
            and pass_segments
            and pass_loci
            and pass_weight
            and pass_genus
        )
        decisions[(idx1, idx2)] = {
            "containment1": containment1,
            "containment2": containment2,
            "mutual_ambiguity": mutual_ambiguity,
            "same_genus": same_genus,
            "pass_reads": pass_reads,
            "pass_segments": pass_segments,
            "shared_loci": shared_loci,
            "pass_loci": pass_loci,
            "pass_weight": pass_weight,
            "pass_genus": pass_genus,
            "pass_containment": pass_containment,
            "base_qualified": base_qualified,
            "qualified": base_qualified and pass_containment,
        }

    ranked_species = sorted(
        range(len(idx_to_species)),
        key=lambda idx: (
            -species_evidence[idx].total_weight,
            -species_evidence[idx].independent_weight,
            -species_evidence[idx].read_count,
            -len(species_evidence[idx].independent_signature_ids),
            _taxid_sort_key(idx_to_species[idx]),
        ),
    )

    assigned: Set[int] = set()
    used_pairs: Set[Tuple[int, int]] = set()
    groups: Dict[int, List[str]] = {}
    species_to_group: Dict[str, int] = {}
    group_meta: Dict[int, Dict[str, object]] = {}
    shadows: Dict[str, str] = {}

    for anchor_idx in ranked_species:
        if anchor_idx in assigned:
            continue
        anchor_evidence = species_evidence[anchor_idx]
        if anchor_evidence.total_weight <= 0.0:
            continue

        candidates = []
        for candidate_idx in ranked_species:
            if candidate_idx == anchor_idx or candidate_idx in assigned:
                continue
            pair_key = (
                (anchor_idx, candidate_idx)
                if anchor_idx < candidate_idx
                else (candidate_idx, anchor_idx)
            )
            decision = decisions.get(pair_key)
            if not decision or not decision["base_qualified"]:
                continue

            if anchor_idx == pair_key[0]:
                anchor_containment = float(decision["containment1"])
                candidate_containment = float(decision["containment2"])
            else:
                anchor_containment = float(decision["containment2"])
                candidate_containment = float(decision["containment1"])

            mutual_ambiguity = float(decision["mutual_ambiguity"])
            if not (
                candidate_containment >= config.min_shadow_containment
                or mutual_ambiguity >= config.min_mutual_ambiguity
            ):
                continue

            candidates.append(
                (
                    candidate_idx,
                    candidate_containment,
                    anchor_containment,
                    mutual_ambiguity,
                    pair_evidence[pair_key].shared_weight,
                    pair_key,
                )
            )

        if not candidates:
            continue

        candidates.sort(
            key=lambda item: (
                item[1],
                item[3],
                item[4],
                species_evidence[item[0]].total_weight,
            ),
            reverse=True,
        )

        member_indices = [anchor_idx]
        shadow_taxids: List[str] = []
        ambiguous_taxids: List[str] = []
        for candidate_idx, candidate_containment, _, _, _, pair_key in candidates:
            if candidate_idx in assigned:
                continue

            candidate_evidence = species_evidence[candidate_idx]
            is_shadow = (
                candidate_containment >= config.min_shadow_containment
                and candidate_evidence.independent_fraction
                <= config.max_shadow_independent_fraction
                and len(candidate_evidence.independent_signature_ids)
                < config.min_independent_loci_to_retain
                and len(anchor_evidence.independent_signature_ids)
                >= config.min_anchor_independent_loci_for_shadow
                and anchor_evidence.total_weight
                >= candidate_evidence.total_weight
                and anchor_evidence.independent_fraction
                >= candidate_evidence.independent_fraction
            )

            member_indices.append(candidate_idx)
            used_pairs.add(pair_key)
            candidate_taxid = idx_to_species[candidate_idx]
            if is_shadow:
                anchor_taxid = idx_to_species[anchor_idx]
                shadows[candidate_taxid] = anchor_taxid
                shadow_taxids.append(candidate_taxid)
            else:
                ambiguous_taxids.append(candidate_taxid)

        if len(member_indices) < 2:
            continue

        group_id = len(groups) + 1
        assigned.update(member_indices)
        member_taxids = sorted(
            (idx_to_species[idx] for idx in member_indices),
            key=_taxid_sort_key,
        )
        anchor_taxid = idx_to_species[anchor_idx]
        group_type = (
            "anchor_with_shadows"
            if shadow_taxids and not ambiguous_taxids
            else "ambiguous_or_mixed"
        )
        groups[group_id] = member_taxids
        group_meta[group_id] = {
            "anchor_taxid": anchor_taxid,
            "group_type": group_type,
            "shadow_taxids": sorted(shadow_taxids, key=_taxid_sort_key),
            "ambiguous_taxids": sorted(ambiguous_taxids, key=_taxid_sort_key),
        }
        for taxid in member_taxids:
            species_to_group[taxid] = group_id

    edge_rows: List[Dict[str, object]] = []
    for (idx1, idx2), pair in pair_evidence.items():
        decision = decisions[(idx1, idx2)]
        species1 = idx_to_species[idx1]
        species2 = idx_to_species[idx2]
        group1 = species_to_group.get(species1)
        group2 = species_to_group.get(species2)
        denominator = pair.mean_weight_denominator

        shadow_direction = ""
        if shadows.get(species1) == species2:
            shadow_direction = f"{species1}->{species2}"
        elif shadows.get(species2) == species1:
            shadow_direction = f"{species2}->{species1}"

        edge_rows.append(
            {
                "SPECIES1": species1,
                "SPECIES2": species2,
                "GENUS1": idx_to_genus[idx1] or "",
                "GENUS2": idx_to_genus[idx2] or "",
                "SAME_GENUS": int(bool(decision["same_genus"])),
                "SHARED_READS": pair.shared_reads,
                "SHARED_QUERY_SEGMENTS": pair.shared_segments,
                "SHARED_SIGNATURES1": len(
                    pair.loci1 or {first for first, _ in pair.locus_pairs}
                ),
                "SHARED_SIGNATURES2": len(
                    pair.loci2 or {second for _, second in pair.locus_pairs}
                ),
                "SHARED_SIGNATURE_PAIRS": len(pair.locus_pairs),
                "SHARED_INDEPENDENT_LOCI": decision["shared_loci"],
                "WEIGHTED_SUPPORT": pair.shared_weight,
                "RAW_WEIGHTED_SUPPORT": pair.raw_weight,
                "TOTAL_WEIGHT1": species_evidence[idx1].total_weight,
                "TOTAL_WEIGHT2": species_evidence[idx2].total_weight,
                "CONTAINMENT1_BY_2": decision["containment1"],
                "CONTAINMENT2_BY_1": decision["containment2"],
                "MUTUAL_AMBIGUITY": decision["mutual_ambiguity"],
                "ANCHOR_SEGMENTS1": pair.anchor1_segments,
                "ANCHOR_SEGMENTS2": pair.anchor2_segments,
                "MEAN_NORMALIZED_SCORE_GAP": (
                    pair.weighted_score_gap_sum / denominator
                    if denominator > 0.0
                    else 0.0
                ),
                "MEAN_QUERY_OVERLAP": (
                    pair.weighted_overlap_sum / denominator
                    if denominator > 0.0
                    else 0.0
                ),
                "PASS_MIN_READS": int(bool(decision["pass_reads"])),
                "PASS_MIN_SEGMENTS": int(bool(decision["pass_segments"])),
                "PASS_MIN_SIGNATURES": int(bool(decision["pass_loci"])),
                "PASS_MIN_WEIGHT": int(bool(decision["pass_weight"])),
                "PASS_GENUS": int(bool(decision["pass_genus"])),
                "PASS_CONTAINMENT": int(bool(decision["pass_containment"])),
                "QUALIFIED": int(bool(decision["qualified"])),
                "USED_FOR_GROUPING": int((idx1, idx2) in used_pairs),
                "SHADOW_DIRECTION": shadow_direction,
                "GROUP1": group1 if group1 is not None else "",
                "GROUP2": group2 if group2 is not None else "",
            }
        )

    edge_rows.sort(
        key=lambda row: (
            -int(row["USED_FOR_GROUPING"]),
            -int(row["QUALIFIED"]),
            -float(row["WEIGHTED_SUPPORT"]),
            -int(row["SHARED_READS"]),
            _taxid_sort_key(str(row["SPECIES1"])),
            _taxid_sort_key(str(row["SPECIES2"])),
        )
    )
    return GroupResult(
        groups=groups,
        species_to_group=species_to_group,
        group_meta=group_meta,
        shadows=shadows,
        edge_rows=edge_rows,
    )


def write_species_groups(result: GroupResult, out_path: str) -> None:
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(
            "GROUP_ID\tN_SPECIES\tSPECIES_TAXIDS\tANCHOR_TAXID\t"
            "GROUP_TYPE\tSHADOW_TAXIDS\tAMBIGUOUS_TAXIDS\n"
        )
        for group_id, species in result.groups.items():
            meta = result.group_meta[group_id]
            out.write(
                f"{group_id}\t{len(species)}\t{','.join(species)}\t"
                f"{meta['anchor_taxid']}\t{meta['group_type']}\t"
                f"{','.join(meta['shadow_taxids'])}\t"
                f"{','.join(meta['ambiguous_taxids'])}\n"
            )


def write_species_links(edge_rows: List[Dict[str, object]], out_path: str) -> None:
    fields = [
        "SPECIES1", "SPECIES2", "GENUS1", "GENUS2", "SAME_GENUS",
        "SHARED_READS", "SHARED_QUERY_SEGMENTS", "SHARED_SIGNATURES1",
        "SHARED_SIGNATURES2", "SHARED_SIGNATURE_PAIRS",
        "SHARED_INDEPENDENT_LOCI",
        "WEIGHTED_SUPPORT", "RAW_WEIGHTED_SUPPORT", "TOTAL_WEIGHT1",
        "TOTAL_WEIGHT2", "CONTAINMENT1_BY_2", "CONTAINMENT2_BY_1",
        "MUTUAL_AMBIGUITY", "ANCHOR_SEGMENTS1", "ANCHOR_SEGMENTS2",
        "MEAN_NORMALIZED_SCORE_GAP", "MEAN_QUERY_OVERLAP",
        "PASS_MIN_READS", "PASS_MIN_SEGMENTS", "PASS_MIN_SIGNATURES",
        "PASS_MIN_WEIGHT", "PASS_GENUS", "PASS_CONTAINMENT", "QUALIFIED",
        "USED_FOR_GROUPING", "SHADOW_DIRECTION", "GROUP1", "GROUP2",
    ]
    float_fields = {
        "WEIGHTED_SUPPORT", "RAW_WEIGHTED_SUPPORT", "TOTAL_WEIGHT1",
        "TOTAL_WEIGHT2", "CONTAINMENT1_BY_2", "CONTAINMENT2_BY_1",
        "MUTUAL_AMBIGUITY", "MEAN_NORMALIZED_SCORE_GAP",
        "MEAN_QUERY_OVERLAP",
    }
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("\t".join(fields) + "\n")
        for row in edge_rows:
            values = [
                f"{float(row[field_name]):.8g}"
                if field_name in float_fields
                else str(row[field_name])
                for field_name in fields
            ]
            out.write("\t".join(values) + "\n")


def write_species_evidence(
    idx_to_species: List[str],
    idx_to_genus: List[Optional[str]],
    evidence: EvidenceResult,
    groups: GroupResult,
    out_path: str,
) -> None:
    header = [
        "SPECIES_TAXID", "GENUS_TAXID", "DISTINCT_READS", "QUERY_SEGMENTS",
        "TOTAL_WEIGHT", "INDEPENDENT_WEIGHT", "INDEPENDENT_FRACTION",
        "SIGNATURES", "INDEPENDENT_SIGNATURES", "GROUP_ID",
        "IS_GROUP_ANCHOR", "SHADOW_OF", "LINKAGE_INTERPRETATION",
    ]
    rows = []
    for idx, taxid in enumerate(idx_to_species):
        value = evidence.species_evidence[idx]
        group_id = groups.species_to_group.get(taxid)
        anchor_taxid = (
            str(groups.group_meta[group_id]["anchor_taxid"])
            if group_id is not None
            else ""
        )
        shadow_of = groups.shadows.get(taxid, "")
        if shadow_of:
            interpretation = "likely_shadow;combine_with_taxon_level_evidence"
        elif group_id is not None and taxid == anchor_taxid:
            interpretation = "group_anchor"
        elif group_id is not None:
            interpretation = "ambiguous_group_member"
        else:
            interpretation = "ungrouped"
        rows.append(
            (
                taxid,
                idx_to_genus[idx] or "",
                value.read_count,
                value.segment_count,
                value.total_weight,
                value.independent_weight,
                value.independent_fraction,
                len(value.signature_ids),
                len(value.independent_signature_ids),
                group_id if group_id is not None else "",
                int(group_id is not None and taxid == anchor_taxid),
                shadow_of,
                interpretation,
            )
        )

    rows.sort(key=lambda row: (-float(row[4]), _taxid_sort_key(str(row[0]))))
    with open(out_path, "w", encoding="utf-8") as out:
        out.write("\t".join(header) + "\n")
        for row in rows:
            formatted = list(row)
            for idx in (4, 5, 6):
                formatted[idx] = f"{float(formatted[idx]):.8g}"
            out.write("\t".join(map(str, formatted)) + "\n")


def write_degree_histogram(
    histogram: Counter,
    out_path: str,
    degree_column: str,
) -> None:
    with open(out_path, "w", encoding="utf-8") as out:
        out.write(f"{degree_column}\tN_OBSERVATIONS\n")
        for degree in sorted(histogram):
            out.write(f"{degree}\t{histogram[degree]}\n")
