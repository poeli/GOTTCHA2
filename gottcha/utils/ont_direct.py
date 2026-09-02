"""Direct Oxford Nanopore read-level alignment resolution for GOTTCHA2.

The direct-ONT resolver intentionally does *not* duplicate GOTTCHA2's normal
alignment-quality filters.  It uses only four SAM columns with pandas:
QNAME, FLAG, RNAME and CIGAR.  ``process_bam.py`` remains responsible for
matchIdentity, matchFraction and matchLength.

Read consistency is determined in original ONT-read coordinates:

1. merge/union aligned intervals independently for each taxid;
2. map those taxids to species and union the already nonredundant taxon
   intervals again at species level;
3. choose a species winner from union-bp support;
4. retain all structural SAM records belonging to the winning species.

The output SAM is streamed from the input and adds:
  ZQ:i  original-read query start (0-based)
  ZE:i  original-read query end (exclusive)
  ZL:i  reconstructed original-read length, including hard clipping
  ZT:Z  alignment taxid
  ZS:Z  winning species taxid
  ZC:i  present on exactly one retained alignment per original ONT read

Hard-clipped supplementary CIGARs are supported, so minimap2 does not need -Y.
"""
from __future__ import annotations

from dataclasses import asdict, dataclass
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Sequence, Tuple
import csv
import logging
import math
import os
import re

import numpy as np
import pandas as pd

try:
    import taxonomy
except ImportError:  # package import
    import gottcha.utils.taxonomy as taxonomy


_CIGAR_RE = re.compile(r"(\d+)([MIDNSHP=X])")
_QUERY_CONSUME_ALIGNED = frozenset("MI=X")
# H does not consume query according to SAM, but it represents original-read
# sequence omitted from SEQ.  Include it when reconstructing molecule length.
_QUERY_CONSUME_ORIGINAL = frozenset("MI=XSH")
_DIRECT_TAG_PREFIXES = ("ZQ:i:", "ZE:i:", "ZL:i:", "ZT:Z:", "ZS:Z:", "ZC:i:")


@dataclass(frozen=True)
class OntDirectConfig:
    """Read-level species-consistency parameters.

    These are ambiguity-resolution thresholds, not alignment-quality
    thresholds. MI/MF/MG stay downstream in process_bam.py.
    """

    min_species_support: float = 0.65
    species_support_ratio: float = 1.50
    conflict_policy: str = "drop"  # drop | best
    max_candidates_per_read: int = 100

    def validate(self) -> None:
        if not 0 <= self.min_species_support <= 1:
            raise ValueError("min_species_support must be between 0 and 1")
        if self.species_support_ratio < 1:
            raise ValueError("species_support_ratio must be >= 1")
        if self.conflict_policy not in {"drop", "best"}:
            raise ValueError("conflict_policy must be 'drop' or 'best'")
        if self.max_candidates_per_read < 1:
            raise ValueError("max_candidates_per_read must be >= 1")


@dataclass
class OntDirectStats:
    total_alignment_records: int = 0
    structural_records: int = 0
    invalid_structural_records: int = 0
    reads_with_candidates: int = 0
    single_species_reads: int = 0
    multi_species_reads: int = 0
    ambiguous_reads: int = 0
    forced_best_reads: int = 0
    capped_candidates_removed: int = 0
    kept_reads: int = 0
    kept_alignments: int = 0

    def to_dict(self) -> Dict[str, int]:
        return asdict(self)


@dataclass(frozen=True)
class _Hit:
    idx: int
    taxid: str
    species_taxid: str
    qstart: int
    qend: int
    qlen: int
    flag: int

    @property
    def qspan(self) -> int:
        return max(0, self.qend - self.qstart)

    @property
    def is_primary(self) -> bool:
        return (self.flag & (0x100 | 0x800)) == 0


AMBIGUOUS_COLUMNS = [
    "QNAME",
    "ACTION",
    "REASON",
    "READ_LENGTH",
    "INPUT_ALIGNMENT_COUNT",
    "CONSIDERED_ALIGNMENT_COUNT",
    "CANDIDATE_CAP_REMOVED",
    "TAXON_COUNT",
    "SPECIES_COUNT",
    "TOP_SPECIES",
    "TOP_SPECIES_NAME",
    "TOP_SPECIES_BP",
    "RUNNER_UP_SPECIES",
    "RUNNER_UP_SPECIES_NAME",
    "RUNNER_UP_BP",
    "TOP_FRACTION",
    "TOP_TO_RUNNER_RATIO",
    "SPECIES_SUPPORT_BP",
    "TAXON_SUPPORT_BP",
]


@lru_cache(maxsize=250_000)
def _parse_cigar_layout(cigar: str) -> Optional[Tuple[int, int, int]]:
    """Return query start/end in alignment orientation plus original length."""
    if not cigar or cigar == "*":
        return None
    ops = [(int(n), op) for n, op in _CIGAR_RE.findall(cigar)]
    if not ops or "".join(f"{n}{op}" for n, op in ops) != cigar:
        return None

    q_aligned = sum(n for n, op in ops if op in _QUERY_CONSUME_ALIGNED)
    qlen = sum(n for n, op in ops if op in _QUERY_CONSUME_ORIGINAL)
    if q_aligned <= 0 or qlen <= 0:
        return None

    leading_clip = 0
    for n, op in ops:
        if op in {"S", "H"}:
            leading_clip += n
        else:
            break

    qstart = leading_clip
    qend = qstart + q_aligned
    if qstart < 0 or qend <= qstart or qend > qlen:
        return None
    return qstart, qend, qlen


def _taxid_from_ref(rname: str) -> Optional[str]:
    # Existing GOTTCHA2 signature names place taxid in the penultimate field.
    parts = str(rname).rsplit("|", 2)
    if len(parts) != 3 or not parts[-2]:
        return None
    return str(parts[-2])


@lru_cache(maxsize=250_000)
def _default_species_taxid(taxid: str) -> str:
    """Return nearest species taxid; fall back to the input taxid if absent."""
    try:
        species = taxonomy.taxid2taxidOnRank(str(taxid), "species")
    except Exception:
        species = ""
    if species in (None, "", "unknown", 0, "0"):
        return str(taxid)
    return str(species)


@lru_cache(maxsize=250_000)
def _tax_name(taxid: str) -> str:
    try:
        name = taxonomy.taxid2name(str(taxid))
    except Exception:
        return ""
    return "" if name is None else str(name)


def _union_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    clean = sorted((int(s), int(e)) for s, e in intervals if int(e) > int(s))
    if not clean:
        return []
    merged: List[Tuple[int, int]] = []
    cur_s, cur_e = clean[0]
    for s, e in clean[1:]:
        if s <= cur_e:
            cur_e = max(cur_e, e)
        else:
            merged.append((cur_s, cur_e))
            cur_s, cur_e = s, e
    merged.append((cur_s, cur_e))
    return merged


def _interval_length(intervals: Iterable[Tuple[int, int]]) -> int:
    return sum(e - s for s, e in intervals)


def _taxon_union_intervals(hits: Sequence[_Hit]) -> Dict[str, List[Tuple[int, int]]]:
    by_taxon: Dict[str, List[Tuple[int, int]]] = {}
    for hit in hits:
        by_taxon.setdefault(hit.taxid, []).append((hit.qstart, hit.qend))
    return {taxid: _union_intervals(intervals) for taxid, intervals in by_taxon.items()}


def _species_union_intervals(
    taxon_intervals: Dict[str, List[Tuple[int, int]]],
    taxon_to_species: Dict[str, str],
) -> Dict[str, List[Tuple[int, int]]]:
    by_species: Dict[str, List[Tuple[int, int]]] = {}
    for taxid, intervals in taxon_intervals.items():
        species = taxon_to_species[taxid]
        # These intervals are already nonredundant for the taxon. Union again at
        # species level to avoid double counting the same ONT-read bases across
        # multiple strains/taxids within one species.
        by_species.setdefault(species, []).extend(intervals)
    return {species: _union_intervals(intervals) for species, intervals in by_species.items()}


def _support_lengths(intervals_by_key: Dict[str, List[Tuple[int, int]]]) -> Dict[str, int]:
    return {key: _interval_length(intervals) for key, intervals in intervals_by_key.items()}


def _rank_hit(hit: _Hit) -> Tuple[int, int, int]:
    # No MAPQ/AS/NM is used. Prefer wider structural evidence, then primary,
    # then earlier SAM order for deterministic behavior.
    return (hit.qspan, 1 if hit.is_primary else 0, -hit.idx)


def _choose_winner(
    species_support: Dict[str, int], cfg: OntDirectConfig
) -> Tuple[Optional[str], Optional[str], float, float, bool]:
    """Return winner, runner-up, winner fraction, ratio, and clear flag."""
    if not species_support:
        return None, None, 0.0, 0.0, False
    ranked = sorted(species_support.items(), key=lambda kv: (-kv[1], kv[0]))
    winner, winner_bp = ranked[0]
    if len(ranked) == 1:
        return winner, None, 1.0, math.inf, True
    runner, runner_bp = ranked[1]
    total = sum(bp for _, bp in ranked)
    fraction = winner_bp / total if total else 0.0
    ratio = math.inf if runner_bp <= 0 else winner_bp / runner_bp
    clear = (
        fraction >= cfg.min_species_support
        and ratio >= cfg.species_support_ratio
    )
    return winner, runner, fraction, ratio, clear


def _cap_hits(hits: Sequence[_Hit], limit: int) -> Tuple[List[_Hit], int]:
    if len(hits) <= limit:
        return list(hits), 0
    ranked = sorted(hits, key=_rank_hit, reverse=True)
    kept = ranked[:limit]
    kept.sort(key=lambda h: h.idx)
    return kept, len(hits) - limit


def _resolve_read_hits(
    hits: Sequence[_Hit], cfg: OntDirectConfig
) -> Tuple[List[_Hit], Dict[str, object]]:
    """Resolve one ONT read using union-bp support, taxon then species."""
    input_count = len(hits)
    working, capped = _cap_hits(hits, cfg.max_candidates_per_read)
    qlen = max((h.qlen for h in working), default=0)

    taxon_intervals = _taxon_union_intervals(working)
    taxon_to_species: Dict[str, str] = {}
    for h in working:
        taxon_to_species[h.taxid] = h.species_taxid
    species_intervals = _species_union_intervals(taxon_intervals, taxon_to_species)
    taxon_support = _support_lengths(taxon_intervals)
    species_support = _support_lengths(species_intervals)

    winner, runner, fraction, ratio, clear = _choose_winner(species_support, cfg)
    ambiguous = bool(winner is not None and len(species_support) > 1 and not clear)

    if winner is None:
        kept: List[_Hit] = []
        action = "drop"
        reason = "no_species_support"
    elif ambiguous and cfg.conflict_policy == "drop":
        kept = []
        action = "drop"
        reason = "ambiguous_species"
    else:
        kept = [h for h in working if h.species_taxid == winner]
        kept.sort(key=lambda h: h.idx)
        if ambiguous:
            action = "keep_best"
            reason = "forced_best"
        else:
            action = "keep"
            reason = "clear_species"

    diag: Dict[str, object] = {
        "input_count": input_count,
        "considered_count": len(working),
        "capped": capped,
        "qlen": qlen,
        "taxon_support": taxon_support,
        "species_support": species_support,
        "winner": winner,
        "runner": runner,
        "winner_fraction": fraction,
        "ratio": ratio,
        "ambiguous": ambiguous,
        "action": action,
        "reason": reason,
    }
    return kept, diag


def _load_structural_sam(
    samfile: str | os.PathLike[str],
    species_lookup: Callable[[str], str] = _default_species_taxid,
) -> Tuple[pd.DataFrame, int]:
    """Load only QNAME, FLAG, RNAME and CIGAR with pandas."""
    if Path(samfile).stat().st_size == 0:
        return pd.DataFrame(columns=["QNAME", "FLAG", "REF", "CIGAR"]), 0
    try:
        df = pd.read_csv(
            samfile,
            sep="\t",
            header=None,
            usecols=[0, 1, 2, 5],
            names=["QNAME", "FLAG", "REF", "CIGAR"],
            dtype={
                "QNAME": "str",
                "FLAG": "uint16",
                "REF": "str",
                "CIGAR": "str",
            },
            low_memory=False,
            memory_map=True,
        )
    except pd.errors.EmptyDataError:
        df = pd.DataFrame(columns=["QNAME", "FLAG", "REF", "CIGAR"])
    total = len(df)
    if total == 0:
        return df, total

    bad_flag = (df["FLAG"].astype("uint16") & (0x4 | 0x200 | 0x400)).astype(bool)
    valid = (~bad_flag) & df["REF"].ne("*") & df["CIGAR"].ne("*")
    df = df.loc[valid].copy()
    if df.empty:
        return df, total

    df["TAXID"] = df["REF"].map(_taxid_from_ref)
    df = df.loc[df["TAXID"].notna()].copy()
    if df.empty:
        return df, total

    layouts = df["CIGAR"].map(_parse_cigar_layout)
    ok = layouts.notna()
    df = df.loc[ok].copy()
    layouts = layouts.loc[ok]
    if df.empty:
        return df, total

    parsed = np.asarray(layouts.tolist(), dtype=np.int64)
    qstart = parsed[:, 0]
    qend = parsed[:, 1]
    qlen = parsed[:, 2]
    reverse = (df["FLAG"].to_numpy(dtype=np.uint16) & 0x10) != 0
    orig_start = np.where(reverse, qlen - qend, qstart)
    orig_end = np.where(reverse, qlen - qstart, qend)

    df["QSTART"] = orig_start.astype(np.int32, copy=False)
    df["QEND"] = orig_end.astype(np.int32, copy=False)
    df["QLEN"] = qlen.astype(np.int32, copy=False)
    df["QSPAN"] = (orig_end - orig_start).astype(np.int32, copy=False)
    df = df.loc[
        (df["QSPAN"] > 0)
        & (df["QSTART"] >= 0)
        & (df["QEND"] <= df["QLEN"])
    ].copy()
    if df.empty:
        return df, total

    unique_taxids = pd.unique(df["TAXID"])
    species_map = {str(t): str(species_lookup(str(t))) for t in unique_taxids}
    df["SPECIES"] = df["TAXID"].map(species_map)
    return df, total


def _group_to_hits(group: pd.DataFrame) -> List[_Hit]:
    return [
        _Hit(
            idx=int(row.Index),
            taxid=str(row.TAXID),
            species_taxid=str(row.SPECIES),
            qstart=int(row.QSTART),
            qend=int(row.QEND),
            qlen=int(row.QLEN),
            flag=int(row.FLAG),
        )
        for row in group.itertuples(index=True)
    ]


def _remove_existing_direct_tags(fields: List[str]) -> List[str]:
    return [f for f in fields if not f.startswith(_DIRECT_TAG_PREFIXES)]


def _format_support(support: Dict[str, int]) -> str:
    ranked = sorted(support.items(), key=lambda kv: (-kv[1], kv[0]))
    return ";".join(f"{taxid}:{bp}" for taxid, bp in ranked)


def _ratio_text(value: float) -> str:
    return "inf" if math.isinf(value) else f"{value:.6g}"


def _ambiguous_row(qname: str, diag: Dict[str, object]) -> Dict[str, object]:
    species_support = diag["species_support"]
    taxon_support = diag["taxon_support"]
    assert isinstance(species_support, dict)
    assert isinstance(taxon_support, dict)
    winner = diag["winner"]
    runner = diag["runner"]
    return {
        "QNAME": qname,
        "ACTION": diag["action"],
        "REASON": diag["reason"],
        "READ_LENGTH": diag["qlen"],
        "INPUT_ALIGNMENT_COUNT": diag["input_count"],
        "CONSIDERED_ALIGNMENT_COUNT": diag["considered_count"],
        "CANDIDATE_CAP_REMOVED": diag["capped"],
        "TAXON_COUNT": len(taxon_support),
        "SPECIES_COUNT": len(species_support),
        "TOP_SPECIES": winner or "",
        "TOP_SPECIES_NAME": _tax_name(str(winner)) if winner else "",
        "TOP_SPECIES_BP": species_support.get(winner, 0) if winner else 0,
        "RUNNER_UP_SPECIES": runner or "",
        "RUNNER_UP_SPECIES_NAME": _tax_name(str(runner)) if runner else "",
        "RUNNER_UP_BP": species_support.get(runner, 0) if runner else 0,
        "TOP_FRACTION": f"{float(diag['winner_fraction']):.6g}",
        "TOP_TO_RUNNER_RATIO": _ratio_text(float(diag["ratio"])),
        "SPECIES_SUPPORT_BP": _format_support(species_support),
        "TAXON_SUPPORT_BP": _format_support(taxon_support),
    }


def resolve_sam(
    input_sam: str | os.PathLike[str],
    output_sam: str | os.PathLike[str],
    config: Optional[OntDirectConfig] = None,
    ambiguous_tsv: Optional[str | os.PathLike[str]] = None,
    species_lookup: Callable[[str], str] = _default_species_taxid,
) -> OntDirectStats:
    """Resolve intact ONT reads to one winning species and write curated SAM."""
    cfg = config or OntDirectConfig()
    cfg.validate()
    stats = OntDirectStats()

    input_path = Path(input_sam)
    output_path = Path(output_sam)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    logging.info("Loading direct-ONT structural SAM fields with pandas...")
    df, total = _load_structural_sam(input_path, species_lookup=species_lookup)
    stats.total_alignment_records = total
    stats.structural_records = len(df)
    stats.invalid_structural_records = total - len(df)

    keep = np.zeros(total, dtype=np.bool_)
    count_once = np.zeros(total, dtype=np.bool_)
    qstart_all = np.full(total, -1, dtype=np.int32)
    qend_all = np.full(total, -1, dtype=np.int32)
    qlen_all = np.full(total, -1, dtype=np.int32)
    species_all = np.empty(total, dtype=object)
    species_all[:] = ""

    tsv_handle = None
    tsv_writer = None
    if ambiguous_tsv:
        ambiguous_path = Path(ambiguous_tsv)
        ambiguous_path.parent.mkdir(parents=True, exist_ok=True)
        tsv_handle = ambiguous_path.open("w", encoding="utf-8", newline="")
        tsv_writer = csv.DictWriter(tsv_handle, fieldnames=AMBIGUOUS_COLUMNS, delimiter="\t")
        tsv_writer.writeheader()

    try:
        if not df.empty:
            original_idx = df.index.to_numpy(dtype=np.int64)
            qstart_all[original_idx] = df["QSTART"].to_numpy(dtype=np.int32)
            qend_all[original_idx] = df["QEND"].to_numpy(dtype=np.int32)
            qlen_all[original_idx] = df["QLEN"].to_numpy(dtype=np.int32)
            species_all[original_idx] = df["SPECIES"].to_numpy(dtype=object)

            stats.reads_with_candidates = int(df["QNAME"].nunique())
            for qname, group in df.groupby("QNAME", sort=False):
                hits = _group_to_hits(group)
                kept_hits, diag = _resolve_read_hits(hits, cfg)
                stats.capped_candidates_removed += int(diag["capped"])

                species_count = len(diag["species_support"])
                if species_count <= 1:
                    stats.single_species_reads += 1
                else:
                    stats.multi_species_reads += 1

                if bool(diag["ambiguous"]):
                    stats.ambiguous_reads += 1
                    if diag["action"] == "keep_best":
                        stats.forced_best_reads += 1
                    if tsv_writer is not None:
                        tsv_writer.writerow(_ambiguous_row(str(qname), diag))

                if not kept_hits:
                    continue

                kept_idx = np.fromiter((h.idx for h in kept_hits), dtype=np.int64)
                keep[kept_idx] = True
                anchor = max(kept_hits, key=_rank_hit)
                count_once[anchor.idx] = True
                stats.kept_reads += 1
                stats.kept_alignments += len(kept_hits)
    finally:
        if tsv_handle is not None:
            tsv_handle.close()

    logging.info("Writing resolved direct-ONT SAM...")
    with input_path.open("r", encoding="utf-8") as fin, output_path.open("w", encoding="utf-8") as fout:
        buffer: List[str] = []
        for idx, line in enumerate(fin):
            if idx >= total or not keep[idx]:
                continue
            fields = _remove_existing_direct_tags(line.rstrip("\n").split("\t"))
            if len(fields) < 6:
                continue
            taxid = _taxid_from_ref(fields[2]) or "0"
            species = str(species_all[idx]) or taxid
            fields.extend(
                [
                    f"ZQ:i:{int(qstart_all[idx])}",
                    f"ZE:i:{int(qend_all[idx])}",
                    f"ZL:i:{int(qlen_all[idx])}",
                    f"ZT:Z:{taxid}",
                    f"ZS:Z:{species}",
                ]
            )
            if count_once[idx]:
                fields.append("ZC:i:1")
            buffer.append("\t".join(fields) + "\n")
            if len(buffer) >= 1000:
                fout.writelines(buffer)
                buffer.clear()
        if buffer:
            fout.writelines(buffer)

    return stats


# Backward-compatible name used by earlier direct-ONT prototypes.
def curate_ont_sam(
    input_sam: str | os.PathLike[str],
    output_sam: str | os.PathLike[str],
    config: Optional[OntDirectConfig] = None,
    temp_parent: Optional[str | os.PathLike[str]] = None,
    ambiguous_tsv: Optional[str | os.PathLike[str]] = None,
) -> OntDirectStats:
    del temp_parent
    return resolve_sam(
        input_sam=input_sam,
        output_sam=output_sam,
        config=config,
        ambiguous_tsv=ambiguous_tsv,
    )
