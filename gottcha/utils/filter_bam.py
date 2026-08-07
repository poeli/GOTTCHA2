#!/usr/bin/env python3
"""Write alignments that meet logged GOTTCHA matching criteria to a BAM file."""

from __future__ import annotations

import argparse
import math
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional, TextIO

import pysam


@dataclass(frozen=True)
class MatchCriteria:
    """Alignment thresholds recorded in a GOTTCHA log file."""

    min_identity: float
    min_fraction: float
    min_length: int


@dataclass(frozen=True)
class FilterStats:
    """Counts produced while filtering a BAM file."""

    total: int
    kept: int
    filtered: int


_CRITERION_PATTERNS = {
    "min_identity": re.compile(
        r"\bmin(?:imum)?\s+match\s+identity\s*:\s*(\S+)", re.IGNORECASE
    ),
    "min_fraction": re.compile(
        r"\bmin(?:imum)?\s+match\s+fraction\s*:\s*(\S+)", re.IGNORECASE
    ),
    "min_length": re.compile(
        r"\bmin(?:imum)?\s+match\s+length\s*:\s*(\S+)", re.IGNORECASE
    ),
}


def load_criteria_from_log(log_path: str) -> MatchCriteria:
    """Load the most recently logged match thresholds from a GOTTCHA log."""
    values = {}
    with open(log_path, encoding="utf-8", errors="replace") as log_file:
        for line_number, line in enumerate(log_file, start=1):
            for name, pattern in _CRITERION_PATTERNS.items():
                match = pattern.search(line)
                if match:
                    values[name] = (match.group(1), line_number)

    missing = [name for name in _CRITERION_PATTERNS if name not in values]
    if missing:
        labels = ", ".join(name.replace("_", " ") for name in missing)
        raise ValueError(f"Missing matching criteria in {log_path}: {labels}")

    try:
        min_identity = float(values["min_identity"][0])
        min_fraction = float(values["min_fraction"][0])
        min_length = int(values["min_length"][0])
    except ValueError as exc:
        raise ValueError(f"Invalid matching criterion in {log_path}: {exc}") from exc

    if not math.isfinite(min_identity) or not 0 <= min_identity <= 1:
        raise ValueError(f"Minimum match identity must be between 0 and 1: {min_identity}")
    if not math.isfinite(min_fraction) or not 0 <= min_fraction <= 1:
        raise ValueError(f"Minimum match fraction must be between 0 and 1: {min_fraction}")
    if min_length < 0:
        raise ValueError(f"Minimum match length must be non-negative: {min_length}")

    return MatchCriteria(min_identity, min_fraction, min_length)


def alignment_meets_criteria(
    alignment: pysam.AlignedSegment,
    reference_length: int,
    criteria: MatchCriteria,
) -> bool:
    """Return whether an alignment passes the three GOTTCHA match thresholds."""
    # AlignedSegment.alen is the length used by process_bam.py. It consumes
    # reference bases (including deletions) and excludes insertions/soft clips.
    aligned_length = alignment.alen

    if criteria.min_identity > 0 and alignment.has_tag("NM"):
        if aligned_length <= 0:
            return False
        identity = 1 - (alignment.get_tag("NM") / aligned_length)
        if identity < criteria.min_identity:
            return False

    if criteria.min_fraction > 0:
        query_length = alignment.query_length or 0
        query_fraction = aligned_length / query_length if query_length > 0 else 0.0
        reference_fraction = (
            aligned_length / reference_length if reference_length > 0 else 0.0
        )
        if max(query_fraction, reference_fraction) < criteria.min_fraction:
            return False

    return aligned_length >= criteria.min_length


def alignment_is_qualified(
    alignment: pysam.AlignedSegment,
    reference_length: int,
    criteria: MatchCriteria,
    min_mapq: int = 0,
    include_secondary: bool = False,
    include_supplementary: bool = False,
    include_duplicates: bool = False,
    include_qcfail: bool = False,
) -> bool:
    """Apply the standard GOTTCHA BAM filters and logged match thresholds."""
    if alignment.is_unmapped:
        return False
    if not include_secondary and alignment.is_secondary:
        return False
    if not include_supplementary and alignment.is_supplementary:
        return False
    if not include_duplicates and alignment.is_duplicate:
        return False
    if not include_qcfail and alignment.is_qcfail:
        return False
    if alignment.mapping_quality < min_mapq:
        return False

    return alignment_meets_criteria(alignment, reference_length, criteria)


def filter_bam(
    input_bam: str,
    output_bam: str,
    criteria: MatchCriteria,
    min_mapq: int = 0,
    threads: int = 1,
    include_secondary: bool = False,
    include_supplementary: bool = False,
    include_duplicates: bool = False,
    include_qcfail: bool = False,
) -> FilterStats:
    """Stream an input BAM and write only qualified alignments to the output BAM."""
    input_path = Path(input_bam).resolve()
    output_path = Path(output_bam).resolve()
    if input_path == output_path:
        raise ValueError("Input and output BAM paths must be different")
    if threads < 1:
        raise ValueError("Threads must be at least 1")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    total = 0
    kept = 0

    with pysam.AlignmentFile(str(input_path), "rb", threads=threads) as source:
        reference_lengths = dict(zip(source.references, source.lengths))
        with pysam.AlignmentFile(
            str(output_path), "wb", template=source, threads=threads
        ) as destination:
            for alignment in source.fetch(until_eof=True):
                total += 1
                reference_length = reference_lengths.get(alignment.reference_name, 0)
                if alignment_is_qualified(
                    alignment,
                    reference_length,
                    criteria,
                    min_mapq=min_mapq,
                    include_secondary=include_secondary,
                    include_supplementary=include_supplementary,
                    include_duplicates=include_duplicates,
                    include_qcfail=include_qcfail,
                ):
                    destination.write(alignment)
                    kept += 1

    return FilterStats(total=total, kept=kept, filtered=total - kept)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Filter a BAM using the minimum match identity, fraction, and length "
            "recorded in a GOTTCHA log file."
        )
    )
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("log", help="GOTTCHA log containing match criteria")
    parser.add_argument("-o", "--output", required=True, help="Output BAM file")
    parser.add_argument(
        "--min-mapq", type=int, default=0, help="Minimum MAPQ (default: 0)"
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=1,
        help="HTSlib compression/decompression threads (default: 1)",
    )
    parser.add_argument(
        "--include-secondary", action="store_true", help="Include secondary alignments"
    )
    parser.add_argument(
        "--include-supplementary",
        action="store_true",
        help="Include supplementary alignments",
    )
    parser.add_argument(
        "--include-duplicates",
        action="store_true",
        help="Include duplicate-marked alignments",
    )
    parser.add_argument(
        "--include-qcfail", action="store_true", help="Include QC-failed alignments"
    )
    return parser


def main(argv: Optional[List[str]] = None, stderr: Optional[TextIO] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    err = stderr or sys.stderr

    try:
        criteria = load_criteria_from_log(args.log)
        stats = filter_bam(
            args.bam,
            args.output,
            criteria,
            min_mapq=args.min_mapq,
            threads=args.threads,
            include_secondary=args.include_secondary,
            include_supplementary=args.include_supplementary,
            include_duplicates=args.include_duplicates,
            include_qcfail=args.include_qcfail,
        )
    except (OSError, ValueError, pysam.utils.SamtoolsError) as exc:
        print(f"ERROR: {exc}", file=err)
        return 2

    print(
        f"Loaded criteria: identity>={criteria.min_identity}, "
        f"fraction>={criteria.min_fraction}, length>={criteria.min_length}",
        file=err,
    )
    print(
        f"Wrote {stats.kept:,} of {stats.total:,} alignments to {args.output} "
        f"({stats.filtered:,} filtered)",
        file=err,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
