#!/usr/bin/env python3
"""
sam_to_bam.py

Convert a headerless GOTTCHA2 SAM file to a coordinate-sorted, indexed BAM.

Optimizations over the previous implementation:
- Streams the SAM to collect unique RNAME/RNEXT references instead of loading
  columns into pandas.
- Excludes '*'/ '=' from both reference columns.
- Builds only the @SQ records actually required by the alignments.
- Pipes "header + SAM" directly into `samtools sort`, eliminating the
  intermediate BAM and its compression/decompression cycle.
- Uses `samtools index -@` directly; pysam is not required.

Expected GOTTCHA2 reference-name format:
    <name>|<start>|<end>|...

where <start> and <end> are used to derive the @SQ LN value.
"""

import argparse
import logging
import os
import subprocess
import sys
import tempfile
from collections import OrderedDict
from typing import Iterable, MutableMapping, Union


# Large userspace read buffer helps when scanning very large SAM files.
SAM_BUFFER_SIZE = 16 * 1024 * 1024


def _add_reference(
    refs: MutableMapping[bytes, None],
    ref: bytes,
) -> None:
    """Add a valid SAM reference name while preserving first-seen order."""
    if ref not in (b"*", b"="):
        refs.setdefault(ref, None)


def collect_references(input_sam: str) -> list[bytes]:
    """
    Extract unique references from both RNAME (column 3) and RNEXT (column 7).

    The input is read as bytes to avoid unnecessary Unicode decoding and object
    creation. Only the first seven SAM fields are parsed.

    Parameters
    ----------
    input_sam
        Headerless or header-containing SAM path.

    Returns
    -------
    list[bytes]
        Unique reference names in first-seen order.
    """
    refs: "OrderedDict[bytes, None]" = OrderedDict()

    with open(input_sam, "rb", buffering=SAM_BUFFER_SIZE) as sam_fh:
        for line_number, line in enumerate(sam_fh, start=1):
            if not line:
                continue

            # Be tolerant if an input SAM still contains headers.
            if line.startswith(b"@"):
                continue

            # We need fields 3 and 7 only. maxsplit=7 avoids parsing the rest
            # of long SAM records (SEQ, QUAL, optional tags, etc.).
            fields = line.rstrip(b"\r\n").split(b"\t", 7)

            if len(fields) < 7:
                raise ValueError(
                    f"Malformed SAM record at line {line_number}: "
                    f"expected at least 7 tab-delimited fields"
                )

            _add_reference(refs, fields[2])  # RNAME
            _add_reference(refs, fields[6])  # RNEXT

    return list(refs.keys())


def reference_length(ref: bytes) -> int:
    """
    Derive fragment length from a GOTTCHA2 reference name.

    Expected form includes:
        ...|START|END|...

    Returns END - START + 1.
    """
    parts = ref.split(b"|")

    if len(parts) < 3:
        ref_text = ref.decode("utf-8", errors="replace")
        raise ValueError(
            f"Cannot derive reference length from '{ref_text}': "
            "expected a reference name containing '|start|end|'"
        )

    try:
        start = int(parts[1])
        end = int(parts[2])
    except ValueError as exc:
        ref_text = ref.decode("utf-8", errors="replace")
        raise ValueError(
            f"Cannot derive reference length from '{ref_text}': "
            "start/end fields are not integers"
        ) from exc

    length = end - start + 1
    if length <= 0:
        ref_text = ref.decode("utf-8", errors="replace")
        raise ValueError(
            f"Invalid reference coordinates in '{ref_text}': "
            f"start={start}, end={end}"
        )

    return length


def write_header(refs: Iterable[bytes], header_file: str) -> None:
    """Write a minimal SAM header for the references used by the input SAM."""
    with open(header_file, "wb", buffering=1024 * 1024) as out:
        # The records being supplied to samtools sort are not yet sorted.
        # samtools sort will set SO:coordinate in the BAM it produces.
        out.write(b"@HD\tVN:1.6\tSO:unsorted\n")

        for ref in refs:
            length = reference_length(ref)
            out.write(
                b"@SQ\tSN:"
                + ref
                + b"\tLN:"
                + str(length).encode("ascii")
                + b"\n"
            )


def _wait_for_pipeline(
    producer: subprocess.Popen,
    consumer: subprocess.Popen,
    producer_name: str,
    consumer_cmd: list[str],
) -> None:
    """
    Wait for a two-process pipeline and raise an informative error on failure.
    """
    consumer_rc = consumer.wait()
    producer_rc = producer.wait()

    if consumer_rc != 0:
        raise subprocess.CalledProcessError(consumer_rc, consumer_cmd)

    # A producer may receive SIGPIPE when the consumer exits early. If the
    # consumer succeeded, however, a non-zero producer exit is still abnormal.
    if producer_rc != 0:
        raise subprocess.CalledProcessError(producer_rc, producer_name)


def convert_sam_to_bam(
    input_sam: str,
    output_bam: str,
    threads: int = 4,
    quiet: bool = False,
    compression_level: int = 1,
) -> None:
    """
    Convert SAM to coordinate-sorted BAM with a minimal reference header.

    The SAM is scanned once to determine the references needed for @SQ records.
    The generated header and original SAM are then concatenated as a stream and
    fed directly to samtools sort. No intermediate BAM is written.

    Parameters
    ----------
    input_sam
        Input SAM file. It may be headerless; any existing header lines are
        ignored when collecting references.
    output_bam
        Output coordinate-sorted BAM path.
    threads
        Number of threads passed to samtools sort/index.
    quiet
        Suppress samtools stderr output.
    compression_level
        BAM compression level for samtools sort (0-9). Level 1 is a good
        speed/size tradeoff for GOTTCHA2 intermediate/output BAMs.
    """
    if threads < 1:
        raise ValueError("threads must be >= 1")

    if not 0 <= compression_level <= 9:
        raise ValueError("compression_level must be between 0 and 9")

    if not os.path.isfile(input_sam):
        raise FileNotFoundError(f"Input SAM does not exist: {input_sam}")

    logging.info("Extracting references from %s...", input_sam)
    refs = collect_references(input_sam)
    logging.info("Unique references extracted: %d", len(refs))

    if not refs:
        logging.warning(
            "No mapped references were found in %s. "
            "A BAM with no @SQ records will be attempted.",
            input_sam,
        )

    output_dir = os.path.dirname(os.path.abspath(output_bam))
    os.makedirs(output_dir, exist_ok=True)

    stderr_target: Union[int, None]
    stderr_target = subprocess.DEVNULL if quiet else None

    with tempfile.TemporaryDirectory(prefix="gottcha2_sam2bam_") as temp_dir:
        header_file = os.path.join(temp_dir, "header.sam")

        logging.info("Creating minimal SAM header...")
        write_header(refs, header_file)

        # Stream:
        #   header.sam + input.sam -> samtools sort -> output.bam
        #
        # This avoids:
        #   samtools view -b -> temp.bam -> samtools sort
        logging.info("Sorting directly from SAM to BAM...")

        cat_cmd = ["cat", header_file, input_sam]
        sort_cmd = ["samtools", "sort", 
            "-@", str(threads),
            "-l", str(compression_level),
            "-o", output_bam,
            "-"
        ]

        logging.debug("Producer command: %s", " ".join(cat_cmd))
        logging.debug("Sort command: %s", " ".join(sort_cmd))

        producer = subprocess.Popen(
            cat_cmd,
            stdout=subprocess.PIPE,
            stderr=stderr_target,
        )

        try:
            consumer = subprocess.Popen(
                sort_cmd,
                stdin=producer.stdout,
                stderr=stderr_target,
            )
        except Exception:
            producer.kill()
            producer.wait()
            raise

        # Close the parent's copy so SIGPIPE is delivered correctly if
        # samtools sort exits.
        assert producer.stdout is not None
        producer.stdout.close()

        _wait_for_pipeline(
            producer,
            consumer,
            producer_name="cat",
            consumer_cmd=sort_cmd,
        )

    logging.info("Creating BAM index...")

    index_cmd = [
        "samtools",
        "index",
        "-@",
        str(threads),
        output_bam,
    ]
    logging.debug("Index command: %s", " ".join(index_cmd))

    try:
        subprocess.run(
            index_cmd,
            check=True,
            stderr=stderr_target,
        )
    except subprocess.CalledProcessError as exc:
        # Preserve the old behavior: BAM conversion can still succeed even if
        # indexing fails.
        logging.warning(
            "Could not create BAM index for %s (exit code %s)",
            output_bam,
            exc.returncode,
        )

    logging.info("Conversion complete: %s", output_bam)


def main(argv: list[str]) -> None:
    parser = argparse.ArgumentParser(
        prog="gottcha2 sam2bam",
        description=(
            "Convert a GOTTCHA2 SAM to a coordinate-sorted BAM using only "
            "references observed in RNAME/RNEXT."
        ),
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input SAM file",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=True,
        help="Output coordinate-sorted BAM file",
    )
    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads for samtools sort/index (default: 4)",
    )
    parser.add_argument(
        "-q",
        "--quiet",
        action="store_true",
        help="Suppress samtools warning/error output",
    )
    parser.add_argument(
        "-l",
        "--compression-level",
        type=int,
        default=1,
        choices=range(0, 10),
        metavar="0-9",
        help="BAM compression level used by samtools sort (default: 1)",
    )

    args = parser.parse_args(argv)

    try:
        convert_sam_to_bam(
            input_sam=args.input,
            output_bam=args.output,
            threads=args.threads,
            quiet=args.quiet,
            compression_level=args.compression_level,
        )
        logging.info("SAM to BAM conversion successful.")
    except Exception as exc:
        logging.error("Error converting file: %s", exc)
        sys.exit(1)


if __name__ == "__main__":
    main(sys.argv[1:])