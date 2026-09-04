from pathlib import Path
import re
import subprocess
import logging
from typing import List, Tuple
import pandas as pd
import logging
import pandas as pd

def minimap2(
    reads: List,
    db: str,
    threads: int,
    mm_options: str,
    presetx: str,
    samfile: Path,
    logfile: Path,
    allow_secondary: bool = False,
    max_secondary: int = 50,
    secondary_ratio: float = 0.5,
) -> Tuple[int, str, int, bool]:
    """Map reads to the GOTTCHA2 signature reference with minimap2.

    Direct ONT mode requests primary, supplementary and a bounded set of
    secondary candidate alignments. The command does not request supplementary
    soft-clipping or minimap2 split-index prefix handling; the direct resolver
    reconstructs original-read coordinates from S/H CIGAR
    clipping and does not use MAPQ for species consistency.
    """
    input_file = " ".join(reads)
    mapped_re = re.compile(r"mapped (\d+) sequences")
    multi_part_index_flag = False
    input_read_count = 0

    opts = [f"-x {presetx}"]
    if mm_options and mm_options.strip():
        opts.append(mm_options.strip())
    opts.extend(["-a", "--eqx", "--sam-hit-only"])
    if allow_secondary:
        opts.extend([
            f"-N{max(0, int(max_secondary))}",
            "--secondary=yes",
            f"-p{float(secondary_ratio):g}",
        ])
    else:
        opts.extend(["-N20", "--secondary=no"])

    mm2_cmd = f"minimap2 {' '.join(opts)} -t{threads} {db} {input_file}"
    filter_cmd = ['samtools', 'view', '-x', 'SA']

    with samfile.open("w", encoding="utf-8") as out_f:
        mm2 = subprocess.Popen(
            mm2_cmd,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )
        filter = subprocess.Popen(
            filter_cmd,
            stdin=mm2.stdout,
            stdout=out_f,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
        )

        mm2.stdout.close()
        with logfile.open("a", encoding="utf-8") as f:
            for line in mm2.stderr:
                if "For a multi-part index" in line:
                    multi_part_index_flag = True
                match = mapped_re.search(line)
                if match:
                    logging.debug(line)
                    input_read_count += int(match.group(1))
                f.write(line)

        mm2.stderr.close()
        filter.stderr.close()
        rc_mm = mm2.wait()
        filter.wait()

    return rc_mm, mm2_cmd, input_read_count, multi_part_index_flag


def post_processing_sam(samfile: Path, samfile_temp: Path) -> Tuple[bool, int, int]:
    """
    Removing multiple hits from the SAM file by keeping only the best alignment for each read.

    Parameters:
        samfile (str): Path to the SAM file
        samfile_temp (str): Path to the temporary SAM file with only the best alignments

    Returns:
        Tuple[bool, int, int]: (
            multiple_hits_removed (bool): False if no multiple hits were found, True if multiple hits were removed,
            total_alignments (int): Total number of alignments in the SAM file,
            top_score_hits (int): Number of top score hits after filtering
        )
    """
    logging.info(f'Loading the SAM file...')

    df = pd.read_csv(samfile,
                     sep='\t',
                     header=None,
                     usecols=[0, 1, 13],
                     names=['QNAME', 'FLAG', 'AS'],
                     converters={
                         'AS': lambda x: x.replace('AS:i:', '')
                     },
                     dtype={'QNAME': 'str', 'FLAG': 'uint16'}
    )

    aln_count = len(df)
    logging.info(f'Total alignments in SAM file: {aln_count}')

    df[['AS']] = df[['AS']].astype('int16')

    logging.info(f'Filtering non-primary hits...')
    # for each row, if the flag bitwise AND with 256 (not primary alignment) or 2048 (supplementary), then remove them from the df
    df = df[~(df['FLAG'] & (256|2048)).astype(bool)]
    logging.info(f'After removing non-primary hits: {len(df)}')

    logging.info(f'Identifying top score hits...')
    # if FLAG bitwise AND with 128 (second in pair), append '/2' to the QNAME
    idx = (df['FLAG'] & 128).astype(bool)
    df.loc[idx, 'QNAME'] = df.loc[idx, 'QNAME'] + '/2'

    # get the index with the best alignment score for each read
    idxmax = df.groupby('QNAME')['AS'].idxmax()
    logging.info(f'Total top score hits: {len(idxmax)}')

    if len(idxmax) == aln_count:
        logging.info(f'No multiple hits found. Keeping the original SAM file.')
        return False, aln_count, aln_count
    else:
        # Create a set of indices for faster lookup
        idxmax_set = set(idxmax.values)
        del idxmax

        logging.info(f'Writing top score hits...')
        with samfile_temp.open("w", encoding="utf-8") as fout, samfile.open("r", encoding="utf-8") as fin:
            lines_to_write = []
            for idx, line in enumerate(fin):
                if idx in idxmax_set:
                    lines_to_write.append(line)
                    if len(lines_to_write) >= 1000:
                        fout.writelines(lines_to_write)
                        lines_to_write.clear()
                        logging.debug(f'Written {idx} alignments...')

            if lines_to_write:
                fout.writelines(lines_to_write)
        logging.info(f'{len(idxmax_set)} hits written to {samfile_temp}.')

        return True, aln_count, len(idxmax_set)
