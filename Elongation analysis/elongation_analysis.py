# -*- coding: utf-8 -*-
"""Generate per‑gene elongation ratios (bin2 + bin3)/(bin1).

The script **no longer performs any statistical testing** – it simply
produces three output tables per BAM file:
  • raw bin counts
  • filtered counts (≥20 reads total)
  • filtered counts + normalised bins + elongation ratio
Those CSVs can be imported directly into DESeq2 / limma‑voom.

A few critical bug‑fixes were applied compared with the previous version
(see comments inline):
  – bins are reversed on minus‑strand CDS so that *bin1* is always the
    5′‑most third of the ORF.
  – reads that only *touch* a boundary are ignored (overlap must be > 0).
  – counting is keyed by **transcript_id** so isoforms are kept separate
    (change `use_transcript_id` to False if you prefer gene‑level).

Usage
-----
python elongation_ratio.py config.yaml
"""
from __future__ import annotations

import os
import sys
import yaml
import logging
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np
import pandas as pd
import pysam

##############################################################################
# logging & helpers
##############################################################################
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)

def load_config(cfg: str | os.PathLike) -> dict:
    with open(cfg, "r") as f:
        return yaml.safe_load(f)

##############################################################################
# directory scaffold
##############################################################################

def create_output_dirs(out: str) -> None:
    for sub in ("bin_counts", "filtered_counts", "normalized_bin_counts"):
        os.makedirs(Path(out, sub), exist_ok=True)

##############################################################################
# core functionality
##############################################################################

def divide_cds_into_bins(start: int, end: int) -> List[Tuple[int, int]]:
    """Return three equal‑length (start, end) tuples covering [start, end)."""
    edges = np.linspace(start, end, num=4, dtype=int)
    return [(int(edges[i]), int(edges[i + 1])) for i in range(3)]


def count_reads_in_bins(
    bam_path: str | os.PathLike,
    annotation_bed: str | os.PathLike,
    min_len: int,
    max_len: int,
    use_transcript_id: bool = True,
) -> Dict[str, List[int]]:
    """Return dict[id] -> [bin1, bin2, bin3] read counts."""
    counts: Dict[str, List[int]] = {}
    bam = pysam.AlignmentFile(bam_path, "rb")

    with open(annotation_bed) as bed:
        for line in bed:
            chrom, start, end, gene, _, strand = line.rstrip().split("\t")
            start, end = int(start), int(end)
            # always have bin1 = 5′, regardless of strand
            bins = divide_cds_into_bins(start, end)
            if strand == "-":
                bins = bins[::-1]

            bin_hits = [0, 0, 0]
            for read in bam.fetch(chrom, start, end):
                if not (min_len <= read.query_length <= max_len):
                    continue
                r_start = read.reference_start
                r_end = read.reference_end  # exclusive
                r_len = r_end - r_start
                for idx, (b_start, b_end) in enumerate(bins):
                    overlap = min(r_end, b_end) - max(r_start, b_start)
                    if overlap > 0 and overlap >= 0.5 * r_len:
                        bin_hits[idx] += 1
                        break  # count read once only

            if use_transcript_id:
                key = gene  # already transcript in BED; change if needed
            else:
                key = gene.split(";")[0]  # keep gene only (expects "gene;tx")

            if key in counts:
                counts[key] = [counts[key][i] + bin_hits[i] for i in range(3)]
            else:
                counts[key] = bin_hits
    bam.close()
    return counts


def normalize_and_write(
    counts: Dict[str, List[int]],
    *,
    outfile_raw: str,
    outfile_filtered: str,
    outfile_norm: str,
    min_total: int = 20,
) -> None:
    """Save three CSVs: raw, filtered, and normalised with elongation ratio."""
    df = pd.DataFrame.from_dict(counts, orient="index", columns=["bin1", "bin2", "bin3"])
    df.index.name = "id"
    df.to_csv(outfile_raw)

    filt = df[df.sum(axis=1) >= min_total].copy()
    filt.to_csv(outfile_filtered)

    totals = filt[["bin1", "bin2", "bin3"]].sum(axis=1)
    for b in ("bin1", "bin2", "bin3"):
        filt[f"norm_{b}"] = filt[b] / totals

    filt["elongation_ratio"] = (filt[["bin2", "bin3"]].sum(axis=1) + 1e-6) / (filt["bin1"] + 1e-6)
    filt.to_csv(outfile_norm)

##############################################################################
# entry‑point
##############################################################################

def main(cfg_file: str):
    cfg = load_config(cfg_file)
    outdir = Path(cfg["output_dir"]).resolve()
    create_output_dirs(outdir)

    ann = cfg["annotation_file"]
    min_len = cfg["min_read_length"]
    max_len = cfg["max_read_length"]

    for group, bams in cfg["samples"].items():
        logging.info("Processing condition: %s", group)
        grp_raw   = Path(outdir, "bin_counts", group)
        grp_filt  = Path(outdir, "filtered_counts", group)
        grp_norm  = Path(outdir, "normalized_bin_counts", group)
        for p in (grp_raw, grp_filt, grp_norm):
            p.mkdir(parents=True, exist_ok=True)

        for bam in bams:
            bam_path = Path(bam).resolve()
            logging.info("  Counting reads in %s", bam_path.name)
            counts = count_reads_in_bins(bam_path, ann, min_len, max_len)

            stem = bam_path.with_suffix("").name
            normalize_and_write(
                counts,
                outfile_raw=grp_raw   / f"{stem}_counts.csv",
                outfile_filtered=grp_filt  / f"{stem}_filtered.csv",
                outfile_norm=grp_norm / f"{stem}_norm.csv",
            )
            logging.info("  ↳ finished %s", stem)

    logging.info("All samples processed. You can now import the *_norm.csv files into DESeq2 / limma‑voom for downstream statistics.")

##############################################################################
if __name__ == "__main__":
    if len(sys.argv) != 2:
        logging.error("Usage: python elongation_ratio.py <config.yaml>")
        sys.exit(1)
    main(sys.argv[1])
