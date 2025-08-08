import yaml
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from multiprocessing import Pool, cpu_count
import os
import re
from pathlib import Path

# Load Config.yaml
with open("Config.yaml", "r") as config_file:
    config = yaml.safe_load(config_file)

# Settings
control_condition = config['control_condition']
max_threads = cpu_count() if config['max_threads'] == 'auto' else int(config['max_threads'])
cds_bed = config['annotations']['map_CDS']
utr_bed = config['annotations']['map_3UTR']

# Create output directories
os.makedirs("SCRT_counts", exist_ok=True)

# Function to count reads using bedtools intersect
def count_reads(sample, bed_file):
    intersect_cmd = f"bedtools intersect -a {bed_file} -b {sample} -c"
    result = subprocess.run(intersect_cmd, shell=True, capture_output=True, text=True)
    counts = {}
    for line in result.stdout.strip().split("\n"):
        fields = line.split("\t")
        gene = fields[3].split(';')[0]  # Extract gene name
        count = int(fields[-1])
        counts[gene] = counts.get(gene, 0) + count
    return counts

# Function to process each sample
def process_sample(sample, sample_type):
    cds_counts = count_reads(sample, cds_bed)
    utr_counts = count_reads(sample, utr_bed)
    output_file = f"SCRT_counts/{sample_type}_{os.path.basename(sample).replace('.bed', '')}_counts.csv"
    with open(output_file, "w") as out:
        out.write("Gene,CDS_counts,UTR_counts,RRS\n")
        for gene in set(cds_counts.keys()).union(utr_counts.keys()):
            cds_count = cds_counts.get(gene, 0)
            utr_count = utr_counts.get(gene, 0)
            rrs = (utr_count + 1e-6) / (cds_count + 1e-6)
            out.write(f"{gene},{cds_count},{utr_count},{rrs}\n")
    return output_file

# Parallel processing of samples
def process_all_samples():
    sample_files = []
    for condition, samples in config['samples'].items():
        for sample in samples:
            sample_files.append((sample, condition))
    with Pool(max_threads) as pool:
        result_files = pool.starmap(process_sample, sample_files)
    return result_files

# Merge and filter to build a matrix of RRS values
def tidy_sample_name(csv_path: str) -> str:
    # "control_Sham1_counts.csv"  ->  "control_Sham1"
    stem = os.path.basename(csv_path).replace("_counts.csv", "")
    stem = re.sub(r"\.bed$", "", stem)
    return stem

def build_rrs_matrix(csv_files, min_reads=100, out_csv="SCRT_counts/RRS_matrix.csv"):
    dfs      = []
    for fp in csv_files:
        df         = pd.read_csv(fp)
        df["Total"] = df["CDS_counts"] + df["UTR_counts"]
        df          = df[df["Total"] >= min_reads]                  # sample-specific filter
        sample_name = tidy_sample_name(fp)
        dfs.append(df[["Gene", "RRS"]].rename(columns={"RRS": sample_name}))

    # successive outer-joins, then drop genes missing in any sample
    merged = dfs[0]
    for d in dfs[1:]:
        merged = merged.merge(d, on="Gene", how="outer")

    merged.dropna(inplace=True)  # keep genes that passed the filter in *all* samples
    merged.to_csv(out_csv, index=False)
    print(f"✓ merged matrix written to {out_csv}")

# Main execution
if __name__ == "__main__":
    print("Processing samples...")
    result_files = process_all_samples()
    build_rrs_matrix(result_files)   

