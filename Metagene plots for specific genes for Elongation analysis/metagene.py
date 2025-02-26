import os
import yaml
import sys
import pandas as pd
import matplotlib.pyplot as plt
import pysam
from collections import defaultdict
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed

# Load config.yaml
with open("config.yaml", "r") as f:
    config = yaml.safe_load(f)

output_dir = config["output_dir"]
annotation_file = config["annotation_file"]
gene_list_file = config["gene_list"]
max_threads = config["max_threads"] if config["max_threads"] != "auto" else os.cpu_count()

# Ensure output directory exists
os.makedirs(output_dir, exist_ok=True)

def get_gene_coordinates(gene_name, annotation_file):
    """Parse the GTF file to get the coordinates of the specified gene."""
    with open(annotation_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "transcript":  # Focus on transcript entries
                attributes = fields[8]
                if f'gene_name "{gene_name}"' in attributes:  # Search for gene_name pattern
                    return fields[0], int(fields[3]), int(fields[4])  # Return chromosome, start, end
    raise ValueError(f"Gene {gene_name} not found in {annotation_file}.")

def get_coverage(bam_file, chrom, start, end):
    """Extract coverage for the specified region using pysam."""
    coverage = np.zeros(end - start + 1)
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        for pileupcolumn in bam.pileup(chrom, start, end, truncate=True):
            if start <= pileupcolumn.reference_pos <= end:
                coverage[pileupcolumn.reference_pos - start] = pileupcolumn.nsegments
    return coverage

def normalize_coverage(coverage, bam_file):
    """Normalize coverage to total mapped reads in the BAM file."""
    with pysam.AlignmentFile(bam_file, "rb") as bam:
        total_reads = bam.mapped
    return coverage / total_reads if total_reads > 0 else coverage

def average_coverage(group_coverages):
    """Compute the average coverage across all samples in a group."""
    return np.mean(group_coverages, axis=0)

def plot_metagene(normalized_coverages, gene_name, output_file):
    """Generate the metagene plot."""
    plt.figure(figsize=(12, 6))
    plt.rcParams.update({
        'font.size': 16,        # General font size
        'axes.titlesize': 20,   # Title font size
        'axes.labelsize': 18,   # Axis label font size
        'xtick.labelsize': 14,  # X-axis tick label size
        'ytick.labelsize': 14,  # Y-axis tick label size
        'legend.fontsize': 14   # Legend font size
    })
    for group, coverage in normalized_coverages.items():
        plt.plot(np.linspace(0, 100, len(coverage)), coverage, label=group)
    plt.xlabel("Normalized Gene Position (%)")
    plt.ylabel("Average Normalized Coverage")
    plt.title(f"Metagene Plot for {gene_name}")
    plt.legend()
    plt.grid(True)
    plt.savefig(output_file)
    plt.close()

def process_gene(gene_name):
    """Process a single gene to generate its metagene plot."""
    print(f"Processing gene: {gene_name}")
    try:
        chrom, start, end = get_gene_coordinates(gene_name, annotation_file)
    except ValueError as e:
        print(e)
        return

    group_coverages = defaultdict(list)

    # Process each group in the config
    for group, bam_files in config["samples"].items():
        for bam_file in bam_files:
            print(f"  Processing {bam_file} for group {group}...")
            coverage = get_coverage(bam_file, chrom, start, end)
            normalized = normalize_coverage(coverage, bam_file)
            group_coverages[group].append(normalized)
    
    # Compute average normalized coverage for each group
    averaged_coverages = {group: average_coverage(coverages) for group, coverages in group_coverages.items()}

    # Plot and save the metagene plot
    output_file = os.path.join(output_dir, f"{gene_name}_metagene_plot.png")
    plot_metagene(averaged_coverages, gene_name, output_file)
    print(f"  Metagene plot saved at {output_file}")

def main():
    # Read the list of genes from the gene list file
    with open(gene_list_file, "r") as f:
        genes = [line.strip() for line in f if line.strip()]
    
    print(f"Found {len(genes)} genes to process.")
    
    # Use ThreadPoolExecutor for parallel processing
    with ThreadPoolExecutor(max_workers=max_threads) as executor:
        future_to_gene = {executor.submit(process_gene, gene): gene for gene in genes}
        
        for future in as_completed(future_to_gene):
            gene = future_to_gene[future]
            try:
                future.result()
            except Exception as exc:
                print(f"Gene {gene} generated an exception: {exc}")

if __name__ == "__main__":
    main()
