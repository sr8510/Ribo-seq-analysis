import yaml
import subprocess
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from multiprocessing import Pool, cpu_count
import os

# Load Config.yaml
with open("Config.yaml", "r") as config_file:
    config = yaml.safe_load(config_file)

# Settings
control_condition = config['control_condition']
max_threads = cpu_count() if config['max_threads'] == 'auto' else int(config['max_threads'])
cds_bed = config['annotations']['map_CDS']
utr_bed = config['annotations']['map_3UTR']

# Create output directories
os.makedirs("counts", exist_ok=True)
os.makedirs("results", exist_ok=True)

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
    output_file = f"counts/{sample_type}_{os.path.basename(sample).replace('.bed', '')}_counts.csv"
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

# Function to perform statistical analysis for each condition vs control
def statistical_analysis(result_files):
    sample_data = {}
    for result_file in result_files:
        df = pd.read_csv(result_file)
        condition = os.path.basename(result_file).split('_')[0]
        if condition not in sample_data:
            sample_data[condition] = {}
        for _, row in df.iterrows():
            gene = row['Gene']
            if gene not in sample_data[condition]:
                sample_data[condition][gene] = []
            sample_data[condition][gene].append(row['RRS'])
    
    for condition in sample_data.keys():
        if condition == control_condition:
            continue
        gene_stats = {}
        for gene in set(g for cond in sample_data for g in sample_data[cond]):
            control_rrs = sample_data[control_condition].get(gene, [0])
            exp_rrs = sample_data[condition].get(gene, [0])
            
            # Check if all values are identical
            if len(set(control_rrs + exp_rrs)) == 1:
                p_value = 1  # Assign p-value of 1 for identical values
            else:
                _, p_value = mannwhitneyu(control_rrs, exp_rrs, alternative='two-sided')
            
            # Calculate fold change and then log2 transformation
            fold_change = (np.mean(exp_rrs) + 1e-6) / (np.mean(control_rrs) + 1e-6)
            log2fc = np.log2(fold_change)
            
            gene_stats[gene] = {
                "log2FC": log2fc,
                "p_value": p_value,
                "avg_RRS_experiment": np.mean(exp_rrs),
                "avg_RRS_control": np.mean(control_rrs)
            }
        
        stats_df = pd.DataFrame.from_dict(gene_stats, orient='index')
        stats_df.to_csv(f"results/statistical_analysis_{condition}_vs_control.csv", index_label="Gene")

# Main execution
if __name__ == "__main__":
    print("Processing samples...")
    result_files = process_all_samples()
    print("Performing statistical analysis...")
    statistical_analysis(result_files)
    print("Analysis complete. Results saved in 'results/statistical_analysis.csv'")

