import os
import yaml
import pysam
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from multiprocessing import Pool, cpu_count
import logging
from tqdm import tqdm  # For progress indication

# Configure logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_config(config_file):
    with open(config_file, 'r') as file:
        return yaml.safe_load(file)

# Create output directory structure
def create_output_dirs(output_dir):
    os.makedirs(f"{output_dir}/bin_counts", exist_ok=True)
    os.makedirs(f"{output_dir}/filtered_counts", exist_ok=True)
    os.makedirs(f"{output_dir}/normalized_bin_counts", exist_ok=True)
    os.makedirs(f"{output_dir}/statistical_analysis/approach_1", exist_ok=True)
    os.makedirs(f"{output_dir}/statistical_analysis/approach_2", exist_ok=True)

# Function to divide CDS into 9 evenly spaced bins
def divide_cds_into_bins(start, end):
    bin_edges = np.linspace(start, end, num=4, dtype=int)  # Divide into 3 bins
    return [(bin_edges[i], bin_edges[i+1]) for i in range(3)]

# Count reads in each bin allowing partial overlap (counts towards the earlier bin if spanning multiple bins)
def count_reads_in_bins(bam_file, annotation_file, min_read_length, max_read_length):
    counts = {}
    bam = pysam.AlignmentFile(bam_file, 'rb')
    with open(annotation_file, 'r') as bed_file:
        for line in bed_file:
            chrom, start, end, gene, _, strand = line.strip().split('\t')
            start, end = int(start), int(end)
            bins = divide_cds_into_bins(start, end)
            bin_counts = [0] * 3
            for read in bam.fetch(chrom, start, end):
                if min_read_length <= read.query_length <= max_read_length:
                    read_start = read.reference_start
                    read_end = read.reference_end
                    for i, (bin_start, bin_end) in enumerate(bins):
                        overlap = min(read_end, bin_end) - max(read_start, bin_start)
                        if overlap > 0.5 * (read_end - read_start):  # Count read if it overlaps at least 50%
                            bin_counts[i] += 1
                            break  # Assign to the first matching bin
            counts[gene] = bin_counts
    bam.close()
    return counts

# Function to normalize bin counts
def normalize_bin_counts(df):
    total_counts = df.iloc[:, :3].sum(axis=1)
    for i in range(3):
        df[f'normalized_bin{i+1}'] = df[f'bin{i+1}'] / total_counts
    return df

# Save counts to CSV and normalize bin counts
def save_and_normalize_counts(counts, output_file, filtered_output_file, normalized_output_file):
    df = pd.DataFrame.from_dict(counts, orient='index', columns=[f"bin{i+1}" for i in range(3)])
    
    # Explicit summation to avoid incorrect sums
    df['sum_bins_2_3'] = df[['bin2', 'bin3']].sum(axis=1)
    df['elongation_ratio'] = (df['sum_bins_2_3'] + 1e-6) / (df['bin1'] + 1e-6)
    df.to_csv(output_file)
    
    # Filter genes with total counts < 20 across all bins
    filtered_df = df[df.iloc[:, :9].sum(axis=1) >= 20]
    filtered_df.to_csv(filtered_output_file)
    
    # Normalize bin counts and save
    normalized_df = normalize_bin_counts(filtered_df.copy())
    normalized_df.to_csv(normalized_output_file)

# Compute statistics for Mann-Whitney test
def compute_stat(args):
    gene, control_data, group_data = args
    try:
        control_values = control_data.loc[gene, 'elongation_ratio'].values
        exp_values = group_data.loc[gene, 'elongation_ratio'].values
        
        # Debugging log to inspect input values
        logging.info(f"Gene: {gene}, Control Values: {control_values}, Exp Values: {exp_values}")
                
        # Skip if not enough data
        if len(control_values) <= 1 or len(exp_values) <= 1:
            p_value = np.nan
            log2fc = np.nan
        else:
            _, p_value = mannwhitneyu(control_values, exp_values, alternative='two-sided')
            log2fc = np.log2(np.mean(exp_values) / np.mean(control_values)) if np.mean(control_values) > 0 else np.nan
        
        return [gene, log2fc, p_value]
    except Exception as e:
        logging.error(f"Error processing gene {gene}: {e}")
        return [gene, np.nan, np.nan]


# Statistical analysis for Approach 1
def statistical_analysis_approach_1(control_data, group_data, output_file):
    genes = control_data.index.intersection(group_data.index)
    all_results = []
    for gene in genes:
        try:
            control_values = control_data.loc[gene, 'elongation_ratio']
            exp_values = group_data.loc[gene, 'elongation_ratio']
            if np.isscalar(control_values):
                control_values = [control_values]
            if np.isscalar(exp_values):
                exp_values = [exp_values]
            _, p_value = mannwhitneyu(control_values, exp_values, alternative='two-sided')
            log2fc = np.log2(np.mean(exp_values) / np.mean(control_values)) if np.mean(control_values) > 0 else np.nan
            all_results.append([gene, log2fc, p_value])
        except Exception as e:
            logging.error(f"Error processing gene {gene} in Approach 1: {e}")
            all_results.append([gene, np.nan, np.nan])
    
    results_df = pd.DataFrame(all_results, columns=['gene', 'log2FC', 'p.value'])
    results_df = results_df.groupby('gene', as_index=False).mean()  # Aggregate by gene to remove duplicates
    results_df['FDR'] = multipletests(results_df['p.value'].fillna(1), method='fdr_bh')[1]
    results_df.to_csv(output_file, index=False)

#statistical approach 2 --> This approach was not used for downstream analysis
def process_bin(args):
    file, bin_name = args
    df = pd.read_csv(file, index_col=0)
    results = []
    for gene in df.index:
        control_values = df.loc[gene, bin_name]
        exp_values = df.loc[gene, bin_name]
        logging.info(f"Gene: {gene}, Bin: {bin_name}, Control: {control_values}, Exp: {exp_values}")
        try:
            if np.isnan(control_values).any() or np.isnan(exp_values).any() or np.mean(control_values) == 0 or np.mean(exp_values) == 0:
                p_value = np.nan
                log2fc = np.nan
            elif isinstance(control_values, (float, int)) and isinstance(exp_values, (float, int)) and control_values == exp_values:
                p_value = 1.0
                log2fc = 0.0
            else:
                if isinstance(control_values, (float, int)):
                    control_values = [control_values]
                if isinstance(exp_values, (float, int)):
                    exp_values = [exp_values]
                _, p_value = mannwhitneyu(control_values, exp_values, alternative='two-sided')
                log2fc = np.log2(np.mean(exp_values) / np.mean(control_values)) if np.mean(control_values) > 0 else np.nan
            results.append([gene, bin_name, p_value, log2fc])
        except Exception as e:
            logging.error(f"Error processing gene {gene} at {bin_name} in Approach 2: {e}")
            logging.info(f"Total number of tasks: {len(tasks)}")
            results.append([gene, bin_name, np.nan, np.nan])
    return results

def statistical_analysis_approach_2(normalized_files, output_file):
    all_results = []
    bins = [f'normalized_bin{i+1}' for i in range(3)]
    with Pool(processes=cpu_count()) as pool:
        tasks = [(file, bin_name) for file in normalized_files for bin_name in bins]
        for result in tqdm(pool.imap(process_bin, tasks), total=len(tasks), desc="Statistical Analysis (Approach 2)"):
            all_results.extend(result)

    results_df = pd.DataFrame(all_results, columns=['gene', 'bin', 'p.value', 'log2FC'])
    results_df['FDR'] = multipletests(results_df['p.value'].fillna(1), method='fdr_bh')[1]
    results_df.to_csv(output_file, index=False)


# Main function
def main(config_file):
    config = load_config(config_file)
    create_output_dirs(config['output_dir'])
    sample_groups = config['samples']
    min_read_length = config['min_read_length']
    max_read_length = config['max_read_length']
    annotation_file = config['annotation_file']
    output_dir = config['output_dir']
    
    # Step 1: Count reads in bins for each group
    for group, bam_files in sample_groups.items():
        group_output_dir = f"{output_dir}/bin_counts/{group}"
        filtered_output_dir = f"{output_dir}/filtered_counts/{group}"
        normalized_output_dir = f"{output_dir}/normalized_bin_counts/{group}"
        os.makedirs(group_output_dir, exist_ok=True)
        os.makedirs(filtered_output_dir, exist_ok=True)
        os.makedirs(normalized_output_dir, exist_ok=True)
        logging.info(f"Counting reads for group: {group}")
        
        for bam_file in bam_files:
            logging.info(f"Processing {bam_file}...")
            counts = count_reads_in_bins(bam_file, annotation_file, min_read_length, max_read_length)
            output_file = f"{group_output_dir}/{os.path.basename(bam_file).replace('.bam', '_counts.csv')}"
            filtered_output_file = f"{filtered_output_dir}/{os.path.basename(bam_file).replace('.bam', '_filtered_counts.csv')}"
            normalized_output_file = f"{normalized_output_dir}/{os.path.basename(bam_file).replace('.bam', '_normalized_counts.csv')}"
            save_and_normalize_counts(counts, output_file, filtered_output_file, normalized_output_file)
            logging.info(f"Saved count file: {output_file}, filtered count file: {filtered_output_file}, and normalized count file: {normalized_output_file}")
        
        # File monitoring
        logging.info(f"Finished counting for group: {group}. Count files located in: {group_output_dir}, filtered files in: {filtered_output_dir}, and normalized files in: {normalized_output_dir}")
    
    # Step 2: Statistical Analysis - Approach 1
    control_files = [f"{output_dir}/filtered_counts/control/{file}" for file in os.listdir(f"{output_dir}/filtered_counts/control") if file.endswith(".csv")]
    if not control_files:
        logging.error("No control sample CSV files found. Please check your read counting step.")
        sys.exit(1)
    control_data = pd.concat([pd.read_csv(file, index_col=0) for file in control_files])
    
    for group in sample_groups.keys():
        if group == config['control_condition']:
            continue
        group_files = [f"{output_dir}/filtered_counts/{group}/{file}" for file in os.listdir(f"{output_dir}/filtered_counts/{group}") if file.endswith(".csv")]
        group_data = pd.concat([pd.read_csv(file, index_col=0) for file in group_files])
        output_file = f"{output_dir}/statistical_analysis/approach_1/{group}_vs_control.csv"
        logging.info(f"Performing Approach 1 statistical analysis for group: {group}")
        statistical_analysis_approach_1(control_data, group_data, output_file)
        logging.info(f"Saved Approach 1 results for group: {group} to {output_file}")
    
    # Step 3: Statistical Analysis - Approach 2 --> Ignore the output of this approach
    for group in sample_groups.keys():
        if group == config['control_condition']:
            continue
        group_files = [f"{output_dir}/normalized_bin_counts/{group}/{file}" for file in os.listdir(f"{output_dir}/normalized_bin_counts/{group}") if file.endswith(".csv")]
        output_file = f"{output_dir}/statistical_analysis/approach_2/{group}_vs_control.csv"
        logging.info(f"Performing Approach 2 statistical analysis for group: {group}")
        statistical_analysis_approach_2(group_files, output_file)
        logging.info(f"Saved Approach 2 results for group: {group} to {output_file}")

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 2:
        logging.error("Usage: python elongation_analysis.py <config.yaml>")
        sys.exit(1)
    main(sys.argv[1])
