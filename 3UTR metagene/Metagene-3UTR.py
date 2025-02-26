import yaml
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from multiprocessing import Pool, cpu_count
import logging
import os

# Configure logging
logging.basicConfig(
    filename="metagene_debug.log",
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
)

# Load Config.yaml
with open("Config.yaml") as file:
    config = yaml.safe_load(file)

# Function to read bedgraph and calculate total coverage (for RPM normalization)
def read_bedgraph(file):
    logging.info(f"Reading {file}")
    df = pd.read_csv(file, sep="\t", header=None, names=["chrom", "start", "end", "value"])
    df["coverage"] = df["value"]
    df["position"] = (df["start"] + df["end"]) // 2  # Midpoint of each region
    df["length"] = df["end"] - df["start"]
    return df, df["coverage"].sum()

# Normalize by RPM (Reads Per Million)
def normalize_rpm(df, total_coverage):
    if total_coverage > 0:
        df["value"] = (df["coverage"] / total_coverage) * 1e6
    return df

# Corrected UTR length calculation and normalization of position as a percentage of total 3'UTR length
def normalize_position(df):
    utr_length = df["length"].sum()  # Sum of all region lengths
    if utr_length > 0:
        df["relative_position"] = ((df["position"] - df["start"].min()) / utr_length) * 100
    else:
        df["relative_position"] = np.nan
    return df

# Function to calculate average normalized coverage for each group with extended coverage
def calculate_average_coverage(group_files):
    normalized_dfs = []

    for file in group_files:
        df, total_coverage = read_bedgraph(file)
        df = normalize_rpm(df, total_coverage)
        df = normalize_position(df)

        if "relative_position" in df.columns:
            normalized_dfs.append(df[["relative_position", "value"]])

    if len(normalized_dfs) == 0:
        logging.warning("No valid data for group normalization")
        return pd.DataFrame(columns=["relative_position", "value", "group"])

    # Merge and calculate the mean of normalized values by relative position
    merged_df = pd.concat(normalized_dfs).groupby("relative_position").mean().reset_index()
    merged_df["value"] = merged_df["value"].fillna(0)

    # Group by bins of 5% intervals
    merged_df["bin"] = pd.cut(merged_df["relative_position"], bins=np.linspace(0, 100, 21), labels=np.arange(0, 100, 5))
    binned_df = merged_df.groupby("bin").mean().reset_index()

    # Replace zero-filled bins with interpolation only
    binned_df["value"] = binned_df["value"].replace(0, np.nan).interpolate()

    # Log bin counts beyond 30% to monitor coverage
    bin_counts = binned_df[binned_df["bin"] > 30]
    logging.info(f"Bin counts beyond 30%: {bin_counts.shape[0]}")

    return binned_df

# Multi-processing function for group-level averaging
def process_group(group):
    logging.info(f"Processing group: {group}")
    group_files = config["samples"][group]
    avg_coverage = calculate_average_coverage(group_files)
    avg_coverage["group"] = group
    return avg_coverage

# Plotting function with reduced smoothing
def plot_metagene(data, output_file):
    logging.info("Plotting metagene data")
    plt.figure(figsize=(12, 6))

    for group, group_data in data.items():
        group_data["smoothed_value"] = group_data["value"].rolling(window=5, min_periods=1).mean()
        plt.plot(group_data["bin"], group_data["smoothed_value"], label=group, lw=2)

    plt.axvline(x=0, color='red', linestyle='--', label="Stop Codon")
    plt.xlabel("Relative Position in 3'UTR (%) (Binned in 5% intervals)", fontsize=14)
    plt.ylabel("Normalized Coverage (RPM)", fontsize=14)
    plt.title("Normalized Metagene Plot for Stop Codon and 3'UTR (Smoothed)", fontsize=16)
    plt.xticks(np.arange(0, 101, 10))
    plt.legend()
    plt.grid(True)

    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    plt.savefig(output_file, dpi=300)
    plt.close()
    logging.info(f"Plot saved to {output_file}")

# Main script
if __name__ == "__main__":
    logging.info("Starting metagene plot generation")

    num_threads = cpu_count()
    logging.info(f"Using {num_threads} CPU threads for parallel processing")

    try:
        with Pool(processes=num_threads) as pool:
            results = pool.map(process_group, config["samples"].keys())

        combined_results = {result["group"].iloc[0]: result for result in results}

        output_file = os.path.join(config["Files"]["output_dir"], "final_restored_metagene_plot.png")
        plot_metagene(combined_results, output_file)

        logging.info("Metagene plot generation completed successfully")

    except Exception as e:
        logging.error(f"An error occurred: {e}")
        raise
