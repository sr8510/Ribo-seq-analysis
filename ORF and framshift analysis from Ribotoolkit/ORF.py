import os
import yaml
import pandas as pd
import numpy as np
import logging

# Setup logging
log_file = "script_debug.log"
logging.basicConfig(filename=log_file, level=logging.DEBUG, format="%(asctime)s - %(levelname)s - %(message)s")

logging.info("Script started.")

# Load configuration file
with open("Config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Output directories
output_dir = "results"
active_translation_dir = os.path.join(output_dir, "active_translation")
frameshift_analysis_dir = os.path.join(output_dir, "frameshift_analysis")

for directory in [
    active_translation_dir, 
    os.path.join(active_translation_dir, "normalized_data"),
    os.path.join(active_translation_dir, "differential_analysis"),
    os.path.join(active_translation_dir, "volcano_plots"),
    frameshift_analysis_dir, 
    os.path.join(frameshift_analysis_dir, "normalized_data"),
    os.path.join(frameshift_analysis_dir, "differential_analysis"),
    os.path.join(frameshift_analysis_dir, "volcano_plots"),
]:
    os.makedirs(directory, exist_ok=True)

# Extract annotation data
annotation_dict = {}

def process_file(file_path):
    df = pd.read_csv(file_path, sep=",")

    # Deduplicate rows that have both the same ORF_ID and the same AAseq
    df.drop_duplicates(subset=["ORF_ID", "AAseq"], keep="first", inplace=True)

    # Now store annotation data after the deduplication
    for _, row in df.iterrows():
        orf_id = row["ORF_ID"]
        annotation_dict[orf_id] = row[[
            "transcript_id", "tx_len", "utr5_len", "cds_len", "utr3_len",
            "ORF_type", "transcript_type", "gene_id", "gene_name", "gene_type",
            "chrom", "strand", "ORF_length", "ORF_tstart", "ORF_tstop",
            "ORF_gstart", "ORF_gstop", "AAseq"
        ]].to_dict()

    return df

# Process all files
sample_groups = config["samples"]
control_group = config["control_condition"]

all_samples = {group: [process_file(file) for file in files] for group, files in sample_groups.items()}

# Function to normalize Psites_sum_cov to CPM
def normalize_to_cpm(df):
    df["Psites_CPM"] = (df["Psites_sum_cov"] / df["Psites_sum_cov"].sum()) * 1e6
    return df

# Normalize all samples
normalized_data = {
    group: [normalize_to_cpm(df) for df in files]
    for group, files in all_samples.items()
}

# Save normalized data for active translation analysis
for group, files in normalized_data.items():
    for i, df in enumerate(files):
        df.to_csv(os.path.join(active_translation_dir, "normalized_data", f"{group}_{i+1}.csv"), index=False)

# Merge normalized CPM values into one matrix
# Collect all normalized samples in a list.
long_cpm_list = []

for group, files in normalized_data.items():
    for i, df in enumerate(files):
        sample_name = f"{group}_{i+1}"  # e.g. "control_1", "1hr_2", etc.

        # Create a temp copy that has 3 columns: ORF_ID, value, sample
        tmp = df[["ORF_ID", "Psites_CPM"]].copy()
        tmp.columns = ["ORF_ID", "value"]
        tmp["sample"] = sample_name

        long_cpm_list.append(tmp)

# Now concatenate everything into one "long" DataFrame
long_cpm_df = pd.concat(long_cpm_list, ignore_index=True)

long_cpm_df = long_cpm_df.drop_duplicates(subset=["ORF_ID", "sample"])

# Pivot from "long" to "wide" format, with one column per sample
merged_cpm = long_cpm_df.pivot_table(
    index="ORF_ID",   # Just use ORF_ID
    columns="sample",
    values="value",
    aggfunc="sum",
    fill_value=0
)

# Fill missing values with 0
merged_cpm = merged_cpm.fillna(0).reset_index()

# Save merged CPM matrix
merged_cpm.to_csv(os.path.join(active_translation_dir, "merged_counts.csv"), index=False)

print("Merged CPM counts saved.")

# Compute Frameshift Ratio (FR)
def compute_frameshift_ratio(df):
    df["FR"] = (df["Psites_sum_frame1"] + df["Psites_sum_frame2"]) / df["Psites_sum_frame0"].replace(0, np.nan)
    return df

# Apply FR computation
fr_data = {group: [compute_frameshift_ratio(df) for df in files] for group, files in all_samples.items()}

# Save FR data for frameshift analysis
for group, files in fr_data.items():
    for i, df in enumerate(files):
        df.to_csv(os.path.join(frameshift_analysis_dir, "normalized_data", f"{group}_{i+1}.csv"), index=False)

# Collect FR data into a long DataFrame for pivot
long_FR_list = []

# ***IMPORTANT: Use fr_data here, not normalized_data***
for group, files in fr_data.items():
    for i, df in enumerate(files):
        sample_name = f"{group}_{i+1}"

        # Either rename "FR" -> "value" so pivot uses values="value",
        # or pivot with values="FR". Let's rename to "value":
        tmp = df[["ORF_ID", "FR"]].copy()
        tmp.columns = ["ORF_ID", "value"]
        tmp["sample"] = sample_name

        long_FR_list.append(tmp)

long_FR_df = pd.concat(long_FR_list, ignore_index=True)

long_FR_df = long_FR_df.drop_duplicates(subset=["ORF_ID", "sample"])

merged_FR = long_FR_df.pivot(
    index="ORF_ID", 
    columns="sample", 
    values="value"
).fillna(0).reset_index()

merged_FR = merged_FR.fillna(0).reset_index()
merged_FR.to_csv(os.path.join(frameshift_analysis_dir, "merged_FR.csv"), index=False)
print("Merged Frameshift Ratios saved.")

# Save annotation dictionary
annotation_df = pd.DataFrame.from_dict(annotation_dict, orient="index")
annotation_df.to_csv(os.path.join(output_dir, "annotation_data.csv"))

print("Preprocessing complete. Run the R script for Limma analysis.")
