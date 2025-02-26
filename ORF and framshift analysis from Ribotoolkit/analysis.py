#!/usr/bin/env python3

import os
import argparse
import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.multitest import multipletests

def main():
    parser = argparse.ArgumentParser(description="Analyze differentially expressed ORFs.")
    parser.add_argument("-s", "--stats_file", required=True, help="Path to Limma output CSV.")
    parser.add_argument("-o", "--output_prefix", required=True, help="Output file prefix (including directory).")
    args = parser.parse_args()
    
    # 1. Load Limma stats data
    df = pd.read_csv(args.stats_file)
    
    # Some Limma output files have "X" or unnamed columns at the start; adjust if needed.
    # ensure we have columns of interest:
    needed_cols = {"ORF_ID", "ORF_type", "transcript_type", "AAseq", "logFC", "adj.P.Val"}
    # If the first column is named "X", sometimes "ORF_ID" is in that column. 
    # Adjust if your file's first column is actually the ORF IDs.
    if "ORF_ID" not in df.columns and "X" in df.columns:
        df.rename(columns={"X": "ORF_ID"}, inplace=True)
    
    missing = needed_cols - set(df.columns)
    if missing:
        raise ValueError(f"Missing required columns in {args.stats_file}: {missing}")
    
    # 2. Filter for significant differentially expressed ORFs
    sig_mask = (df["adj.P.Val"] < 0.05) & (df["logFC"].abs() > 1)
    df_filt = df[sig_mask].copy()
    
    # 3. Summaries by ORF_type and transcript_type
    orf_type_counts = df_filt["ORF_type"].value_counts().reset_index()
    orf_type_counts.columns = ["ORF_type", "count"]
    
    transcript_type_counts = df_filt["transcript_type"].value_counts().reset_index()
    transcript_type_counts.columns = ["transcript_type", "count"]
    
    # Save these
    out_dir = os.path.dirname(args.output_prefix)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir)
    
    orf_type_counts.to_csv(f"{args.output_prefix}_ORFtype_summary.csv", index=False)
    transcript_type_counts.to_csv(f"{args.output_prefix}_TranscriptType_summary.csv", index=False)
    
    # 4. Prepare for amino-acid analysis
    #    Ignore non-standard chars (x or *), skip ORFs with no valid AAseq.
    standard_aas = list("ACDEFGHIKLMNPQRSTVWY")  # 20 standard residues
    
    def clean_sequence(seq):
        # Remove any lower/upper X or '*' or whitespace
        seq = seq.upper().replace("*", "").replace("X", "")
        return seq.strip()
    
    df_filt["AAseq"] = df_filt["AAseq"].astype(str).apply(clean_sequence)
    df_filt = df_filt[df_filt["AAseq"].str.len() > 0]  # Drop if sequence is empty
    
    # 5. Count each amino acid and compute z-scores within each sequence
    zscore_rows = []
    
    for i, row in df_filt.iterrows():
        orf_id = row["ORF_ID"]
        seq = row["AAseq"]
        logfc = row["logFC"]
        group = "Upregulated" if logfc > 1 else "Downregulated"
        
        # Count each standard amino acid in this sequence
        counts = [seq.count(aa) for aa in standard_aas]
        # Compute mean & std dev among the 20 amino acids
        mean_aa = np.mean(counts)
        std_aa  = np.std(counts, ddof=1)  # sample std dev
        
        if std_aa == 0:
            # If they are all zero or all identical, z-scores would be inf or 0
            # We can skip or keep them as 0. We'll keep them as 0 in this example.
            zscores = [0]*len(standard_aas)
        else:
            zscores = [(c - mean_aa)/std_aa for c in counts]
        
        # Build row for output
        row_dict = {
            "ORF_ID": orf_id,
            "Group": group
        }
        # Add each amino acid's z-score
        for aa, z in zip(standard_aas, zscores):
            row_dict[aa] = z
        
        zscore_rows.append(row_dict)
    
    # Create z-score matrix
    zscore_df = pd.DataFrame(zscore_rows, columns=["ORF_ID","Group"]+standard_aas)
    
    # Save the z-score matrix
    zscore_df.to_csv(f"{args.output_prefix}_AAzscore_matrix.csv", index=False)
    
    # 6. Mann-Whitney U test: compare up vs down for each AA's z-score
    #    We'll do BH correction across the 20 tests.
    up_mask = zscore_df["Group"] == "Upregulated"
    down_mask = zscore_df["Group"] == "Downregulated"
    
    stats_list = []
    
    for aa in standard_aas:
        up_values = zscore_df.loc[up_mask, aa].values
        down_values = zscore_df.loc[down_mask, aa].values
        
        # Mean z-scores in each group
        up_mean = up_values.mean() if len(up_values) > 0 else float("nan")
        down_mean = down_values.mean() if len(down_values) > 0 else float("nan")

        # Only run test if we have at least 2 values in each group
        if len(up_values) < 2 or len(down_values) < 2:
            u_stat, pval = (np.nan, np.nan)
        else:
            u_stat, pval = st.mannwhitneyu(up_values, down_values, alternative="two-sided")
        
        stats_list.append([aa, up_mean, down_mean, u_stat, pval])
    
    stats_df = pd.DataFrame(stats_list, columns=["AminoAcid", "UpMeanZ", "DownMeanZ", "U_stat", "p_value"])
    
    # multiple testing correction
    valid_mask = ~stats_df["p_value"].isna()
    corrections = multipletests(stats_df.loc[valid_mask, "p_value"], method="fdr_bh")
    stats_df.loc[valid_mask, "adj_p_value"] = corrections[1]
    
    # Save test results
    stats_df.to_csv(f"{args.output_prefix}_AAs_MannWhitney.csv", index=False)
    
    print("Analysis complete.")
    print(f"Filtered ORFs: {df_filt.shape[0]}")
    print(f"Outputs:\n  {args.output_prefix}_ORFtype_summary.csv\n  {args.output_prefix}_TranscriptType_summary.csv\n  {args.output_prefix}_AAzscore_matrix.csv\n  {args.output_prefix}_AAs_MannWhitney.csv")

if __name__ == "__main__":
    main()
