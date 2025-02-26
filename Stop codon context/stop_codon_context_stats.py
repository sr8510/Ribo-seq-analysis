import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from collections import Counter, defaultdict
import scipy.stats as stats
from scipy.stats import zscore
from scipy.stats import fisher_exact, mannwhitneyu
import numpy as np

# Define isoacceptor groups
isoacceptors = {
    "TTT": "Phe", "TTC": "Phe",
    "TTA": "Leu", "TTG": "Leu", "CTT": "Leu", "CTC": "Leu", "CTA": "Leu", "CTG": "Leu",
    "ATT": "Ile", "ATC": "Ile", "ATA": "Ile",
    "ATG": "Met",
    "GTT": "Val", "GTC": "Val", "GTA": "Val", "GTG": "Val",
    "TCT": "Ser", "TCC": "Ser", "TCA": "Ser", "TCG": "Ser", "AGT": "Ser", "AGC": "Ser",
    "CCT": "Pro", "CCC": "Pro", "CCA": "Pro", "CCG": "Pro",
    "ACT": "Thr", "ACC": "Thr", "ACA": "Thr", "ACG": "Thr",
    "GCT": "Ala", "GCC": "Ala", "GCA": "Ala", "GCG": "Ala",
    "TAT": "Tyr", "TAC": "Tyr",
    "CAT": "His", "CAC": "His",
    "CAA": "Gln", "CAG": "Gln",
    "AAT": "Asn", "AAC": "Asn",
    "AAA": "Lys", "AAG": "Lys",
    "GAT": "Asp", "GAC": "Asp",
    "GAA": "Glu", "GAG": "Glu",
    "TGT": "Cys", "TGC": "Cys",
    "TGG": "Trp",
    "CGT": "Arg", "CGC": "Arg", "CGA": "Arg", "CGG": "Arg", "AGA": "Arg", "AGG": "Arg",
    "GGT": "Gly", "GGC": "Gly", "GGA": "Gly", "GGG": "Gly"
}

# Load the dataset
df = pd.read_csv("stop_codon_context.csv")

# Separate Upregulated and Downregulated Genes
upregulated = df[df["Direction"] == "Upregulated"]
downregulated = df[df["Direction"] == "Downregulated"]

# Define frames
codon_positions = ["-1", "-2", "-3", "-4", "-5"]

# Compute codon frequencies for each frame in both groups
up_codon_freq = {pos: Counter(upregulated[pos]) for pos in codon_positions}
down_codon_freq = {pos: Counter(downregulated[pos]) for pos in codon_positions}

# Convert to a single DataFrame format
up_codon_df = pd.DataFrame(up_codon_freq).fillna(0)
down_codon_df = pd.DataFrame(down_codon_freq).fillna(0)

# Normalize by total occurrences (percentage)
for pos in codon_positions:
    up_codon_df[pos] = up_codon_df[pos] / up_codon_df[pos].sum() * 100
    down_codon_df[pos] = down_codon_df[pos] / down_codon_df[pos].sum() * 100

# Save the processed codon frequency data
up_codon_df.to_csv("upregulated_codon_frequencies.csv")
down_codon_df.to_csv("downregulated_codon_frequencies.csv")

# Fisher's Exact Test for each codon at each position
fisher_results = []
for pos in codon_positions:
    for codon in set(up_codon_df.index).union(down_codon_df.index):
        up_count = up_codon_freq[pos].get(codon, 0)
        down_count = down_codon_freq[pos].get(codon, 0)
        total_up = sum(up_codon_freq[pos].values())
        total_down = sum(down_codon_freq[pos].values())

        # Construct contingency table
        table = [[up_count, total_up - up_count], [down_count, total_down - down_count]]

        # Perform Fisher's exact test
        _, p_value = stats.fisher_exact(table, alternative='two-sided')

        fisher_results.append([pos, codon, up_count, down_count, p_value])

# Convert Fisher test results to DataFrame
fisher_df = pd.DataFrame(fisher_results, columns=["Frame", "Codon", "Upregulated_Count", "Downregulated_Count", "p_value"])
fisher_df["Significant"] = fisher_df["p_value"] < 0.05
fisher_df.to_csv("fisher_codon_comparison.csv", index=False)

# Analyze isoacceptor frequencies per codon
up_isoacceptor_freq = defaultdict(lambda: {pos: 0 for pos in codon_positions})
down_isoacceptor_freq = defaultdict(lambda: {pos: 0 for pos in codon_positions})

for codon, amino_acid in isoacceptors.items():
    for pos in codon_positions:
        total_up = sum(up_codon_df[pos].get(iso_codon, 0) for iso_codon in isoacceptors if isoacceptors[iso_codon] == amino_acid)
        total_down = sum(down_codon_df[pos].get(iso_codon, 0) for iso_codon in isoacceptors if isoacceptors[iso_codon] == amino_acid)
        up_isoacceptor_freq[codon][pos] = up_codon_df[pos].get(codon, 0) / total_up if total_up > 0 else 0
        down_isoacceptor_freq[codon][pos] = down_codon_df[pos].get(codon, 0) / total_down if total_down > 0 else 0

up_isoacceptor_df = pd.DataFrame(up_isoacceptor_freq).T
down_isoacceptor_df = pd.DataFrame(down_isoacceptor_freq).T

up_isoacceptor_df.to_csv("upregulated_isoacceptor_frequencies.csv")
down_isoacceptor_df.to_csv("downregulated_isoacceptor_frequencies.csv")


# Statistical analysis: Fisher’s Exact Test
stats_results = []
for pos in codon_positions:
    frame_results = []
    for codon in set(up_codon_df.index).union(set(down_codon_df.index)):
        up_count = up_codon_df[pos].get(codon, 0)
        down_count = down_codon_df[pos].get(codon, 0)
        total_up = up_codon_df[pos].sum()
        total_down = down_codon_df[pos].sum()
        
        # Avoid invalid Fisher's test scenarios
        if total_up == 0 or total_down == 0 or (up_count == 0 and down_count == 0):
            p_value = 1  # Default p-value when one group has zero occurrences
        else:
            contingency_table = [[up_count, down_count],
                                 [total_up - up_count, total_down - down_count]]
            _, p_value = fisher_exact(contingency_table)
        
        frame_results.append({"Frame": pos, "Codon": codon, "p_value": p_value})
    
    stats_results.extend(frame_results)

stats_df = pd.DataFrame(stats_results)
stats_df["Significant"] = stats_df["p_value"] < 0.05
stats_df.to_csv("isoacceptors_statistics_per_codon.csv", index=False)

# Convert isoacceptor frequencies to z-scores for each frame
for pos in codon_positions:
    up_isoacceptor_df[pos] = stats.zscore(up_isoacceptor_df[pos].astype(float), nan_policy="omit")
    down_isoacceptor_df[pos] = stats.zscore(down_isoacceptor_df[pos].astype(float), nan_policy="omit")

# Save the z-score transformed isoacceptor frequency data
up_isoacceptor_df.to_csv("upregulated_isoacceptor_zscores.csv")
down_isoacceptor_df.to_csv("downregulated_isoacceptor_zscores.csv")

# Generate isoacceptor bar plots using z-scores for each frame separately
for pos in codon_positions:
    plt.figure(figsize=(15, 6))

    # Get the number of codons for spacing
    indices = np.arange(len(up_isoacceptor_df.index))
    bar_width = 0.4  # Adjust width to avoid overlap

    # Plot side-by-side bars
    plt.bar(indices - bar_width/2, up_isoacceptor_df[pos], bar_width, color="blue", label="Upregulated")
    plt.bar(indices + bar_width/2, down_isoacceptor_df[pos], bar_width, color="red", label="Downregulated")

    plt.ylabel("Isoacceptor Z-score")
    plt.xlabel("Amino Acid Isoacceptor")
    plt.title(f"Isoacceptor Frequency Comparison (Z-score) at Frame {pos}")
    plt.xticks(indices, up_isoacceptor_df.index, rotation=90)
    plt.legend()
    plt.tight_layout()
    
    # Save the corrected plot
    plt.savefig(f"isoacceptor_comparison_zscore_frame_{pos}.png")
    plt.show()

# Generate individual bar plots for each frame using codon percentages
for pos in codon_positions:
    combined_df = pd.concat([up_codon_df[pos], down_codon_df[pos]], axis=1).fillna(0)

    plt.figure(figsize=(12, 6))  # Enlarged figure
    ax = combined_df.plot(kind="bar", color=["blue", "red"], alpha=0.6)

    plt.ylabel("Percentage (%)", fontsize=14)  # Larger font
    plt.xlabel("Codon", fontsize=14)
    plt.title(f"Codon Frequency Comparison at Frame {pos}: Upregulated vs Downregulated Genes", fontsize=10)

    plt.xticks(rotation=90, fontsize=12)  # Improve readability
    plt.yticks(fontsize=12)
    plt.legend(["Upregulated", "Downregulated"], fontsize=12)

    plt.tight_layout()
    
    # Save each frame plot separately
    plt.savefig(f"codon_comparison_frame_{pos}.png", dpi=300)  # Higher resolution
    plt.show()

# heatmap for isoacceptors
# Prepare data for heatmap
heatmap_data = []
for pos in codon_positions:
    for regulation in ["Upregulated", "Downregulated"]:
        if regulation == "Upregulated":
            df = up_isoacceptor_df
        else:
            df = down_isoacceptor_df
        
        for codon in df.index:
            heatmap_data.append([codon, f"{pos} ({regulation})", df.at[codon, pos]])

# Convert to DataFrame
heatmap_df = pd.DataFrame(heatmap_data, columns=["Codon", "Frame", "Z-score"])
heatmap_pivot = heatmap_df.pivot(index="Codon", columns="Frame", values="Z-score")

# Save the heatmap data as a CSV file
heatmap_pivot.to_csv("isoacceptor_heatmap_data.csv")

# Generate heatmap
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_pivot, cmap="coolwarm", center=0, annot=False, fmt=".1f", linewidths=0.5)

plt.xlabel("Frames (Upregulated & Downregulated)")
plt.ylabel("Isoacceptor Codon")
plt.title("Isoacceptor Frequency Z-score Heatmap")
plt.xticks(rotation=45, ha="right")
plt.yticks(rotation=0)
plt.tight_layout()

# Save the heatmap
plt.savefig("isoacceptor_heatmap.png")
plt.show()

# Print significant Fisher test results
print("Significant Codon Differences:")
print(fisher_df[fisher_df["Significant"]])
