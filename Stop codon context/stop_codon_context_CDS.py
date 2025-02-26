import pandas as pd
from Bio import SeqIO, SeqRecord, Seq
import argparse
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

# Argument parser
def parse_arguments():
    parser = argparse.ArgumentParser(description="Analyze stop codon context using CDS FASTA sequences.")
    parser.add_argument("-i", "--input", required=True, help="Path to the input CSV file.")
    parser.add_argument("-f", "--fasta", required=True, help="Path to the CDS FASTA file.")
    parser.add_argument("-o", "--output", required=True, help="Output directory for the results.")
    parser.add_argument("-t", "--threads", type=int, default=4, help="Number of threads to use (default: 4).")
    return parser.parse_args()

# Function to get upstream codons
def get_upstream_codons(sequence, stop_codon_pos, num_codons=5):
    codons = []
    for i in range(1, num_codons + 1):
        start = stop_codon_pos - (i * 3)
        if start < 0:
            codons.append("NNN")
        else:
            codons.append(sequence[start:start + 3])
    return codons[::-1]  # Reverse to maintain -1 to -5 order

def identify_stop_codon_and_context(cds, num_codons=5):
    cds = cds.strip().replace("\n", "")  # Ensure no extra whitespace or newlines

    if len(cds) < 3:
        return "NNN", ["NNN"] * num_codons  # Return NNN if CDS is too short

    # Iterate backwards in steps of 3 to find a valid stop codon
    for i in range(len(cds) - 3, -1, -3):
        stop_codon = cds[i:i + 3]
        if stop_codon in {"TAA", "TAG", "TGA"}:
            stop_codon_pos = i
            upstream_codons = get_upstream_codons(cds, stop_codon_pos, num_codons)
            return stop_codon, upstream_codons

    # If no valid stop codon is found, return NNN
    return "NNN", ["NNN"] * num_codons

# Function to process each FASTA record
def process_fasta_record(record, significant_genes):
    # Extract gene name from 'gene_symbol:' in the description
    match = re.search(r"gene_symbol:([^\s]+)", record.description)
    if match:
        gene_name = match.group(1).lower()  # Convert to lowercase for matching
        significant_genes['Gene'] = significant_genes['Gene'].str.lower()
        if gene_name in significant_genes['Gene'].values:
            cds = str(record.seq)
            return gene_name, cds
    return None

def select_longest_transcripts(cds_dict):
    """Select the longest transcript for each gene."""
    longest_transcripts = {}
    for gene, transcripts in cds_dict.items():
        longest_transcript = max(transcripts, key=lambda x: len(x[1]))  # x[1] is the CDS sequence
        longest_transcripts[gene] = longest_transcript[1]  # Store the longest CDS
    return longest_transcripts

def main():
    args = parse_arguments()
    input_csv = args.input
    fasta_file = args.fasta
    output_dir = args.output
    num_threads = args.threads
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Read the input CSV file
    data = pd.read_csv(input_csv)
    significant_genes = data[(data['p_value'] < 0.05) & (abs(data['log2FC']) > 1)].copy()
    significant_genes['Direction'] = significant_genes['log2FC'].apply(lambda x: 'Upregulated' if x > 0 else 'Downregulated')

    # Parse the FASTA file and collect all transcripts for each gene using multithreading
    cds_dict = {}
    with ThreadPoolExecutor(max_workers=num_threads) as executor:
        futures = [executor.submit(process_fasta_record, record, significant_genes) for record in SeqIO.parse(fasta_file, "fasta")]
        for future in as_completed(futures):
            result = future.result()
            if result:
                gene_name, cds = result
                if gene_name not in cds_dict:
                    cds_dict[gene_name] = []
                cds_dict[gene_name].append((gene_name, cds))

    # Select the longest transcript for each gene
    longest_transcripts = select_longest_transcripts(cds_dict)

    # Output FASTA files and stop codon context for longest transcripts
    upregulated_records = []
    downregulated_records = []
    output_data = []

    for gene, cds in longest_transcripts.items():
        stop_codon, upstream_codons = identify_stop_codon_and_context(cds)
        direction = significant_genes[significant_genes['Gene'].str.lower() == gene]['Direction'].values[0]
        record = SeqRecord.SeqRecord(Seq.Seq(cds), id=gene, description=f"{gene} {direction}")
        
        if direction == 'Upregulated':
            upregulated_records.append(record)
        else:
            downregulated_records.append(record)
        
        output_data.append([gene, direction, stop_codon] + upstream_codons)

    # Write FASTA files
    SeqIO.write(upregulated_records, os.path.join(output_dir, "upregulated_genes.fasta"), "fasta")
    SeqIO.write(downregulated_records, os.path.join(output_dir, "downregulated_genes.fasta"), "fasta")

    # Save the stop codon context to CSV
    stop_codon_context_path = os.path.join(output_dir, "stop_codon_context.csv")
    output_df = pd.DataFrame(output_data, columns=["Gene", "Direction", "Stop codon", "-1", "-2", "-3", "-4", "-5"])
    output_df.to_csv(stop_codon_context_path, index=False)

    print(f"Analysis complete. Results saved to {stop_codon_context_path}.")


if __name__ == "__main__":
    main()
