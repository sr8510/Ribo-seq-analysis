from Bio import SeqIO
from collections import defaultdict
import sys
from Bio.Seq import Seq

def parse_gene_list(gene_list_file):
    """Reads the list of genes to extract."""
    with open(gene_list_file, 'r') as f:
        return set(line.strip() for line in f)  # Keep gene names as provided

def parse_transcriptome_fasta(fasta_file):
    """Parses the transcriptome FASTA file and returns a dictionary of gene transcripts."""
    gene_transcripts = defaultdict(list)
    
    for record in SeqIO.parse(fasta_file, "fasta"):
        header = record.description
        seq = str(record.seq)  # Keep sequence case as provided
        header_parts = header.split()
        
        # Extract gene name from the first parentheses in the header
        gene_name = None
        for part in header_parts:
            if "(" in part and ")" in part:
                gene_name = part.strip("(),")
                break
        
        if not gene_name and len(header_parts) > 2:
            gene_name = header_parts[-2].strip("(),")
        elif not gene_name:
            gene_name = header_parts[1]
        
        gene_transcripts[gene_name].append((header, seq))
    
    return gene_transcripts

def get_longest_transcript(transcripts):
    """Returns the longest transcript for a given gene."""
    return max(transcripts, key=lambda x: len(x[1]))

def extract_cds(sequence):
    """Extracts the CDS sequence by identifying the first ATG and stopping at the first stop codon."""
    start_idx = sequence.find("ATG")
    if start_idx == -1:
        return None, None  # No valid start codon
    
    stop_codons = {"TAA", "TAG", "TGA"}
    cds_seq = ""
    for i in range(start_idx, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        cds_seq += codon
        if codon in stop_codons:
            return cds_seq, i + 3  # Return sequence and stop position
    
    return None, None  # No stop codon found

def extract_extended_cds(sequence, stop_pos):
    """Extracts the extended CDS by continuing past the first stop codon until a second stop codon is found."""
    stop_codons = {"TAA", "TAG", "TGA"}
    extended_seq = ""
    
    for i in range(stop_pos, len(sequence) - 2, 3):
        codon = sequence[i:i+3]
        if codon in stop_codons:
            return extended_seq  # Return extended CDS sequence without first stop codon
        extended_seq += codon
    
    return None  # No second stop codon found

def translate_sequence(dna_seq):
    """Translates a DNA sequence to an amino acid sequence."""
    return str(Seq(dna_seq).translate(to_stop=False))

def process_sequences(fasta_file, gene_list_file):
    """Processes the sequences and writes CDS and extended CDS to output files."""
    gene_list = parse_gene_list(gene_list_file)
    transcripts = parse_transcriptome_fasta(fasta_file)
    
    cds_records = []
    extended_cds_records = []
    cds_protein_records = []
    extended_cds_protein_records = []
    
    for gene in gene_list:
        if gene in transcripts:
            header, sequence = get_longest_transcript(transcripts[gene])
            cds_seq, stop_pos = extract_cds(sequence)
            
            if cds_seq:
                cds_records.append(f">{header}\n{cds_seq}")
                cds_protein_records.append(f">{header}\n{translate_sequence(cds_seq)}")
                extended_seq = extract_extended_cds(sequence, stop_pos)
                
                if extended_seq:
                    extended_cds = cds_seq + extended_seq
                    extended_cds_records.append(f">{header}\n{extended_cds}")
                    extended_cds_protein_records.append(f">{header}\n{translate_sequence(extended_cds)}")
            else:
                print(f"Warning: No valid CDS found for gene {gene}")
        else:
            print(f"Warning: Gene {gene} not found in FASTA file")
    
    # Write to FASTA files
    if cds_records:
        with open("CDS_sequences.fasta", "w") as f:
            f.write("\n".join(cds_records) + "\n")
        with open("CDS_proteins.fasta", "w") as f:
            f.write("\n".join(cds_protein_records) + "\n")
    else:
        print("No CDS sequences found.")
    
    if extended_cds_records:
        with open("Extended_CDS_sequences.fasta", "w") as f:
            f.write("\n".join(extended_cds_records) + "\n")
        with open("Extended_CDS_proteins.fasta", "w") as f:
            f.write("\n".join(extended_cds_protein_records) + "\n")
    else:
        print("No extended CDS sequences found.")
    
    print("Processing complete.")

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python extract_cds.py <transcriptome_fasta> <gene_list_txt>")
        sys.exit(1)
    
    fasta_file = sys.argv[1]
    gene_list_file = sys.argv[2]
    process_sequences(fasta_file, gene_list_file)