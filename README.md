"""
# Ribo-seq-analysis
‘’’
How to generate annotation files used in the analysis:
(For this analysis, we used gtf downloaded from UCSC, CDS fasta downloaded from Ensembl, and transcriptome fasta downloaded from NCBI Refseq)

First, prepare the annotation file for the generation of bed files:

gffread annotation.gtf -T -o annotation.gtf

gtf2bed --gtf annotation.gtf --bed annotation.bed

#generation of CDS annotation

awk '$3 == "CDS" {
    match($0, /gene_name "([^"]+)"/, gname);
    if (gname[1] != "") print $1"\t"$4-1"\t"$5"\t"gname[1]"\t.\t"$7
}' annotation.gtf > annotation_CDS.bed

#generation for CDS annotation with exclusion of the first 18 and last 15 nucleotides:

awk '$3 == "CDS" {
    print $1 "\t" ($4-1) "\t" ($4+17) "\t.\t.\t" $7;  # First 18 nt
    print $1 "\t" ($5-15) "\t" $5 "\t.\t.\t" $7;      # Last 15 nt
}' annotation.gtf > exclude_regions.bed

bedtools subtract -a annotation_CDS.bed -b exclude_regions.bed > filtered_CDS.bed

#generation of 3'UTR annotation non-overlaping with CDS:

awk -F'\t' '
$3 == "CDS" {
    match($9, /gene_id "([^"]+)"/, gid);
    match($9, /transcript_id "([^"]+)"/, tid);
    gene_id = (gid[1] != "") ? gid[1] : "NA";
    transcript_id = (tid[1] != "") ? tid[1] : "NA";
    transcript_key = gene_id ";" transcript_id;
    if (!last_CDS_end[transcript_key] || $5 > last_CDS_end[transcript_key]) {
        last_CDS_end[transcript_key] = $5;
        last_CDS_chr[transcript_key] = $1;
        last_CDS_strand[transcript_key] = $7;
    }
}
END {
    for (transcript_key in last_CDS_end) {
        chr = last_CDS_chr[transcript_key];
        start = last_CDS_end[transcript_key];
        end = (last_CDS_strand[transcript_key] == "+") ? start + 100 : start - 100;
        if (last_CDS_strand[transcript_key] == "+") {
            print chr "\t" start "\t" end "\t" transcript_key "\t.\t+";
        } else {
            print chr "\t" end "\t" start "\t" transcript_key "\t.\t-";
        }
    }
}' annotation.gtf > initial_3UTR.bed

bedtools subtract -a initial_3UTR.bed -b annotation_CDS.bed > 3UTR.bed

bedtools intersect -a 3UTR.bed -b filtered_CDS.bed -wa -wb > overlaps_check.txt #check for overlapping

#generation of stop codon annotation:

awk -F'\t' '
$3 == "CDS" {
    match($9, /gene_id "([^"]+)"/, gid);
    match($9, /transcript_id "([^"]+)"/, tid);
    gene_id = (gid[1] != "") ? gid[1] : "NA";
    transcript_id = (tid[1] != "") ? tid[1] : "NA";
    transcript_key = gene_id ";" transcript_id;
    if (!last_CDS_end[transcript_key] || $5 > last_CDS_end[transcript_key]) {
        last_CDS_end[transcript_key] = $5;
        last_CDS_chr[transcript_key] = $1;
        last_CDS_strand[transcript_key] = $7;
    }
}
END {
    for (transcript_key in last_CDS_end) {
        print last_CDS_chr[transcript_key] "\t" (last_CDS_end[transcript_key]-3) "\t" last_CDS_end[transcript_key] "\t" transcript_key "\t.\t" last_CDS_strand[transcript_key];
    }
}' annotation.gtf > stop_codons.bed

#cleaning of Ensembl CDS fasta file:

grep -v "^>" Mus_musculus.GRCm39.cds.all.fa | grep -E "[^ACGTacgt]" #check for strange characters

#If NN appears then there is an issue

#How to clean:

awk '/^>/{print (NR==1 ? "" : "\n") $0; next}{printf "%s", $0}END{print "\n"}' Mus_musculus.GRCm39.cds.all.fa  > out.fasta

grep ">" out.fasta

sed 's/N/A/g' out.fasta > out_clean.fasta #use this file.

**How to run scripts:**

**SCRT:**
python RRS.py

**Stop codon context:**
python stop_codon_context_CDS.py -i “SCRT statistics file from SCRT analysis” -f “CDS fasta” -t “Number of threads” -o “output directory

**Stop codon context statistical analysis:**
python stop_codon_context_stats.py

**Metagene around 3’UTR:**
python Metagene-3UTR.py

**Elongation analysis:**
python elongation_analysis.py config config.yaml

**Metagene for specific genes CDS (linked to elongation analysis):**
python metagene.py

**Protein prediction for SCRT events:**
python extract_cds.py transcriptome.fasta gene_list.txt

**Gene-level frameshift analysis:**
python frameshift.py Config.yaml 

**Active ORF analysis (using Ribotoolkit input):**
python ORF.py #this script uses the translatedORF csv file output of the single sample analysis of Ribotoolkit as input. 

python analysis.py -s limma.csv -o outdir/name #for analyzing amino acids sequences and types of ORFs. This script is used on the Limma output. Make sure to use the annotation output file from the ORF.py script to annotate the Limma output so the necessary info are included in the final file


"""


