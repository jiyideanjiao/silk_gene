# amino acid composition

import sys
import csv
from Bio import SeqIO

def calculate_aa_percentage(sequence):
    """Calculate the percentage of each amino acid in the given sequence."""
    aa_counts = {aa: 0 for aa in "ARNDCQEGHILKMFPSTWYV"}
    total_length = len(sequence)
    
    for aa in sequence:
        if aa.upper() in aa_counts:
            aa_counts[aa.upper()] += 1
    
    aa_percentages = {aa: (count / total_length) * 100 for aa, count in aa_counts.items()}
    return aa_percentages

def process_fasta(input_file, output_file):
    """Process the input FASTA file and write results to a CSV file."""
    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        header = ["gene name", "protein length"] + list("ARNDCQEGHILKMFPSTWYV")
        csv_writer.writerow(header)

        for record in SeqIO.parse(input_file, "fasta"):
            gene_name = record.id
            sequence = str(record.seq)
            protein_length = len(sequence)
            aa_percentages = calculate_aa_percentage(sequence)

            row = [gene_name, protein_length] + [aa_percentages[aa] for aa in "ARNDCQEGHILKMFPSTWYV"]
            csv_writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python AA_info.py input.fa output.csv")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_csv = sys.argv[2]
    process_fasta(input_fasta, output_csv)
