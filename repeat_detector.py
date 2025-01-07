# detect the types and numbers of repeats in a silk protein

import sys
import csv
from Bio import SeqIO
from collections import defaultdict

# dipeptide repeat
def count_non_overlapping_repeats(sequence, dipeptide):
    """Count non-overlapping occurrences of a specific dipeptide in the sequence."""
    count = 0
    i = 0
    while i <= len(sequence) - 2:
        if sequence[i:i + 2] == dipeptide:
            count += 1
            i += 2  # Move by 2 to prevent overlapping
        else:
            i += 1
    return count

def find_all_gly_dipeptide_repeats(sequence):
    """Find all glycine-containing dipeptides and their non-overlapping counts."""
    gly_dipeptide_counts = defaultdict(int)
    protein_length = len(sequence)

    # Generate all unique dipeptides with Glycine (G) in at least one position.
    for i in range(protein_length - 1):
        dipeptide = sequence[i:i + 2]
        if 'G' in dipeptide:
            gly_dipeptide_counts[dipeptide] = count_non_overlapping_repeats(sequence, dipeptide)

    return gly_dipeptide_counts, protein_length

def process_fasta(input_file, output_file):
    """Process the input FASTA file and write results to a CSV file."""
    all_gly_repeats_set = set()
    gene_repeat_info = {}

    for record in SeqIO.parse(input_file, "fasta"):
        gene_name = record.id
        sequence = str(record.seq)
        gly_dipeptide_counts, protein_length = find_all_gly_dipeptide_repeats(sequence)
        gene_repeat_info[gene_name] = (gly_dipeptide_counts, protein_length)
        all_gly_repeats_set.update(gly_dipeptide_counts.keys())

    all_gly_repeats = sorted(list(all_gly_repeats_set))

    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        header = ["gene name", "protein length"] + all_gly_repeats + ["most frequent repeat", "repeat count"]
        csv_writer.writerow(header)

        for gene_name, (gly_dipeptide_counts, protein_length) in gene_repeat_info.items():
            # Determine the most frequent dipeptide and its count
            if gly_dipeptide_counts:
                most_frequent_repeat = max(gly_dipeptide_counts, key=gly_dipeptide_counts.get)
                max_count = gly_dipeptide_counts[most_frequent_repeat]
            else:
                most_frequent_repeat = "None"
                max_count = 0

            row = [gene_name, protein_length] + [gly_dipeptide_counts.get(dipeptide, 0) for dipeptide in all_gly_repeats]
            row += [most_frequent_repeat, max_count]
            csv_writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python repeat_detector.py input.fas dipeptide_repeat.csv")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_csv = sys.argv[2]
    process_fasta(input_fasta, output_csv)

# tripeptide repeat
def count_non_overlapping_repeats(sequence, tripeptide):
    """Count non-overlapping occurrences of a specific tripeptide in the sequence."""
    count = 0
    i = 0
    while i <= len(sequence) - 3:
        if sequence[i:i + 3] == tripeptide:
            count += 1
            i += 3  # Move by 3 to prevent overlapping
        else:
            i += 1
    return count

def find_all_gly_tripeptide_repeats(sequence):
    """Find all glycine-containing tripeptides and their non-overlapping counts."""
    gly_tripeptide_counts = defaultdict(int)
    protein_length = len(sequence)

    # Generate all unique tripeptides with Glycine (G) in at least one position.
    for i in range(protein_length - 2):
        tripeptide = sequence[i:i + 3]
        if 'G' in tripeptide:
            gly_tripeptide_counts[tripeptide] = count_non_overlapping_repeats(sequence, tripeptide)

    return gly_tripeptide_counts, protein_length

def process_fasta(input_file, output_file):
    """Process the input FASTA file and write results to a CSV file."""
    all_gly_repeats_set = set()
    gene_repeat_info = {}

    for record in SeqIO.parse(input_file, "fasta"):
        gene_name = record.id
        sequence = str(record.seq)
        gly_tripeptide_counts, protein_length = find_all_gly_tripeptide_repeats(sequence)
        gene_repeat_info[gene_name] = (gly_tripeptide_counts, protein_length)
        all_gly_repeats_set.update(gly_tripeptide_counts.keys())

    all_gly_repeats = sorted(list(all_gly_repeats_set))

    with open(output_file, 'w', newline='') as csvfile:
        csv_writer = csv.writer(csvfile)
        header = ["gene name", "protein length"] + all_gly_repeats + ["most frequent repeat", "repeat count"]
        csv_writer.writerow(header)

        for gene_name, (gly_tripeptide_counts, protein_length) in gene_repeat_info.items():
            # Determine the most frequent tripeptide and its count
            if gly_tripeptide_counts:
                most_frequent_repeat = max(gly_tripeptide_counts, key=gly_tripeptide_counts.get)
                max_count = gly_tripeptide_counts[most_frequent_repeat]
            else:
                most_frequent_repeat = "None"
                max_count = 0

            row = [gene_name, protein_length] + [gly_tripeptide_counts.get(tripeptide, 0) for tripeptide in all_gly_repeats]
            row += [most_frequent_repeat, max_count]
            csv_writer.writerow(row)

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python repeat_detector.py input.fa tripeptide_repeat.csv")
        sys.exit(1)

    input_fasta = sys.argv[1]
    output_csv = sys.argv[2]
    process_fasta(input_fasta, output_csv)
