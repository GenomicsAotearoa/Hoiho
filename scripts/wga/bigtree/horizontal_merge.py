#!/usr/bin/env python3
import glob
import os
from Bio import SeqIO

def merge_fasta_files(file_pattern, merged_filename, partition_filename):
    # Find all FASTA files matching the pattern and sort them
    fasta_files = sorted(glob.glob(file_pattern))
    if not fasta_files:
        print("No files found with pattern:", file_pattern)
        return

    # Determine the union of all taxon names across files and store records for each file
    taxa_set = set()
    file_records = {}  # Map filename to its list of SeqIO records
    for fasta in fasta_files:
        records = list(SeqIO.parse(fasta, "fasta"))
        file_records[fasta] = records
        for record in records:
            taxa_set.add(record.id)
    taxa_list = sorted(taxa_set)

    # Initialize a dictionary to hold the concatenated sequences for each taxon
    concatenated_seqs = {taxon: "" for taxon in taxa_list}

    # List to hold partition info: each entry is (partition_name, start, end)
    partitions = []
    current_start = 1

    for fasta in fasta_files:
        records = file_records[fasta]
        if not records:
            continue
        # Create a dictionary for quick lookup: taxon -> sequence for the current file
        current_dict = {record.id: str(record.seq) for record in records}
        # Assume all sequences in the file have the same length
        aln_length = len(records[0].seq)
        current_end = current_start + aln_length - 1

        # Use the filename (without extension) as the partition name
        partition_name = os.path.splitext(os.path.basename(fasta))[0]
        partitions.append((partition_name, current_start, current_end))

        # For each taxon, append its sequence from the current file or gaps if missing
        for taxon in taxa_list:
            if taxon in current_dict:
                concatenated_seqs[taxon] += current_dict[taxon]
            else:
                concatenated_seqs[taxon] += '-' * aln_length

        current_start = current_end + 1

    # Write the merged alignment to a FASTA file
    with open(merged_filename, 'w') as out_f:
        for taxon in taxa_list:
            out_f.write(f">{taxon}\n{concatenated_seqs[taxon]}\n")

    # Write the partition file (format: DNA, partition_name = start-end)
    with open(partition_filename, 'w') as part_f:
        for partition_name, start, end in partitions:
            part_f.write(f"DNA, {partition_name} = {start}-{end}\n")

    print(f"Merged FASTA written to {merged_filename}")
    print(f"Partition file written to {partition_filename}")

if __name__ == "__main__":
    # Adjust the file pattern if needed
    file_pattern = "reheader_out-*.fasta"
    merged_filename = "merged.fasta"
    partition_filename = "partitions.txt"

    merge_fasta_files(file_pattern, merged_filename, partition_filename)
