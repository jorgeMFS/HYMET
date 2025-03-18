#!/usr/bin/env python3

"""
process_sequences.py

This script processes sequencing data by converting a FASTA file into two paired-end FASTQ files with properly formatted headers and actual quality scores.

Usage:
    python process_sequences.py <input_file> [options]

Example:
    python process_sequences.py random_reference.fasta --read_length 300 --overlap 50 --output_l fastq_l.fastq --output_r fastq_r.fastq

Additional Options:
    --include_read_id    Include read identifier in the "+" separator line (default: False)
"""

import argparse
import os
import sys
import random
from Bio import SeqIO
from Bio.Seq import Seq

def is_fasta(file_path):
    """
    Determines if a file is in FASTA format based on its content.

    Parameters:
        file_path (str): Path to the input file.

    Returns:
        bool: True if FASTA, False otherwise.
    """
    try:
        with open(file_path, 'r') as f:
            first_char = f.read(1)
            return first_char == '>'
    except Exception as e:
        print(f"Error reading file {file_path}: {e}", file=sys.stderr)
        sys.exit(1)

def generate_quality_scores(length):
    """
    Generates a random quality score string of a given length with high-quality scores.

    Parameters:
        length (int): The length of the sequence.

    Returns:
        str: A string representing quality scores.
    """
    # Generate quality scores between Phred 30 and 40
    phred_scores = [random.randint(30, 40) for _ in range(length)]
    quality_scores = ''.join([chr(score + 33) for score in phred_scores])
    return quality_scores

def fasta_to_paired_fastq(fasta_file, read_length=150, overlap=50, output_l='fastq_l.fastq', output_r='fastq_r.fastq', include_read_id=False):
    """
    Converts a FASTA file to paired-end FASTQ files with properly formatted headers and overlapping reads.

    Parameters:
        fasta_file (str): Path to the input FASTA file.
        read_length (int): Length of each read.
        overlap (int): Number of bases overlapping between forward and reverse reads.
        output_l (str): Output FASTQ file for forward reads.
        output_r (str): Output FASTQ file for reverse reads.
        include_read_id (bool): Whether to include read identifier in the "+" separator line.
    """
    # Calculate insert size based on read length and overlap
    insert_size = 2 * read_length - overlap

    # Define fixed fields for the header (customize as needed)
    instrument = "M06023"
    run_number = "50"
    flowcell_id = "000000000-C2KVH"
    lane = "1"

    # Initialize read pair counter for unique positions
    tile = 1101
    x_pos = 1000
    y_pos = 1000

    try:
        with open(fasta_file, 'r') as fasta, \
             open(output_l, 'w') as fastq_l, \
             open(output_r, 'w') as fastq_r:

            for record in SeqIO.parse(fasta, 'fasta'):
                seq = str(record.seq).upper()
                seq_len = len(seq)

                if seq_len < insert_size:
                    print(f"Warning: Sequence {record.id} is shorter than the calculated insert size ({insert_size}). Skipping.", file=sys.stderr)
                    continue

                # Calculate step size for sliding window
                step_size = read_length - overlap

                for i in range(0, seq_len - insert_size + 1, step_size):
                    # Extract forward and reverse reads
                    forward_seq = seq[i:i + read_length]
                    reverse_seq = seq[i + insert_size - read_length:i + insert_size]
                    reverse_seq_rc = str(Seq(reverse_seq).reverse_complement())

                    if len(forward_seq) < read_length or len(reverse_seq_rc) < read_length:
                        continue  # Skip incomplete reads at the end

                    # Generate quality strings
                    forward_quality = generate_quality_scores(read_length)
                    reverse_quality = generate_quality_scores(read_length)

                    # Generate formatted headers following Illumina standards
                    read_id = f"{instrument}:{run_number}:{flowcell_id}:{lane}:{tile}:{x_pos}:{y_pos}"
                    read_id_forward = f"{read_id} 1:N:0:{record.id}"
                    read_id_reverse = f"{read_id} 2:N:0:{record.id}"

                    # Write to FASTQ files
                    if include_read_id:
                        fastq_l.write(f"@{read_id_forward}\n{forward_seq}\n+{read_id_forward}\n{forward_quality}\n")
                        fastq_r.write(f"@{read_id_reverse}\n{reverse_seq_rc}\n+{read_id_reverse}\n{reverse_quality}\n")
                    else:
                        fastq_l.write(f"@{read_id_forward}\n{forward_seq}\n+\n{forward_quality}\n")
                        fastq_r.write(f"@{read_id_reverse}\n{reverse_seq_rc}\n+\n{reverse_quality}\n")

                    # Update positions for uniqueness
                    x_pos += 1
                    y_pos += 1

        print(f"Paired-end FASTQ files generated from FASTA:\n - {output_l}\n - {output_r}")
    except Exception as e:
        print(f"Error during FASTA to FASTQ conversion: {e}", file=sys.stderr)
        sys.exit(1)

def main():
    parser = argparse.ArgumentParser(description='Convert a FASTA file to paired-end FASTQ files with properly formatted headers and actual quality scores.')

    parser.add_argument('input_file', help='Input FASTA file')
    parser.add_argument('--read_length', type=int, default=150, help='Length of each read (default: 150)')
    parser.add_argument('--overlap', type=int, default=50, help='Overlap between forward and reverse reads (default: 50)')
    parser.add_argument('--output_l', default='fastq_l.fastq', help='Output FASTQ file for forward reads (default: fastq_l.fastq)')
    parser.add_argument('--output_r', default='fastq_r.fastq', help='Output FASTQ file for reverse reads (default: fastq_r.fastq)')
    parser.add_argument('--include_read_id', action='store_true', help='Include read identifier in the "+" separator line (default: False)')

    args = parser.parse_args()

    input_file = args.input_file
    output_l = args.output_l
    output_r = args.output_r
    include_read_id = args.include_read_id

    if not os.path.isfile(input_file):
        print(f"Error: Input file '{input_file}' does not exist.", file=sys.stderr)
        sys.exit(1)

    if is_fasta(input_file):
        print(f"Detected FASTA input: {input_file}")
        fasta_to_paired_fastq(
            fasta_file=input_file,
            read_length=args.read_length,
            overlap=args.overlap,
            output_l=output_l,
            output_r=output_r,
            include_read_id=include_read_id
        )
    else:
        print(f"Error: The input file '{input_file}' is not in FASTA format.", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
