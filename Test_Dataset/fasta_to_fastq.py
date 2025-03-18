import sys
import argparse

# Character to use for quality
quality_char = 'I'

def convert_fasta_to_fastq(input_file, output_file):
    with open(input_file, 'r') as fasta, open(output_file, 'w') as fastq:
        sequence = ""
        for line in fasta:
            line = line.strip()
            if line.startswith('>'):
                if sequence:  # Write the previous sequence if it exists
                    fastq.write(f'{sequence}\n')  # Sequence line
                    fastq.write('+\n')  # Separator
                    fastq.write(f'{"".join([quality_char] * len(sequence))}\n')  # Quality line
                    sequence = ""
                # Convert FASTA header to FASTQ header
                fastq.write(f'@{line[1:]}\n')  # Skip '>'
            else:
                sequence += line
        
        # Write the last sequence
        if sequence:
            fastq.write(f'{sequence}\n')  # Sequence line
            fastq.write('+\n')  # Separator
            fastq.write(f'{"".join([quality_char] * len(sequence))}\n')  # Quality line

    print("Conversion complete!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Convert FASTA to FASTQ')
    parser.add_argument('-i', '--input', help='Input FASTA file', required=True)
    parser.add_argument('-o', '--output', help='Output FASTQ file', required=True)
    args = parser.parse_args()

    convert_fasta_to_fastq(args.input, args.output)
