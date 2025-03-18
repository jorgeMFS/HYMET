import os
import re
from Bio import SeqIO

def get_directory_path(prompt):
    while True:
        path = input(prompt)
        if os.path.isdir(path):
            return path
        else:
            print("The provided path is not a valid directory. Please try again.")

def main():
    input_directory = get_directory_path("Enter the path to the input directory: ")

    output_file = input("Enter the path to the output file (e.g., combined_archaeagenomes.fasta): ")

    sequences = []

    id_pattern = re.compile(r'N[A-Z]_\S+')

    for filename in os.listdir(input_directory):
        if filename.startswith('GCF_') and filename.endswith('.fna'):
            with open(os.path.join(input_directory, filename), 'r') as handle:
                for record in SeqIO.parse(handle, 'fasta'):
                    match = id_pattern.search(record.description)
                    
                    if match:
                        new_id = match.group(0)
                        record.id = new_id
                        record.description = ''
                        sequences.append(record)
                    else:
                        print(f"Warning: No valid identifier found for {filename}")
                    break  

    SeqIO.write(sequences, output_file, 'fasta')

    print(f"File {output_file} created successfully!")

if __name__ == "__main__":
    main()
