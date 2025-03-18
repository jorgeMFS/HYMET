import os
import glob
from Bio import SeqIO

def process_gcf_file(file_path):
    sequences = list(SeqIO.parse(file_path, "fasta"))
    total_bases = sum(len(seq) for seq in sequences)
    target_bases = int(total_bases * 0.1)  # 10% of total bases
    
    if len(sequences) == 1:
        # Case with only one sequence
        seq = sequences[0]
        trimmed_seq = seq[:target_bases]
        trimmed_seq.id = seq.id
        trimmed_seq.description = f"First 10% segment (1-{target_bases})"
        return [trimmed_seq]
    else:
        # Case with multiple sequences
        trimmed_sequences = []
        bases_per_seq = target_bases // len(sequences)
        remaining_bases = target_bases % len(sequences)
        
        for i, seq in enumerate(sequences):
            if i == len(sequences) - 1:
                bases_to_extract = bases_per_seq + remaining_bases
            else:
                bases_to_extract = bases_per_seq
            
            if len(seq) > bases_to_extract:
                trimmed_seq = seq[:bases_to_extract]
                trimmed_seq.id = seq.id
                trimmed_seq.description = f"First {bases_to_extract} bases"
            else:
                trimmed_seq = seq[:]
                trimmed_seq.id = seq.id
                trimmed_seq.description = "Full sequence (shorter than target segment)"
            trimmed_sequences.append(trimmed_seq)
        
        return trimmed_sequences

def process_gcf_files(input_dir, output_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    for file_path in glob.glob(os.path.join(input_dir, "GCF_*.fna")):
        file_name = os.path.basename(file_path)
        output_path = os.path.join(output_dir, file_name)
        
        trimmed_sequences = process_gcf_file(file_path)
        
        with open(output_path, "w") as output_file:
            SeqIO.write(trimmed_sequences, output_file, "fasta")
        
        print(f"Processed: {file_name}")

def main():
    # Prompt user for input and output directory paths
    input_directory = input("Enter the path to the input directory containing GCF files: ")
    output_directory = input("Enter the path to the output directory: ")

    # Validate paths
    if not os.path.isdir(input_directory):
        print("The provided input path is not a valid directory. Please check and try again.")
        return
    
    process_gcf_files(input_directory, output_directory)

if __name__ == "__main__":
    main()
