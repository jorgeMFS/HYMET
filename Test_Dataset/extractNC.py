import os
import re
import csv
from pathlib import Path
from collections import defaultdict

def extract_identifiers(file_path):
    gc_id = re.search(r'(GCF_\d+\.\d+)', file_path.name)
    gc_id = gc_id.group(1) if gc_id else None

    ids = []
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('>'):
                id_match = re.search(r'(N[A-Z]_\w+\.\d+)', line)
                if id_match:
                    ids.append(id_match.group(1))

    return gc_id, ids

def create_mapping(domain_directory):
    mapping = defaultdict(list)
    print(f"Processing directory: {domain_directory}")
    for file in os.listdir(domain_directory):
        if file.endswith('.fna'):
            file_path = Path(domain_directory) / file
            gc_id, ids = extract_identifiers(file_path)
            if gc_id and ids:
                mapping[gc_id].extend(ids)
                print(f"Mapped: {gc_id} -> {', '.join(ids)}")
            else:
                print(f"Failed to map: {file}")
    print(f"Mapping created: {len(mapping)} entries")
    return mapping

def write_mapping_to_csv(mapping, output_file):
    with open(output_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['GCF_ID', 'Sequence_IDs'])
        for gcf_id, ids in mapping.items():
            writer.writerow([gcf_id, ';'.join(ids)])

def main():
    # Prompt user for the domain directory
    domain_directory = input("Enter the path to the domain directory (e.g., archaea): ")
    
    # Check if the provided path is valid
    if not os.path.isdir(domain_directory):
        print("The provided path is not a valid directory. Please check and try again.")
        return

    output_file = Path(domain_directory) / "gcf_sequence_id_mapping.csv"

    mapping = create_mapping(domain_directory)
    write_mapping_to_csv(mapping, output_file)

    print(f"GCF to Sequence ID mapping has been written to {output_file}")
    print(f"Total number of GCF entries: {len(mapping)}")

if __name__ == "__main__":
    main()
