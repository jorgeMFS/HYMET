import csv
from pathlib import Path
import logging
from Bio import Entrez

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def read_assembly_summary(file_path):
    taxonomy_dict = {}
    with open(file_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 20:
                continue
            assembly_accession = fields[0]
            species_taxid = fields[6]
            organism_name = fields[7]
            ftp_path = fields[19].split('/')[-1]
            taxonomy_dict[ftp_path] = {
                'assembly_accession': assembly_accession,
                'species_taxid': species_taxid,
                'organism_name': organism_name
            }
    return taxonomy_dict

def get_taxonomy(taxid):
    handle = Entrez.efetch(db="taxonomy", id=str(taxid), retmode="xml")
    records = Entrez.read(handle)
    if records:
        lineage = records[0]["LineageEx"]
        full_lineage = {rank["Rank"]: rank["ScientificName"] for rank in lineage}
        full_lineage["species"] = records[0]["ScientificName"]
        return full_lineage
    return None

def create_full_taxonomy_mapping(domain_path, summary_file, output_csv):
    logging.info(f"Processing {domain_path.name} assembly summary")
    taxonomy = read_assembly_summary(summary_file)

    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['accession', 'domain', 'organism_name', 'species_taxid', 'assembly_accession',
                      'Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for file in domain_path.iterdir():
            if file.suffix == '.fna':
                file_prefix = file.stem.split('_genomic')[0]
                if file_prefix in taxonomy:
                    row = {'accession': file_prefix, 'domain': domain_path.name}
                    row.update(taxonomy[file_prefix])
                    
                    full_taxonomy = get_taxonomy(row['species_taxid'])
                    if full_taxonomy:
                        row['Kingdom'] = full_taxonomy.get('superkingdom', '')
                        row['Phylum'] = full_taxonomy.get('phylum', '')
                        row['Class'] = full_taxonomy.get('class', '')
                        row['Order'] = full_taxonomy.get('order', '')
                        row['Family'] = full_taxonomy.get('family', '')
                        row['Genus'] = full_taxonomy.get('genus', '')
                        row['Species'] = full_taxonomy.get('species', '')
                    else:
                        logging.warning(f"No full taxonomy found for {file_prefix}")
                    
                    writer.writerow(row)
                else:
                    logging.warning(f"No taxonomy found for {file_prefix}")

    logging.info(f"Full taxonomy mapping for {domain_path.name} has been written to {output_csv}")

def main():
    Entrez.email = input("Enter your email (required for NCBI E-utilities): ")

    domain_directory_input = input("Enter the path to the directory containing domain .fna files: ")
    summary_file_input = input("Enter the path to the assembly summary file: ")

    domain_path = Path(domain_directory_input)
    summary_file = Path(summary_file_input)

    if not domain_path.is_dir():
        print("The provided path for domain files is not a valid directory. Please check and try again.")
        return
    
    if not summary_file.is_file():
        print("The provided path for the assembly summary file is not a valid file. Please check and try again.")
        return

    output_csv = domain_path / f"{domain_path.name}_taxonomy_mapping.csv"
    
    create_full_taxonomy_mapping(domain_path, summary_file, output_csv)

if __name__ == "__main__":
    main()
