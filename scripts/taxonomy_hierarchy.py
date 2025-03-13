#!/usr/bin/env python3

import os
import csv

def parse_names_dmp(file_path):
    """
    Parse the names.dmp file to create a dictionary mapping TaxID to Scientific Name.
    """
    taxid_to_name = {}
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = parts[0].strip()
            name = parts[1].strip()
            name_class = parts[3].strip("\t|\n")
            if name_class == "scientific name":  # Filter only scientific names
                taxid_to_name[taxid] = name
    return taxid_to_name

def parse_nodes_dmp(file_path):
    """
    Parse the nodes.dmp file to create a dictionary mapping TaxID to (Rank, ParentTaxID).
    """
    taxid_to_hierarchy = {}
    with open(file_path, "r", encoding="utf-8") as f:
        for line in f:
            parts = line.split("\t|\t")
            taxid = parts[0].strip()
            parent_taxid = parts[1].strip()
            rank = parts[2].strip("\t|\n")
            if rank == "no rank" and "strain" in parts[4]:
                rank = "strain"
            taxid_to_hierarchy[taxid] = (rank, parent_taxid)
    return taxid_to_hierarchy

def generate_taxonomy_hierarchy(names_file, nodes_file, output_file):
    """
    Generate the taxonomy_hierarchy.tsv file with TaxID, Scientific Name, Rank, ParentTaxID, and Full Lineage.
    """
    print("Loading data from files...")
    taxid_to_name = parse_names_dmp(names_file)
    taxid_to_hierarchy = parse_nodes_dmp(nodes_file)

    def get_lineage(taxid):
        lineage = []
        current_taxid = taxid
        while current_taxid != "1":
            rank, parent_taxid = taxid_to_hierarchy.get(current_taxid, ("", "1"))
            name = taxid_to_name.get(current_taxid, "Unknown")
            lineage.insert(0, f"{rank}:{name}")
            current_taxid = parent_taxid
        return ";".join(lineage)

    print("Generating taxonomy_hierarchy.tsv...")
    with open(output_file, "w", encoding="utf-8") as out:
        out.write("TaxID\tName\tRank\tParentTaxID\tLineage\n")
        for taxid, (rank, parent_taxid) in taxid_to_hierarchy.items():
            name = taxid_to_name.get(taxid, "Unknown")
            lineage = get_lineage(taxid)
            out.write(f"{taxid}\t{name}\t{rank}\t{parent_taxid}\t{lineage}\n")

    print(f"File generated successfully: {output_file}")

if __name__ == "__main__":
    BASE_PATH = '/mnt/storagelv/home/inesbrancomartins/Tese/tool1_fast'
    TAXONOMY_FILES_DIR = os.path.join(BASE_PATH, 'taxonomy_files')
    SCRIPTS_DIR = os.path.join(BASE_PATH, 'scripts')
    
    # Create scripts directory if it does not exist
    if not os.path.exists(SCRIPTS_DIR):
        os.makedirs(SCRIPTS_DIR)
    
    NAMES_DMP_PATH = os.path.join(TAXONOMY_FILES_DIR, 'names.dmp')
    NODES_DMP_PATH = os.path.join(TAXONOMY_FILES_DIR, 'nodes.dmp')
    OUTPUT_HIERARCHY_FILE_PATH = os.path.join(SCRIPTS_DIR, 'taxonomy_hierarchy.tsv')
    
    # Check if files exist
    if not os.path.exists(NAMES_DMP_PATH):
        raise FileNotFoundError(f"File {NAMES_DMP_PATH} not found.")
    
    if not os.path.exists(NODES_DMP_PATH):
        raise FileNotFoundError(f"File {NODES_DMP_PATH} not found.")

    # Generate taxonomy hierarchy
    generate_taxonomy_hierarchy(NAMES_DMP_PATH, NODES_DMP_PATH, OUTPUT_HIERARCHY_FILE_PATH)
