#!/usr/bin/env python3

import argparse
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

def main():
    parser = argparse.ArgumentParser(
        description="Generate taxonomy hierarchy TSV from NCBI taxonomy dumps."
    )
    parser.add_argument(
        "taxonomy_dir",
        help="Directory containing names.dmp and nodes.dmp files",
    )
    parser.add_argument(
        "output_dir",
        help="Directory where taxonomy_hierarchy.tsv will be written",
    )

    args = parser.parse_args()

    taxonomy_dir = os.path.abspath(args.taxonomy_dir)
    output_dir = os.path.abspath(args.output_dir)

    if not os.path.isdir(taxonomy_dir):
        raise NotADirectoryError(f"Taxonomy directory not found: {taxonomy_dir}")

    os.makedirs(output_dir, exist_ok=True)

    names_dmp_path = os.path.join(taxonomy_dir, 'names.dmp')
    nodes_dmp_path = os.path.join(taxonomy_dir, 'nodes.dmp')
    output_hierarchy_file_path = os.path.join(output_dir, 'taxonomy_hierarchy.tsv')

    if not os.path.exists(names_dmp_path):
        raise FileNotFoundError(f"File {names_dmp_path} not found.")

    if not os.path.exists(nodes_dmp_path):
        raise FileNotFoundError(f"File {nodes_dmp_path} not found.")

    generate_taxonomy_hierarchy(names_dmp_path, nodes_dmp_path, output_hierarchy_file_path)


if __name__ == "__main__":
    main()