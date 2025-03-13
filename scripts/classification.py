#!/usr/bin/env python3
import os
import csv
from collections import defaultdict
import argparse
from multiprocessing import Pool
import logging
import sys

csv.field_size_limit(sys.maxsize)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def load_taxonomy_file(taxonomy_file):
    taxonomy = {}
    with open(taxonomy_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            taxid = row["TaxID"]
            for identifier in row["Identifiers"].split(";"):
                cleaned_id = identifier.strip()
                if cleaned_id:
                    taxonomy[cleaned_id] = taxid
    logging.info(f"Loaded {len(taxonomy)} taxonomy mappings")
    return taxonomy

def load_taxonomy_hierarchy_file(taxonomy_hierarchy_file):
    hierarchy = {}
    with open(taxonomy_hierarchy_file, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            taxid = row["TaxID"]
            hierarchy[taxid] = row["Lineage"].strip()
    logging.info(f"Loaded {len(hierarchy)} taxonomy hierarchies")
    return hierarchy

def parse_mashmap_file(mashmap_file):
    query_map = defaultdict(list)
    ref_counts = defaultdict(int)

    with open(mashmap_file, "r") as f:
        for line in f:
            parts = line.strip().split("\t")
            if len(parts) < 14:  # Ensure enough columns exist
                continue

            query_id = parts[0]       # Query sequence name
            query_len = int(parts[1])  # Query sequence length
            ref_id = parts[5]         # Reference sequence name

            start_query = int(parts[2])  # Query start position
            end_query = int(parts[3])    # Query end position

            align_len = end_query - start_query  # Calculate alignment length

            identity_str = parts[12]  # e.g., "id:f:1"
            identity = float(identity_str.split(":")[2]) if "id:f:" in identity_str else 0.0


            coverage = align_len / query_len if query_len > 0 else 0
            is_exact = (identity >= 0.99) and (coverage >= 0.9)  # Adjust criteria as needed

            query_map[query_id].append((ref_id, coverage, is_exact))
            ref_counts[ref_id] += 1

    logging.info(f"Processed {len(query_map)} queries from Mashmap file")
    return query_map, ref_counts

def determine_taxonomic_level(lineage):
    rank_order = [
        'superkingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'strain'
    ]
    
    current_level = None
    for part in lineage.split(';'):
        part = part.strip()
        if ':' in part:
            rank, name = part.split(':', 1)
            rank = rank.strip().lower()
            if rank in rank_order:
                if current_level is None:
                    current_level = rank
                else:
                    current_index = rank_order.index(current_level)
                    new_index = rank_order.index(rank)
                    if new_index > current_index:
                        current_level = rank
    return current_level if current_level is not None else 'root'

def calculate_weighted_lineage(refs, ref_abundance, taxonomy):
    taxid_weights = defaultdict(float)
    total_weight = 0.0
    
    for ref_id, coverage, _ in refs:
        if ref_id not in taxonomy:
            continue
            
        taxid = taxonomy[ref_id]
        weight = coverage * ref_abundance.get(ref_id, 1)
        taxid_weights[taxid] += weight
        total_weight += weight
    
    return taxid_weights, total_weight

def determine_lca(taxid_weights, total_weight, taxonomy_hierarchy):
    if total_weight == 0:
        return "Unknown", "root", 0.0

    lineages = []
    for taxid, weight in taxid_weights.items():
        if taxid in taxonomy_hierarchy:
            lineage = taxonomy_hierarchy[taxid].split(";")
            lineages.append((lineage, weight/total_weight))

    if not lineages:
        return "Unknown", "root", 0.0

    consensus = {}
    confidence = 1.0
    rank_order = [
        'superkingdom', 'phylum', 'class', 'order',
        'family', 'genus', 'species', 'strain'
    ]

    for rank in rank_order:
        level_counts = defaultdict(float)
        for lineage, weight in lineages:
            for part in lineage:
                if part.startswith(f"{rank}:"):
                    level_counts[part] += weight
                    break
        
        if level_counts:
            best_match, conf = max(level_counts.items(), key=lambda x: x[1])
            consensus[rank] = best_match
            confidence *= conf
        else:
            break  # Stop at first missing rank

    lineage_parts = [consensus.get(rank) for rank in rank_order if consensus.get(rank)]
    if not lineage_parts:
        return "Unknown", "root", 0.0

    full_lineage = ";".join(lineage_parts)
    level = determine_taxonomic_level(full_lineage)
    return full_lineage, level, min(confidence, 1.0)

def process_query(args):
    query, refs, ref_abundance, taxonomy, taxonomy_hierarchy = args
    
    # Check for exact matches first
    exact_matches = [ref for ref, _, is_exact in refs if is_exact and ref in taxonomy]
    if exact_matches:
        taxid = taxonomy[exact_matches[0]]
        if taxid in taxonomy_hierarchy:
            lineage = taxonomy_hierarchy[taxid]
            level = determine_taxonomic_level(lineage)
            return (query, lineage, level, 1.0)

    # Calculate LCA for non-exact matches
    taxid_weights, total_weight = calculate_weighted_lineage(refs, ref_abundance, taxonomy)
    lineage, level, confidence = determine_lca(taxid_weights, total_weight, taxonomy_hierarchy)
    
    return (query, lineage, level, confidence)

def main_process(mashmap_file, taxonomy_file, hierarchy_file, output_file, processes=4):
    taxonomy = load_taxonomy_file(taxonomy_file)
    taxonomy_hierarchy = load_taxonomy_hierarchy_file(hierarchy_file)
    query_map, ref_abundance = parse_mashmap_file(mashmap_file)

    tasks = [(query, refs, ref_abundance, taxonomy, taxonomy_hierarchy) 
             for query, refs in query_map.items()]

    results = []
    with Pool(processes) as pool:
        results = pool.map(process_query, tasks)

    with open(output_file, 'w') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(['Query', 'Lineage', 'Taxonomic Level', 'Confidence'])
        
        classified = 0
        for query, lineage, level, confidence in results:
            if lineage != 'Unknown':
                classified += 1
            writer.writerow([query, lineage, level, f"{confidence:.4f}"])

    logging.info(f"Classification complete. Results saved to {output_file}")
    logging.info(f"Classified: {classified}/{len(results)} ({classified/len(results):.1%})")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Advanced LCA/Best Match Taxonomic Classifier for Mashmap Output")
    parser.add_argument("--mashmap", required=True, help="Input Mashmap file")
    parser.add_argument("--taxonomy", required=True, help="Taxonomy mapping file")
    parser.add_argument("--hierarchy", required=True, help="Taxonomy hierarchy file")
    parser.add_argument("--output", required=True, help="Output TSV file")
    parser.add_argument("--processes", type=int, default=4, help="Number of parallel processes")
    
    args = parser.parse_args()
    
    main_process(
        args.mashmap,
        args.taxonomy,
        args.hierarchy,
        args.output,
        args.processes
    )