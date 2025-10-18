# HYMET Test Dataset Toolkit

This directory contains helper scripts used in the paper to synthesize small test datasets and sanity-check the HYMET pipeline. They are not wired into the main case-study or benchmark runners yet, so here’s how to use them manually.

## Requirements

```bash
pip install wget biopython
```

The scripts expect to run from `HYMET/testdataset/` and create output under `HYMET/data/testdataset/`.

## Workflow Overview

1. `filterGCF.py` – filters downloaded GCF FASTA files (requires Biopython).  
   Usage example: `python filterGCF.py input_dir output_dir`

2. `mutationGCF.py` – introduces synthetic mutations into FASTA sequences.  
   Usage: `python mutationGCF.py input.fna output.fna --mutation-rate 0.01`

3. `extractTaxonomy.py` – pulls taxonomy metadata via NCBI Entrez (set `ENTREZ_EMAIL`).

4. `extractNC.py` – extracts nucleotide segments based on regex patterns.

5. `createDatabase.py` – orchestrates the full mini-database build (downloads refs via wget, filters, mutates).  
   Run only after installing dependencies; ensure output directories exist.

These scripts remain standalone until we integrate them into the main pipeline. For the journal revision we should capture representative commands, outputs, and note any external dependencies (NCBI APIs, etc.) in the supplementary material.
