# ZymoGut D6331 Case Study

Ground-truth validation using the ZymoBIOMICS Gut Microbiome Standard (D6331) to address Reviewer #1's comment.

## Overview

This self-contained case study processes Oxford Nanopore data from the ZymoBIOMICS Gut Microbiome Standard, assembles contigs, runs HYMET classification, and compares against the known ground truth.

### Dataset

| Property | Value |
|----------|-------|
| Product | ZymoBIOMICS Gut Microbiome Standard (D6331) |
| Manufacturer | Zymo Research |
| Composition | 21 strains / 15 species |
| Sequencing | Oxford Nanopore SUP basecalls |
| Source | MicroBench (Kirk3gaard) |

### Expected Species

| Species | Abundance (%) |
|---------|--------------|
| Faecalibacterium prausnitzii | 17.63 |
| Veillonella rogosae | 15.87 |
| Bacteroides fragilis | 9.94 |
| Roseburia hominis | 9.89 |
| Lactobacillus fermentum | 9.63 |
| Bifidobacterium adolescentis | 8.78 |
| Fusobacterium nucleatum | 7.49 |
| Saccharomyces cerevisiae | 6.09 |
| Candida albicans | 5.59 |
| Prevotella corporis | 4.98 |
| Clostridioides difficile | 2.62 |
| Akkermansia muciniphila | 0.97 |
| + 3 low-abundance species | <0.1 |

## Quick Start

```bash
cd case/zymogut

# Full pipeline (~2-3 hours)
./run_zymogut.sh

# Or run phases individually:
./run_zymogut.sh --phase 1     # Download (~10 min, ~5GB)
./run_zymogut.sh --phase 2     # Assembly (~30-60 min)
./run_zymogut.sh --phase 3     # Build truth (~5 min)
./run_zymogut.sh --phase 4     # HYMET classification (~10-30 min)
./run_zymogut.sh --phase 5     # Analysis (~2 min)
./run_zymogut.sh --phase 6 --publish  # Package for paper
```

## Requirements

- HYMET initialized (`bin/hymet init`)
- Flye >= 2.9 (`mamba install -c bioconda flye`)
- ~20GB RAM, ~50GB disk space

## Directory Structure

```
case/zymogut/
├── run_zymogut.sh          # Master orchestrator
├── README.md               # This file
├── config/
│   └── zymogut_d6331.yml   # Ground truth configuration
├── scripts/
│   ├── 01_download.sh      # Download reads + references
│   ├── 02_assemble.sh      # Flye assembly
│   ├── 03_build_truth.sh   # Build ground truth
│   ├── 04_classify.sh      # HYMET classification
│   ├── 05_analyse.sh       # Analysis & figures
│   └── 06_package.sh       # Package results
└── work/                   # Working directory (gitignored)
    ├── downloads/          # Raw data
    ├── assembly/           # Flye output
    ├── truth/              # Ground truth files
    ├── hymet_output/       # HYMET results
    └── results/            # Analysis & figures
```

## Output

After `--publish`, results are integrated into:

- `case/truth/zymogut_d6331/` - Ground truth files
- `results/cases/zymogut/run_*/` - Classification results
- `case/manifest_zymogut.tsv` - Manifest entry

### Key Metrics

| Metric | Description |
|--------|-------------|
| L1 Distance | Sum of absolute abundance differences |
| Bray-Curtis | Ecological dissimilarity |
| Correlation | Pearson correlation |

### Figures for Paper

- `fig_zymogut_abundance_comparison.png` - Bar chart comparison
- `fig_zymogut_correlation.png` - Scatter plot

## Citation

- HYMET (this paper)
- Zymo Research D6331
- Kirk3gaard/MicroBench
