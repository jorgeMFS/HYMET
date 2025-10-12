# HYMET (Hybrid Metagenomic Tool)

[![Conda Version](https://anaconda.org/bioconda/hymet/badges/version.svg)](https://anaconda.org/bioconda/hymet)
[![Downloads](https://anaconda.org/bioconda/hymet/badges/downloads.svg)](https://anaconda.org/bioconda/hymet)
[![Platforms](https://anaconda.org/bioconda/hymet/badges/platforms.svg)](https://anaconda.org/bioconda/hymet)
[![License](https://anaconda.org/bioconda/hymet/badges/license.svg)](https://anaconda.org/bioconda/hymet)
[![Latest Release Date](https://anaconda.org/bioconda/hymet/badges/latest_release_date.svg)](https://anaconda.org/bioconda/hymet)

## Installation and Configuration

Follow the steps below to install and configure **HYMET**

### 1. Installation with Conda (Recommended)

The easiest way to install HYMET is through Bioconda:

```bash
conda install -c bioconda hymet
```

After installation, you will need to download the reference databases as described in the [Reference Sketched Databases](#reference-sketched-databases) section.

### 2. Clone the Repository

Alternatively, you can clone the repository to your local environment:

```bash
git clone https://github.com/ieeta-pt/HYMET.git
cd HYMET
```

### 3. Installation with Docker

If you prefer using Docker, follow these steps:

1. **Build the Docker Image**:
   ```bash
   docker build -t hymet .
   ```

2. **Run the Container**:
   ```bash
   docker run -it hymet
   ```

3. **Inside the Container**:
   - The environment will already be set up with all dependencies installed.
   - Run the tool as needed.

### 4. Installation with Conda Environment File

If you cloned the repository, you can create a Conda environment from the included file:

1. **Create the Conda Environment**:
   ```bash
   conda env create -f environment.yml
   ```

2. **Activate the Environment**:
   ```bash
   conda activate hymet_env
   ```

## Input Requirements

### Input File Format
The tool expects input files in **FASTA format** (`.fna` or `.fasta`). Each file should contain metagenomic sequences with headers in the following format:
```
>sequence_id additional_info
SEQUENCE_DATA
```
- **sequence_id**: A unique identifier for the sequence.
- **additional_info**: Optional metadata (e.g., source organism, length).
- **SEQUENCE_DATA**: The nucleotide sequence.

### Input Directory
Place your input files in the directory specified by the `$input_dir` variable in the `main.pl` script. 

For example, if your input directory contains the following files:
```
input_dir/
├── sample1.fna
├── sample2.fna
└── sample3.fna
```
Each file (`sample1.fna`, `sample2.fna`, etc.) should follow the FASTA format described above.


## Execution
Ensure all scripts have execution permissions:
   ```bash
   chmod +x config.pl
   chmod +x main.pl
   chmod +x scripts/*.sh
   chmod +x scripts/*.py
   ```

After installation and configuration, first you should run the configuration script to download and prepare the taxonomy files, and define the main paths.

```bash
./config.pl
```

Then, you can run the main tool to perform taxonomic identification:

```bash
./main.pl
```

If installed via Conda, you can use:
```bash
hymet-config
hymet
```

## Reference Sketched Databases

The databases required to run the tool are available for download on Google Drive:
- **sketch1.msh.gz**: [Download](https://drive.google.com/drive/folders/1YC0N77UUGinFHNbLpbsucu1iXoLAM6lm?usp=share_link)
- **sketch2.msh.gz**: [Download](https://drive.google.com/drive/folders/1YC0N77UUGinFHNbLpbsucu1iXoLAM6lm?usp=share_link)
- **sketch3.msh.gz**: [Download](https://drive.google.com/drive/folders/1YC0N77UUGinFHNbLpbsucu1iXoLAM6lm?usp=share_link)

### Steps to Use the Databases

1. **Download the Files**:
   - Click on the links above to download the `.gz` files.

2. **Place the Files in the `data/` Directory**:
   - Move the downloaded files to the `data/` directory of the project.

3. **Unzip the Files**:
   - Use the following command to unzip the `.gz` files:
     ```bash
     gunzip data/sketch1.msh.gz
     gunzip data/sketch2.msh.gz
     gunzip data/sketch3.msh.gz
     ```
   - This will extract the files `sketch1.msh`, `sketch2.msh`, and `sketch3.msh` in the `data/` directory.

4. **Verify the Files**:
   - After unzipping, ensure the files are in the correct format and location:
     ```bash
     ls data/
     ```
   - You should see the files `sketch1.msh`, `sketch2.msh`, and `sketch3.msh`.


##  Project Structure

- **config.pl**: Configuration script that downloads and prepares taxonomy files.
- **main.pl**: Main script that runs the taxonomic identification pipeline.
- **scripts/**: Directory containing helper scripts in Perl, Python, and Bash.
  - **mash.sh**: Script to run Mash.
  - **downloadDB.py**: Script to download genomes.
  - **minimap.sh**: Script to run Minimap2.
  - **classification.py**: Script for taxonomic classification.
- **taxonomy_files/**: Directory containing downloaded taxonomy files.
- **data/**: Directory for storing intermediate data.
  - sketch1.msh
  - sketch2.msh
  - sketch3.msh
  - taxonomy_hierarchy.tsv
- **output/**: Directory where final results are saved.

## Example Output

The tool generates a `classified_sequences.tsv` file in the `output/` directory with the following columns:

- **Query**: Identifier of the queried sequence.
- **Lineage**: Identified taxonomic lineage.
- **Taxonomic Level**: Taxonomic level (e.g., species, genus).
- **Confidence**: Classification confidence (0 to 1).

## Test Dataset
- This folder includes scripts to install and prepare all necessary data to replicate the work using our dataset.
  - **Prerequisites**:
    - Before running the scripts in this folder, users need to download the assembly files (`assembly_files.txt`) for each domain from the NCBI FTP site.
  - **Scripts**:
    - **`create_database.py`**: Downloads 10% of the content from each downloaded assembly file and organizes the datasets by domain.
    - **`extractNC.py`**: Maps the content of each Genome Collection File (GCF) with its respective sequence identifiers. It generates a CSV file containing       this mapping, with one column for the GCF and another column for the sequence identifiers (such as NC, NZ, etc.) present in each GCF.
    - **`extractTaxonomy.py`**: Creates a CSV file containing the GCF and its respective taxonomy, among other information.
    - Additional scripts modify the data format and organization, including:
      - Implementing mutations
      - Converting formats (e.g., FASTA to FASTQ)
      - Formatting into paired-end reads
    - **`GCFtocombinedfasta.py`**: Combines all GCFs from each domain into a single FASTA file, separating sequences by identifier. This script is used as input for most of the tools.

## CAMI Benchmarking

HYMET now includes comprehensive benchmarking capabilities against CAMI (Critical Assessment of Metagenome Interpretation) datasets. This allows you to evaluate HYMET's performance against standardized ground truth data.

### Prerequisites for CAMI Benchmarking

1. **CAMI Dataset**: Download a CAMI dataset (e.g., CAMI Challenge datasets)
2. **Required Tools**: Ensure the following tools are available in `/data/tools/`:
   - `hymet2cami_fixed.py` - Converts HYMET output to CAMI format
   - `build_id_map.py` - Builds ID-to-TaxID mappings for fallback classification
   - `mini_classify.py` - Fallback classification method
3. **Taxonomy Database**: Ensure `taxonomy_files/` directory is properly configured

### Quick Start: CAMI Benchmarking

```bash
# Basic usage with default parameters
bash benchmark_cami.sh

# Custom parameters
INPUT_FASTA=/path/to/your/sample.fna \
OUTDIR=/path/to/output \
THREADS=16 \
bash benchmark_cami.sh
```

### Detailed Usage

#### 1. Basic Benchmarking

```bash
# Run with default CAMI sample_0 dataset
cd /data/HYMET
bash benchmark_cami.sh
```

This will:
- Run the complete HYMET pipeline on the input FASTA
- Generate CAMI-format taxonomic profiles
- Evaluate against CAMI ground truth
- Produce comprehensive performance metrics

#### 2. Custom Input and Output

```bash
# Specify custom input and output locations
INPUT_FASTA=/data/cami/your_sample.fna \
OUTDIR=/data/results/your_sample \
THREADS=8 \
bash benchmark_cami.sh
```

#### 3. Custom Ground Truth Files

```bash
# Use custom CAMI ground truth files
TRUTH_CONTIGS_1=/path/to/gsa_mapping.tsv \
TRUTH_PROFILE=/path/to/taxonomic_profile.txt \
bash benchmark_cami.sh
```

### Configuration Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `INPUT_FASTA` | `/data/cami/sample_0.fna` | Input FASTA file to analyze |
| `OUTDIR` | `/data/hymet_out/sample_0` | Output directory for results |
| `THREADS` | `16` | Number of parallel threads |
| `TRUTH_CONTIGS_1` | `/data/cami/sample_0/.../gsa_mapping.tsv` | Primary ground truth contig mapping |
| `TRUTH_CONTIGS_2` | `/data/cami/sample_0/.../gsa_mapping_new.tsv` | Alternative ground truth contig mapping |
| `TRUTH_PROFILE` | `/data/cami/sample_0/taxonomic_profile_0.txt` | Ground truth taxonomic profile |

### Output Files

The benchmarking process generates several output files in the `$OUTDIR/eval/` directory:

#### Performance Metrics
- **`profile_summary.tsv`** - Profile-level metrics by taxonomic rank
  - L1 total variation, Bray-Curtis distance
  - Precision, Recall, F1-score
  - True Positives, False Positives, False Negatives

- **`contigs_per_rank.tsv`** - Contig-level accuracy by taxonomic rank
  - Number of contigs, correct classifications, accuracy percentage

- **`contigs_exact.tsv`** - Exact TaxID matching statistics
  - Usable pairs, exact matches, accuracy percentage

#### Debug Information
- **`_debug_info.txt`** - Paths and configuration used for evaluation
- **`summary.txt`** - Human-readable summary of results

### Understanding the Results

#### Profile-Level Metrics
- **L1 Total Variation**: Lower is better (0-100%)
- **Bray-Curtis Distance**: Lower is better (0-100%)
- **Precision**: True positives / (True positives + False positives)
- **Recall**: True positives / (True positives + False negatives)
- **F1-Score**: Harmonic mean of precision and recall

#### Contig-Level Metrics
- **Exact TaxID Accuracy**: Percentage of contigs with exact taxonomic ID matches
- **Rank-wise Accuracy**: Classification accuracy at each taxonomic level

### Example Output Interpretation

```
# Profile-level metrics (per rank)
superkingdom    L1=50.000  BC=100.000%  P/R/F1=0.0/0.0/0.0% (TP=0, FP=0, FN=1)
phylum          L1=50.000  BC=100.000%  P/R/F1=0.0/0.0/0.0% (TP=0, FP=0, FN=3)
family          L1=61.158  BC=92.331%   P/R/F1=100.0/10.0/18.2% (TP=1, FP=0, FN=9)
species         L1=66.913  BC=98.163%   P/R/F1=33.3/19.0/24.2% (TP=4, FP=8, FN=17)

# Contig-level accuracy
Exact TaxID: 26/188149 (0.01%)
superkingdom    n=188149    acc=100.00%
phylum          n=188149    acc=99.79%
family          n=188149    acc=88.44%
species         n=188149    acc=0.07%
```

### Troubleshooting

#### Common Issues

1. **Missing CAMI Dataset**
   ```bash
   # Ensure CAMI dataset is properly downloaded and extracted
   ls /data/cami/sample_0/
   # Should contain: sample_0.fna, taxonomic_profile_0.txt, etc.
   ```

2. **Missing Tools**
   ```bash
   # Check if required tools exist
   ls /data/tools/hymet2cami_fixed.py
   ls /data/tools/build_id_map.py
   ls /data/tools/mini_classify.py
   ```

3. **Taxonomy Database Issues**
   ```bash
   # Ensure taxonomy database is configured
   ls /data/HYMET/taxonomy_files/
   # Should contain NCBI taxonomy files
   ```

4. **Memory Issues**
   ```bash
   # Reduce thread count for lower memory usage
   THREADS=4 bash benchmark_cami.sh
   ```

#### Performance Optimization

- **For Large Datasets**: Increase `THREADS` for faster processing
- **For Memory Constraints**: Reduce `THREADS` and ensure sufficient disk space
- **For Better Accuracy**: Ensure high-quality input sequences and proper taxonomy database

### Advanced Usage

#### Custom Evaluation Scripts

You can also run the evaluation components separately:

```bash
# Run only HYMET pipeline
bash run_hymet_cami.sh

# Run only evaluation
python3 tools/eval_cami.py \
  --pred-profile /path/to/hymet.sample_0.cami.tsv \
  --truth-profile /path/to/taxonomic_profile_0.txt \
  --pred-contigs /path/to/classified_sequences.tsv \
  --truth-contigs /path/to/gsa_mapping.tsv \
  --outdir /path/to/eval_results
```

#### Batch Processing

For processing multiple samples:

```bash
for sample in sample_0 sample_1 sample_2; do
  INPUT_FASTA="/data/cami/${sample}.fna" \
  OUTDIR="/data/results/${sample}" \
  bash benchmark_cami.sh
done
```

## Support

For questions or issues, please open an [issue](https://github.com/ieeta-pt/HYMET/issues) in the repository.
