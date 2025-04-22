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
git clone https://github.com/inesbmartins02/HYMET.git
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

## Support

For questions or issues, please open an [issue](https://github.com/inesbmartins02/HYMET/issues) in the repository.
