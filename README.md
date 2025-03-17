# HYMET (Hybrid Metagenomic Tool)


## Installation and Configuration

Follow the steps below to install and configure **HYMET**

### 1. Clone the Repository

First, clone the repository to your local environment:

```bash
git clone [[https://github.com/your-username/your-repository.git](https://github.com/inesbmartins02/HYMET.git)]
cd HYMET
```



### 2. Installation with Docker (Recommended)

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



### 3. Installation with Conda

If you prefer using Conda, follow these steps:

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

After installation and configuration, first you should run the configuration script to download and prepare the taxonomy files, and define the main paths.

```bash
./config.pl
```

Then, you can run the main tool to perform taxonomic identification:

```bash
./main.pl
```

---

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


## 8. Support

For questions or issues, please open an [issue]([https://github.com/your-username/your-repository/issues](https://github.com/inesbmartins02/HYMET.git)) in the repository.
