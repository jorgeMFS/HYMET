#!/usr/bin/env python3
import os
import sys
import gzip
import shutil
import csv
import logging
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from collections import defaultdict
from time import sleep

# Global settings
MAX_WORKERS = 64  # Maximum number of threads for parallel downloads
RETRIES = 3       # Maximum number of retries per download
TIMEOUT = 15      # Timeout (in seconds) for each download attempt

# Logging configuration
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def setup_directories(output_dir, cache_dir):
    """
    Creates the output and cache directories if they don't exist.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

class GenomeDownloader:
    def __init__(self, output_dir, cache_dir):
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.assembly_summaries = self.download_assembly_summaries()
        self.failed_downloads = set()
        self.successful_downloads = set()
        self.assembly_data = self.load_assembly_summaries()

    def download_assembly_summaries(self):
        """
        Downloads the assembly summary files from NCBI (RefSeq and GenBank).
        """
        summaries = {}
        urls = [
            ("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", "refseq"),
            ("https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt", "genbank")
        ]
        
        for url, key in urls:
            cache_file = os.path.join(self.cache_dir, f"assembly_summary_{key}.txt")
            if not os.path.exists(cache_file):
                self.download_file_wget(url, cache_file)
            summaries[key] = cache_file
        
        return summaries

    def download_file_wget(self, url, destination, retries=RETRIES):
        """
        Downloads a file using wget.
        """
        for attempt in range(retries):
            try:
                command = [
                    "wget",
                    "-O", destination,
                    "-q",  # Silent mode
                    "--tries=3",
                    "--timeout=15",
                    url
                ]
                subprocess.run(command, check=True)
                return
            except subprocess.CalledProcessError as e:
                logging.warning(f"Attempt {attempt + 1} of {retries} failed for {url}: {e}")
                if attempt < retries - 1:
                    sleep(2 ** attempt)  # Exponential backoff
                else:
                    raise

    def load_assembly_summaries(self):
        """
        Loads data from the assembly summary files.
        """
        assembly_data = {}
        for key, file_path in self.assembly_summaries.items():
            with open(file_path, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) > 19 and parts[19]:
                        gcf = parts[0]
                        assembly_data[gcf] = {
                            'ftp_path': parts[19].replace('ftp://', 'https://'),
                            'organism_name': parts[7],
                            'taxid': parts[5],
                            'file_name': f"{parts[0]}_{parts[1]}.fna"
                        }
        return assembly_data

    def process_identifiers(self, genomes_file):
        """
        Processes the file containing genome identifiers.
        """
        with open(genomes_file) as f:
            return [self.extract_gcf(line.strip()) for line in f if line.strip()]

    def extract_gcf(self, filename):
        """
        Extracts the GCF identifier from a filename.
        """
        parts = filename.split('_')
        return f"{parts[0]}_{parts[1]}"

    def execute_downloads(self, identifiers):
        """
        Executes downloads in parallel using ThreadPoolExecutor.
        """
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(self.download_genome, gcf): gcf for gcf in identifiers}
            for future in as_completed(futures):
                gcf = futures[future]
                try:
                    result = future.result()
                    if not result:
                        self.failed_downloads.add(gcf)
                except Exception as e:
                    logging.error(f"Error downloading {gcf}: {str(e)}")
                    self.failed_downloads.add(gcf)

    def download_genome(self, gcf):
        """
        Downloads a genome using wget.
        """
        metadata = self.assembly_data.get(gcf)
        if not metadata:
            return False

        output_path = os.path.join(self.output_dir, metadata['file_name'])
        if os.path.exists(output_path):
            self.successful_downloads.add(metadata['file_name'])
            return True

        return self.save_genome(metadata, output_path)

    def save_genome(self, metadata, output_path, retries=RETRIES):
        """
        Downloads and saves a genome using wget.
        """
        url = f"{metadata['ftp_path']}/{os.path.basename(metadata['ftp_path'])}_genomic.fna.gz"
        temp_gz = f"{output_path}.gz"

        for attempt in range(retries):
            try:
                command = [
                    "wget",
                    "-O", temp_gz,
                    "-q",
                    "--tries=3",
                    "--timeout=15",
                    url
                ]
                subprocess.run(command, check=True)

                with gzip.open(temp_gz, 'rb') as f_in:
                    with open(output_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(temp_gz)
                self.successful_downloads.add(metadata['file_name'])
                return True
            except subprocess.CalledProcessError as e:
                logging.warning(f"Attempt {attempt + 1} of {retries} failed for {metadata['file_name']}: {e}")
                if attempt < retries - 1:
                    sleep(2 ** attempt)
                else:
                    logging.error(f"Failed to download {metadata['file_name']}: {str(e)}")
                    return False

    def create_detailed_taxonomy_from_directory(self, taxonomy_file):
        """
        Creates a detailed taxonomy file from the downloaded genomes.
        """
        mapping = defaultdict(lambda: {"taxid": "Unknown TaxID", "identifiers": set()})
        
        for file in os.listdir(self.output_dir):
            if file.endswith(".fna"):
                file_path = os.path.join(self.output_dir, file)
                gcf = self.extract_gcf(file)
                with open(file_path, "r") as f:
                    for line in f:
                        if line.startswith(">"):
                            identifier = line.split()[0][1:]  # Remove '>' and get the identifier
                            mapping[gcf]["identifiers"].add(identifier)
                
                metadata = self.assembly_data.get(gcf, {})
                mapping[gcf]["taxid"] = metadata.get('taxid', "Unknown TaxID")
        
        with open(taxonomy_file, "w", newline="") as csvfile:
            writer = csv.writer(csvfile, delimiter="\t")
            writer.writerow(["GCF", "TaxID", "Identifiers"])
            for gcf, data in mapping.items():
                writer.writerow([
                    gcf,
                    data["taxid"],
                    ";".join(data["identifiers"])
                ])

        logging.info(f"Detailed taxonomy file saved at: {taxonomy_file}")

    def concatenate_genomes(self, output_file):
        """
        Concatenates all downloaded genomes into a single file.
        """
        logging.info("Concatenating genomes...")
        with open(output_file, 'w') as out_f:
            for filename in self.successful_downloads:
                file_path = os.path.join(self.output_dir, filename)
                try:
                    with open(file_path, 'r') as in_f:
                        shutil.copyfileobj(in_f, out_f)
                except FileNotFoundError:
                    logging.warning(f"Warning: File {filename} not found, skipping...")
        logging.info(f"Genomes concatenated into {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Usage: python3 download_genomes.py <genomes_file> <output_dir> <taxonomy_file> <cache_dir>")
        sys.exit(1)

    genomes_file = sys.argv[1]
    output_dir = sys.argv[2]
    taxonomy_file = sys.argv[3]
    cache_dir = sys.argv[4]

    setup_directories(output_dir, cache_dir)
    
    downloader = GenomeDownloader(output_dir, cache_dir)
    identifiers = downloader.process_identifiers(genomes_file)
    
    logging.info(f"Starting download of {len(identifiers)} genomes...")
    downloader.execute_downloads(identifiers)
    downloader.create_detailed_taxonomy_from_directory(taxonomy_file)
    
    combined_genomes_file = os.path.join(output_dir, "combined_genomes.fasta")
    downloader.concatenate_genomes(combined_genomes_file)
    
    logging.info("\nSummary:")
    logging.info(f" - Successfully downloaded: {len(downloader.successful_downloads)}")
    logging.info(f" - Failed downloads: {len(downloader.failed_downloads)}")
    logging.info(f" - Combined file: {combined_genomes_file}")