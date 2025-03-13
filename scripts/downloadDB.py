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

# Configurações globais
MAX_WORKERS = 64  # Número máximo de threads para downloads paralelos
RETRIES = 3       # Número máximo de tentativas por download
TIMEOUT = 15      # Tempo limite (em segundos) para cada tentativa de download

# Configuração de logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def configurar_diretorios(output_dir, cache_dir):
    """
    Cria os diretórios de saída e cache, se não existirem.
    """
    os.makedirs(output_dir, exist_ok=True)
    os.makedirs(cache_dir, exist_ok=True)

class GenomeDownloader:
    def __init__(self, output_dir, cache_dir):
        self.output_dir = output_dir
        self.cache_dir = cache_dir
        self.assembly_summaries = self.baixar_assembly_summaries()
        self.failed_downloads = set()
        self.successful_downloads = set()
        self.assembly_data = self.carregar_assembly_summaries()

    def baixar_assembly_summaries(self):
        """
        Baixa os arquivos de sumário de assembly do NCBI (RefSeq e GenBank).
        """
        summaries = {}
        urls = [
            ("https://ftp.ncbi.nlm.nih.gov/genomes/refseq/assembly_summary_refseq.txt", "refseq"),
            ("https://ftp.ncbi.nlm.nih.gov/genomes/genbank/assembly_summary_genbank.txt", "genbank")
        ]
        
        for url, key in urls:
            cache_file = os.path.join(self.cache_dir, f"assembly_summary_{key}.txt")
            if not os.path.exists(cache_file):
                self.baixar_arquivo_wget(url, cache_file)
            summaries[key] = cache_file
        
        return summaries

    def baixar_arquivo_wget(self, url, destino, retries=RETRIES):
        """
        Baixa um arquivo usando wget.
        """
        for attempt in range(retries):
            try:
                comando = [
                    "wget",
                    "-O", destino,
                    "-q",  # Modo silencioso
                    "--tries=3",
                    "--timeout=15",
                    url
                ]
                subprocess.run(comando, check=True)
                return
            except subprocess.CalledProcessError as e:
                logging.warning(f"Tentativa {attempt + 1} de {retries} falhou para {url}: {e}")
                if attempt < retries - 1:
                    sleep(2 ** attempt)  # Exponential backoff
                else:
                    raise

    def carregar_assembly_summaries(self):
        """
        Carrega os dados dos arquivos de sumário de assembly.
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

    def processar_identificadores(self, genomes_file):
        """
        Processa o arquivo de identificadores de genomas.
        """
        with open(genomes_file) as f:
            return [self.extrair_gcf(line.strip()) for line in f if line.strip()]

    def extrair_gcf(self, filename):
        """
        Extrai o identificador GCF de um nome de arquivo.
        """
        parts = filename.split('_')
        return f"{parts[0]}_{parts[1]}"

    def executar_downloads(self, identifiers):
        """
        Executa os downloads em paralelo usando ThreadPoolExecutor.
        """
        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {executor.submit(self.baixar_genoma, gcf): gcf for gcf in identifiers}
            for future in as_completed(futures):
                gcf = futures[future]
                try:
                    result = future.result()
                    if not result:
                        self.failed_downloads.add(gcf)
                except Exception as e:
                    logging.error(f"Erro no download de {gcf}: {str(e)}")
                    self.failed_downloads.add(gcf)

    def baixar_genoma(self, gcf):
        """
        Baixa um genoma usando wget.
        """
        metadata = self.assembly_data.get(gcf)
        if not metadata:
            return False

        output_path = os.path.join(self.output_dir, metadata['file_name'])
        if os.path.exists(output_path):
            self.successful_downloads.add(metadata['file_name'])
            return True

        return self.salvar_genoma(metadata, output_path)

    def salvar_genoma(self, metadata, output_path, retries=RETRIES):
        """
        Baixa e salva um genoma usando wget.
        """
        url = f"{metadata['ftp_path']}/{os.path.basename(metadata['ftp_path'])}_genomic.fna.gz"
        temp_gz = f"{output_path}.gz"

        for attempt in range(retries):
            try:
                comando = [
                    "wget",
                    "-O", temp_gz,
                    "-q",
                    "--tries=3",
                    "--timeout=15",
                    url
                ]
                subprocess.run(comando, check=True)

                with gzip.open(temp_gz, 'rb') as f_in:
                    with open(output_path, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)

                os.remove(temp_gz)
                self.successful_downloads.add(metadata['file_name'])
                return True
            except subprocess.CalledProcessError as e:
                logging.warning(f"Tentativa {attempt + 1} de {retries} falhou para {metadata['file_name']}: {e}")
                if attempt < retries - 1:
                    sleep(2 ** attempt)
                else:
                    logging.error(f"Falha no download de {metadata['file_name']}: {str(e)}")
                    return False

    def create_detailed_taxonomy_from_directory(self, taxonomy_file):
        """
        Cria um arquivo de taxonomia detalhada a partir dos genomas baixados.
        """
        mapping = defaultdict(lambda: {"taxid": "Unknown TaxID", "identifiers": set()})
        
        for file in os.listdir(self.output_dir):
            if file.endswith(".fna"):
                file_path = os.path.join(self.output_dir, file)
                gcf = self.extrair_gcf(file)
                with open(file_path, "r") as f:
                    for line in f:
                        if line.startswith(">"):
                            identifier = line.split()[0][1:]  # Remove o '>' e pega o identificador
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

        logging.info(f"Arquivo de taxonomia detalhada salvo em: {taxonomy_file}")

    def concatenar_genomas(self, output_file):
        """
        Concatena todos os genomas baixados em um único arquivo.
        """
        logging.info("Concatenando genomas...")
        with open(output_file, 'w') as out_f:
            for filename in self.successful_downloads:
                file_path = os.path.join(self.output_dir, filename)
                try:
                    with open(file_path, 'r') as in_f:
                        shutil.copyfileobj(in_f, out_f)
                except FileNotFoundError:
                    logging.warning(f"Aviso: Arquivo {filename} não encontrado, pulando...")
        logging.info(f"Genomas concatenados em {output_file}")

if __name__ == "__main__":
    if len(sys.argv) != 5:
        print("Uso: python3 download_genomes.py <genomes_file> <output_dir> <taxonomy_file> <cache_dir>")
        sys.exit(1)

    genomes_file = sys.argv[1]
    output_dir = sys.argv[2]
    taxonomy_file = sys.argv[3]
    cache_dir = sys.argv[4]

    configurar_diretorios(output_dir, cache_dir)
    
    downloader = GenomeDownloader(output_dir, cache_dir)
    identifiers = downloader.processar_identificadores(genomes_file)
    
    logging.info(f"Iniciando download de {len(identifiers)} genomas...")
    downloader.executar_downloads(identifiers)
    downloader.create_detailed_taxonomy_from_directory(taxonomy_file)
    
    combined_genomes_file = os.path.join(output_dir, "combined_genomes.fasta")
    downloader.concatenar_genomas(combined_genomes_file)
    
    logging.info("\nResumo:")
    logging.info(f" - Baixados com sucesso: {len(downloader.successful_downloads)}")
    logging.info(f" - Falhas: {len(downloader.failed_downloads)}")
    logging.info(f" - Arquivo combinado: {combined_genomes_file}")
