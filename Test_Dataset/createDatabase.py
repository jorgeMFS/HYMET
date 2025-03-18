from os import listdir, path, makedirs, remove
from os.path import isfile, join
import random
import shutil
import wget
import requests
from bs4 import BeautifulSoup
from pathlib import Path
import gzip
import logging
import time

# Paths for input and output directories
# Note: These paths should be adjusted according to your local setup
summariesPath = Path(input("Enter the path to the assembly summaries directory: ")).resolve()
dstPathOfDatabaseSequences = Path(input("Enter the path to the destination directory for database sequences: ")).resolve()

# Seed for randomness
random.seed(0)

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def main():
    selected_domain = input("Enter the domain you want to process (e.g., archaea, bacteria, fungi): ").strip().lower()
    selectedSequences = getSequences(selected_domain)
    
    if not selectedSequences:
        logging.error(f"No sequences found for the domain: {selected_domain}")
        return
    
    _initialize()
    downloadSequences(selectedSequences["from_url"], dstPathOfDatabaseSequences)
    splitSequences(selectedSequences, dstPathOfDatabaseSequences)

def _initialize():
    if not path.exists(dstPathOfDatabaseSequences):
        makedirs(dstPathOfDatabaseSequences)

def splitSequences(selectedSequences, dstPathOfDatabaseSequences):
    logging.info("Starting to split and decompress sequences")
    allFiles = [f for f in listdir(dstPathOfDatabaseSequences) if isfile(join(dstPathOfDatabaseSequences, f))]
    
    for x in selectedSequences:
        for domain in selectedSequences[x]:
            dst = join(dstPathOfDatabaseSequences, domain)
            if not path.exists(dst):
                makedirs(dst)
            for entry in selectedSequences[x][domain]:
                fileName = entry.split("/")[-1]
                matchingFiles = [f for f in allFiles if fileName in f]
                if matchingFiles:
                    fileToMove = matchingFiles[0]
                    source = join(dstPathOfDatabaseSequences, fileToMove)
                    destination = join(dst, fileToMove)
                    shutil.move(source, destination)
                    logging.info(f"Moved {fileToMove} to {domain} folder")
                    
                    if destination.endswith('.gz'):
                        with gzip.open(destination, 'rb') as f_in:
                            with open(destination[:-3], 'wb') as f_out:
                                shutil.copyfileobj(f_in, f_out)
                        remove(destination)  
                        logging.info(f"Decompressed {fileToMove}")
                else:
                    logging.warning(f"File not found for entry: {entry}")
    
    logging.info("Finished splitting and decompressing sequences")

def getSequences(selected_domain):
    # Note: Assembly files need to be downloaded manually
    selectedSequences = {
        "from_url": {},
        "from_db": {}
    }
    
    onlyfiles = [f for f in listdir(summariesPath) if isfile(join(summariesPath, f)) and f.endswith("assembly_summary.txt")]
    
    for fileName in onlyfiles:
        domain = fileName.split("_")[0].lower() 
        if domain != selected_domain:
            continue 
        
        with open(join(summariesPath, fileName), 'r') as fp:
            content = fp.readlines()
            lenOfFile = len(content) - 1  
            numberToDownload = max(1, lenOfFile // 10)  # Calculate 10% of total entries
            
            listOfValues = sorted(random.sample(range(1, lenOfFile + 1), numberToDownload))  
            selectedSequences["from_url"][domain] = _getRandomEntries(listOfValues, content)

    return selectedSequences

def _getRandomEntries(listOfValues, content, singleColumn=False):
    readsOfInterest = []
    for n in listOfValues:
        if singleColumn:
            readsOfInterest.append(content[n].replace("\n", ""))
        else:
            readsOfInterest.append(content[n].split("\t")[19]) 
    
    return readsOfInterest

def downloadSequences(selectedSequences, dst):
    for domain in selectedSequences:
        for entry in selectedSequences[domain]:
            try:
                r = requests.get(entry, timeout=10)
                r.raise_for_status()
                soup = BeautifulSoup(r.text, 'html.parser')
                
                for link in soup.find_all('a'):
                    if "genomic.fna.gz" in link.get('href') and "from_genomic" not in link.get('href'):
                        fLink = join(entry, link.get('href'))
                        logging.info(f"Downloading {fLink} to {dst}")

                        success = False
                        attempts = 0
                        while not success and attempts < 3: 
                            try:
                                wget.download(fLink, out=str(dst))
                                success = True
                            except Exception as e:
                                attempts += 1
                                logging.error(f"Error downloading {fLink}: {e}")
                                time.sleep(5)  

                        if not success:
                            logging.error(f"Failed to download after multiple attempts: {fLink}")
            
            except requests.RequestException as e:
                logging.error(f"Error retrieving links from {entry}: {e}")

if __name__ == "__main__":
    main()
