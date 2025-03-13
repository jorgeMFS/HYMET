#!/bin/bash

# Parameters
INPUT_QUERY="$1"       # Query file (concatenated_input.fasta).
REFERENCE_SET="$2"     # Reference file (combined_genomes.fasta).
OUTPUT_FILE="$3"       # MashMap output file (mashmap.out).
THREADS=8              # Number of threads (adjust as needed).

# Check if arguments were provided
if [ -z "$INPUT_QUERY" ] || [ -z "$REFERENCE_SET" ] || [ -z "$OUTPUT_FILE" ]; then
    echo "Usage: $0 <query_file> <reference_file> <output_file>"
    exit 1
fi

# Check if MashMap is installed
if ! command -v mashmap &> /dev/null; then
    echo "Error: MashMap not found. Please install MashMap and try again."
    exit 1
fi

# Run MashMap
echo "Running MashMap with $THREADS threads..."
mashmap -t "$THREADS" -r "$REFERENCE_SET" -q "$INPUT_QUERY" --perc_identity 85 -s 10000 -o "$OUTPUT_FILE"

# Check if MashMap ran successfully
if [ $? -ne 0 ]; then
    echo "Error running MashMap."
    exit 1
fi

echo "Alignment completed successfully! Results saved in $OUTPUT_FILE."
