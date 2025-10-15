#!/bin/bash

# Accept arguments
INPUT_DIR="$1"          # This argument will no longer be used for .fna files directly.
REFERENCE_SET="$2"      # Path to combined_genomes.fasta.
NT_MMI="$3"             # Path to the minimap2 index (reference.mmi).
RESULTADOS_PAF="$4"     # Path to the alignment results (resultados.paf).

# Create the reference set index
echo "Creating index with minimap2..."
minimap2 -I2g -d "$NT_MMI" "$REFERENCE_SET" # General index

# Check if the index creation was successful
if [ $? -ne 0 ]; then
    echo "Error creating index with minimap2."
    exit 1
fi

# Run alignment using minimap2 (for long reads)
echo "Running alignment with minimap2..."
minimap2 -x asm10 "$NT_MMI" "$INPUT_DIR"/*.fna >"$RESULTADOS_PAF"

# Check if the alignment was successful
if [ $? -ne 0 ]; then
    echo "Error running alignment with minimap2."
    exit 1
fi

echo "Alignment completed successfully! Results saved to $RESULTADOS_PAF."