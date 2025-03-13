#!/bin/bash

# Receives arguments
INPUT_DIR="$1"
MASH_SCREEN="$2"
SCREEN_TAB="$3"
FILTERED_SCREEN="$4"
SORTED_SCREEN="$5"
TOP_HITS="$6"
SELECTED_GENOMES="$7"
THRESHOLD="$8"  # New argument for the threshold

# Runs mash screen
mash screen -p 8 -w -v 0.9 "$MASH_SCREEN" "$INPUT_DIR"/*.fna > "$SCREEN_TAB"

# Filters unique entries by column 5
sort -u -k5,5 "$SCREEN_TAB" > "$FILTERED_SCREEN"

# Sorts in descending order
sort -gr "$FILTERED_SCREEN" > "$SORTED_SCREEN"

# Uses the user-defined threshold
awk -v threshold="$THRESHOLD" '$1 > threshold' "$SORTED_SCREEN" > "$TOP_HITS"

# Extracts column 5 from the best entries
cut -f5 "$TOP_HITS" > "$SELECTED_GENOMES"

echo "Process completed successfully!"