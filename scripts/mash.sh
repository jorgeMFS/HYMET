#!/bin/bash

# Recebe argumentos
INPUT_DIR="$1"
MASH_SCREEN="$2"
SCREEN_TAB="$3"
FILTERED_SCREEN="$4"
SORTED_SCREEN="$5"
TOP_HITS="$6"
SELECTED_GENOMES="$7"
INITIAL_THRESHOLD="$8"

# Passo 1: Executa o mash screen e gera o sorted_screen
mash screen -p 8 -v 0.9 "$MASH_SCREEN" "$INPUT_DIR"/*.fna > "$SCREEN_TAB"
sort -u -k5,5 "$SCREEN_TAB" > "$FILTERED_SCREEN"
sort -gr "$FILTERED_SCREEN" > "$SORTED_SCREEN"

# Passo 2: Ajusta o threshold e seleciona os genomas
num_sequences=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fna" | wc -l)
min_candidates=$(echo "$num_sequences * 3.25" | bc | awk '{printf("%d\n",$1 + 0.5)}')
min_candidates=$(( min_candidates < 5 ? 5 : min_candidates ))

best_threshold=$INITIAL_THRESHOLD
current_threshold=$INITIAL_THRESHOLD
threshold_found=0

echo "===================================="
echo "Número de sequências de entrada: $num_sequences"
echo "Número mínimo de candidatos esperado: $min_candidates"
echo "===================================="

while (( $(echo "$current_threshold >= 0.70" | bc -l) )); do
    count=$(awk -v t="$current_threshold" '$1 > t' "$SORTED_SCREEN" | wc -l)
    
    echo "Testando threshold: $current_threshold"
    echo "Candidatos encontrados: $count"

    if [ "$count" -ge "$min_candidates" ]; then
        best_threshold=$current_threshold
        threshold_found=1
        break
    fi
    
    current_threshold=$(echo "$current_threshold - 0.02" | bc -l)
done

if [ "$threshold_found" -eq 0 ]; then
    best_threshold=0.71
    count=$(awk -v t="$best_threshold" '$1 > t' "$SORTED_SCREEN" | wc -l)
    echo "Nenhum threshold encontrado. Usando 0.70."
fi

# Filtra com o melhor threshold encontrado
awk -v threshold="$best_threshold" '$1 > threshold' "$SORTED_SCREEN" > "$TOP_HITS"
cut -f5 "$TOP_HITS" > "$SELECTED_GENOMES"

echo "===================================="
echo "Threshold final utilizado: $best_threshold"
echo "Candidatos encontrados: $count"
echo "===================================="
