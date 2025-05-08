#!/bin/bash
set -euo pipefail

# ==============================================
# 1. VALIDAÇÃO INICIAL E CONFIGURAÇÃO
# ==============================================

# Verificar número de argumentos
if [ "$#" -ne 8 ]; then
    echo "Usage: $0 <input_dir> <mash_screen> <screen_tab> <filtered_screen> <sorted_screen> <top_hits> <selected_genomes> <initial_threshold>"
    exit 1
fi

# Receber argumentos
INPUT_DIR="$1"
MASH_SCREEN="$2"
SCREEN_TAB="$3"
FILTERED_SCREEN="$4"
SORTED_SCREEN="$5"
TOP_HITS="$6"
SELECTED_GENOMES="$7"
INITIAL_THRESHOLD="$8"

# Criar diretório de output se não existir
mkdir -p "$(dirname "$SCREEN_TAB")" || {
    echo "ERROR: Cannot create output directory for $SCREEN_TAB"
    exit 1
}

# ==============================================
# 2. VERIFICAÇÃO DE ARQUIVOS DE ENTRADA
# ==============================================

# Verificar se o diretório de input existe
if [ ! -d "$INPUT_DIR" ]; then
    echo "ERROR: Input directory $INPUT_DIR does not exist"
    exit 1
fi

# Verificar se existem arquivos .fna no diretório
if [ -z "$(ls -A "$INPUT_DIR"/*.fna 2>/dev/null)" ]; then
    echo "ERROR: No .fna files found in $INPUT_DIR"
    exit 1
fi

# Verificar se o arquivo sketch existe
if [ ! -f "$MASH_SCREEN" ]; then
    echo "ERROR: Mash sketch file $MASH_SCREEN not found"
    exit 1
fi

# ==============================================
# 3. EXECUÇÃO DO MASH SCREEN
# ==============================================

echo "Running mash screen..."
if ! mash screen -p 8 -v 0.9 "$MASH_SCREEN" "$INPUT_DIR"/*.fna > "$SCREEN_TAB"; then
    echo "ERROR: mash screen command failed"
    exit 1
fi

# Verificar se gerou resultados
if [ ! -s "$SCREEN_TAB" ]; then
    echo "WARNING: mash screen produced empty results"
    touch "$SCREEN_TAB" "$FILTERED_SCREEN" "$SORTED_SCREEN" "$TOP_HITS" "$SELECTED_GENOMES"
    exit 0
fi

# ==============================================
# 4. PROCESSAMENTO DOS RESULTADOS
# ==============================================

echo "Processing results..."
sort -u -k5,5 "$SCREEN_TAB" > "$FILTERED_SCREEN"
sort -gr "$FILTERED_SCREEN" > "$SORTED_SCREEN"

# ==============================================
# 5. ANÁLISE DE THRESHOLD
# ==============================================

num_sequences=$(find "$INPUT_DIR" -maxdepth 1 -name "*.fna" | wc -l)
min_candidates=$(echo "$num_sequences * 3.25" | bc | awk '{printf("%d\n",$1 + 0.5)}')
min_candidates=$(( min_candidates < 5 ? 5 : min_candidates ))

best_threshold=$INITIAL_THRESHOLD
current_threshold=$INITIAL_THRESHOLD
threshold_found=0

echo "===================================="
echo "Number of input sequences: $num_sequences"
echo "Minimum expected candidates: $min_candidates"
echo "===================================="

while (( $(echo "$current_threshold >= 0.70" | bc -l) )); do
    count=$(awk -v t="$current_threshold" '$1 > t' "$SORTED_SCREEN" 2>/dev/null | wc -l)
    
    echo "Testing threshold: $current_threshold"
    echo "Candidates found: $count"

    if [ "$count" -ge "$min_candidates" ]; then
        best_threshold=$current_threshold
        threshold_found=1
        break
    fi
    
    current_threshold=$(echo "$current_threshold - 0.02" | bc -l)
done

if [ "$threshold_found" -eq 0 ]; then
    best_threshold=0.71
    count=$(awk -v t="$best_threshold" '$1 > t' "$SORTED_SCREEN" 2>/dev/null | wc -l)
    echo "No suitable threshold found. Using default 0.70."
fi

# ==============================================
# 6. GERAÇÃO DOS ARQUIVOS FINAIS
# ==============================================

awk -v threshold="$best_threshold" '$1 > threshold' "$SORTED_SCREEN" > "$TOP_HITS"
cut -f5 "$TOP_HITS" > "$SELECTED_GENOMES"

echo "===================================="
echo "Final threshold used: $best_threshold"
echo "Candidates found: $count"
echo "===================================="

exit 0