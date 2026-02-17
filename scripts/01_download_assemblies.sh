#!/usr/bin/env bash
set -euo pipefail

# 01_download_assemblies.sh
# Downloads genome FASTA files for all GCA/GCF accessions in lists/assembly_list.txt

PROJECT="$HOME/kp_liver_project"
LIST="$PROJECT/lists/assembly_list.txt"
DL="$PROJECT/downloads"
OUT="$PROJECT/genomes_fasta"
LOG="$PROJECT/logs/01_download_assemblies.log"

mkdir -p "$DL" "$OUT" "$PROJECT/logs"

if ! command -v datasets >/dev/null 2>&1; then
  echo "ERROR: NCBI datasets CLI not found in PATH" | tee -a "$LOG"
  exit 1
fi

if [[ ! -s "$LIST" ]]; then
  echo "ERROR: Missing or empty list: $LIST" | tee -a "$LOG"
  exit 1
fi

echo "Downloading assemblies from: $LIST" | tee "$LOG"

datasets download genome accession \
  --inputfile "$LIST" \
  --include genome \
  --filename "$DL/gca_download.zip" \
  2>&1 | tee -a "$LOG"

unzip -o "$DL/gca_download.zip" -d "$DL/gca_dataset" 2>&1 | tee -a "$LOG"

find "$DL/gca_dataset" -name "*genomic.fna" -exec cp {} "$OUT/" \;

echo "Assembly FASTA count:" $(ls "$OUT"/*.fna 2>/dev/null | wc -l) | tee -a "$LOG"
echo "DONE" | tee -a "$LOG"
