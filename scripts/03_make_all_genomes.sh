#!/usr/bin/env bash
set -euo pipefail

# 03_make_all_genomes.sh
# Collects assembly FASTAs + GenBank isolate FASTAs into genomes_fasta/all_genomes

PROJECT="$HOME/kp_liver_project"
OUT="$PROJECT/genomes_fasta/all_genomes"
LOG="$PROJECT/logs/03_make_all_genomes.log"

mkdir -p "$OUT" "$PROJECT/logs"

echo "Collecting genomes into $OUT" | tee "$LOG"

cp "$PROJECT/genomes_fasta/"*.fna "$OUT/" 2>/dev/null || true
cp "$PROJECT/genomes_fasta/genbank_isolates/"*.fasta "$OUT/" 2>/dev/null || true
cp "$PROJECT/assemblies/"*.fasta "$OUT/" 2>/dev/null || true

echo "Total genomes in all_genomes:" $(ls "$OUT" | wc -l) | tee -a "$LOG"
echo "DONE" | tee -a "$LOG"
