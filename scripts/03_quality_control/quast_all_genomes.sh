#!/bin/bash
#SBATCH --job-name=quast_all
#SBATCH --output=logs/quast_all_%j.log
#SBATCH --error=logs/quast_all_%j.err
#SBATCH --time=4:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu

#=============================================================================
# Run QUAST on all 234 K. pneumoniae genomes
# Author: Maida Jajja
# Date: 2025-02-19
#=============================================================================

set -e

# Activate conda environment
eval "$(conda shell.bash hook)"
conda activate quast_env

# Project directories
PROJECT_DIR="$HOME/kp_liver_project"
GENOMES_DIR="${PROJECT_DIR}/genomes_fasta/all_genomes"
QUAST_OUTPUT="${PROJECT_DIR}/quast_results"

echo "========================================"
echo "Running QUAST on all genomes"
echo "Started: $(date)"
echo "========================================"

# Run QUAST on all genomes at once
quast.py ${GENOMES_DIR}/*.fna ${GENOMES_DIR}/*.fasta \
  -o "$QUAST_OUTPUT" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --silent

echo ""
echo "========================================"
echo "QUAST Complete!"
echo "========================================"
echo "Results: $QUAST_OUTPUT"
echo "HTML report: $QUAST_OUTPUT/report.html"
echo "Text report: $QUAST_OUTPUT/report.txt"
echo "TSV report: $QUAST_OUTPUT/report.tsv"
echo "Completed: $(date)"
echo "========================================"
