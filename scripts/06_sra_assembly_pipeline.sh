#!/bin/bash
#SBATCH --job-name=sra_assembly
#SBATCH --output=logs/sra_assembly_%A_%a.log
#SBATCH --error=logs/sra_assembly_%A_%a.err
#SBATCH --array=1-10
#SBATCH --time=24:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=8
#SBATCH --partition=cpu

#=============================================================================
# Download, QC, and assemble K. pneumoniae genomes from SRA
# Author: Maida Jajja
# Date: 2025-02-17
#=============================================================================

set -e  # Exit on error

# Activate conda environment
source ~/.bashrc
conda activate assembly_env

# Load modules
module load sra-tools/3.0.3-gcc-13.2.0

# Project directories
PROJECT_DIR="$HOME/kp_liver_project"
RAW_READS_DIR="${PROJECT_DIR}/raw_reads"
FASTQC_DIR="${PROJECT_DIR}/fastqc"
ASSEMBLIES_DIR="${PROJECT_DIR}/assemblies"
GENOMES_DIR="${PROJECT_DIR}/genomes_fasta/all_genomes"
LOGS_DIR="${PROJECT_DIR}/logs"

# Create directories
mkdir -p "$RAW_READS_DIR" "$FASTQC_DIR" "$ASSEMBLIES_DIR" "$GENOMES_DIR" "$LOGS_DIR"

# SRA accessions
SRA_LIST=(
  "SRR16202849"
  "SRR16202839"
  "SRR16202829"
  "SRR26938680"
  "SRR26896786"
  "SRR27842369"
  "SRR31742341"
  "SRR31742346"
  "SRR28905830"
  "SRR28905831"
)

# Get SRA accession for this array task
SRA_ACC="${SRA_LIST[$SLURM_ARRAY_TASK_ID-1]}"

echo "========================================"
echo "Processing: $SRA_ACC"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Started: $(date)"
echo "========================================"

#=============================================================================
# Step 1: Download FASTQ from SRA
#=============================================================================

echo "[1/4] Downloading FASTQ files from SRA..."
cd "$RAW_READS_DIR"

# Download using fasterq-dump (faster than fastq-dump)
fasterq-dump "$SRA_ACC" \
  --split-files \
  --threads "$SLURM_CPUS_PER_TASK" \
  --progress

# Check if files were created
if [ ! -f "${SRA_ACC}_1.fastq" ]; then
  echo "ERROR: Download failed for $SRA_ACC"
  exit 1
fi

echo "  ✓ Downloaded ${SRA_ACC}_1.fastq and ${SRA_ACC}_2.fastq"

#=============================================================================
# Step 2: Quality Control with fastp
#=============================================================================

echo "[2/4] Running quality control with fastp..."

QC_DIR="${RAW_READS_DIR}/${SRA_ACC}_qc"
mkdir -p "$QC_DIR"

fastp \
  -i "${SRA_ACC}_1.fastq" \
  -I "${SRA_ACC}_2.fastq" \
  -o "${QC_DIR}/${SRA_ACC}_1_clean.fastq.gz" \
  -O "${QC_DIR}/${SRA_ACC}_2_clean.fastq.gz" \
  --thread "$SLURM_CPUS_PER_TASK" \
  --detect_adapter_for_pe \
  --cut_front \
  --cut_tail \
  --cut_window_size 4 \
  --cut_mean_quality 20 \
  --qualified_quality_phred 20 \
  --unqualified_percent_limit 40 \
  --n_base_limit 5 \
  --length_required 50 \
  --html "${FASTQC_DIR}/${SRA_ACC}_fastp.html" \
  --json "${FASTQC_DIR}/${SRA_ACC}_fastp.json" \
  2>&1 | tee "${LOGS_DIR}/${SRA_ACC}_fastp.log"

echo "  ✓ Quality control complete"
echo "  ✓ Clean reads: ${QC_DIR}/${SRA_ACC}_1_clean.fastq.gz"

# Remove raw FASTQ to save space
rm "${SRA_ACC}_1.fastq" "${SRA_ACC}_2.fastq"

#=============================================================================
# Step 3: Genome Assembly with SPAdes
#=============================================================================

echo "[3/4] Assembling genome with SPAdes..."

ASSEMBLY_DIR="${ASSEMBLIES_DIR}/${SRA_ACC}"

spades.py \
  -1 "${QC_DIR}/${SRA_ACC}_1_clean.fastq.gz" \
  -2 "${QC_DIR}/${SRA_ACC}_2_clean.fastq.gz" \
  -o "$ASSEMBLY_DIR" \
  --threads "$SLURM_CPUS_PER_TASK" \
  --memory $((SLURM_MEM_PER_NODE / 1024)) \
  --careful \
  2>&1 | tee "${LOGS_DIR}/${SRA_ACC}_spades.log"

if [ ! -f "${ASSEMBLY_DIR}/scaffolds.fasta" ]; then
  echo "ERROR: Assembly failed for $SRA_ACC"
  exit 1
fi

echo "  ✓ Assembly complete: ${ASSEMBLY_DIR}/scaffolds.fasta"

#=============================================================================
# Step 4: Copy final assembly to genomes folder
#=============================================================================

echo "[4/4] Copying final assembly..."

cp "${ASSEMBLY_DIR}/scaffolds.fasta" "${GENOMES_DIR}/${SRA_ACC}_assembled.fasta"

echo "  ✓ Final genome: ${GENOMES_DIR}/${SRA_ACC}_assembled.fasta"

#=============================================================================
# Summary
#=============================================================================

echo ""
echo "========================================"
echo "Assembly Complete: $SRA_ACC"
echo "========================================"
echo "Final assembly: ${GENOMES_DIR}/${SRA_ACC}_assembled.fasta"
echo "Assembly stats: ${ASSEMBLY_DIR}/scaffolds.fasta"
echo "QC report: ${FASTQC_DIR}/${SRA_ACC}_fastp.html"
echo "Completed: $(date)"
echo "========================================"
