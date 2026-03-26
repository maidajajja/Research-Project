#!/bin/bash
#SBATCH --job-name=mlst_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=22G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/MLST/mlst_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/MLST/mlst_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate mlst_env

GENOMES=/scratch/users/k22017808/KP_Research_Project/03_FinalGenomes/genomes_fasta/all_genomes
OUTDIR=/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/MLST

mkdir -p "$OUTDIR"

echo "Starting MLST at $(date)"

mlst \
    --scheme klebsiella \
    --threads 6 \
    --nopath \
    "$GENOMES"/*.fasta \
    > "$OUTDIR/mlst_results.tsv"

echo "Finished MLST at $(date)"

awk -F'\t' '{print $3}' "$OUTDIR/mlst_results.tsv" | sort | uniq -c | sort -rn
