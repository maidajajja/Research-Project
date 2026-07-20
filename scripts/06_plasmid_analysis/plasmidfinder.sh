#!/bin/bash
#SBATCH --job-name=plasmidfinder_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/plasmidfinder_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/plasmidfinder_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate mobsuite_env

FNADIR=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES
OUTDIR=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/PlasmidFinder

mkdir -p "$OUTDIR"

echo "Starting PlasmidFinder at $(date)"

for fna in "$FNADIR"/*.fna; do
    sample=$(basename "$fna" .fna)
    mkdir -p "$OUTDIR/${sample}"
    plasmidfinder.py \
        -i "$fna" \
        -o "$OUTDIR/${sample}" \
        -p /users/k22017808/.conda/envs/mobsuite_env/share/plasmidfinder-2.1.6/database \
        -l 0.60 \
        -t 0.95
done

echo "Finished PlasmidFinder at $(date)"
