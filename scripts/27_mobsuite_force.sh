#!/bin/bash
#SBATCH --job-name=mobsuite_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=08:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mobsuite_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mobsuite_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate mobsuite_env

FNADIR=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES
OUTDIR=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/MOBsuite
DBDIR=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mob_db

mkdir -p "$OUTDIR"

echo "Starting MOB-suite at $(date)"

for fna in "$FNADIR"/*.fna; do
    sample=$(basename "$fna" .fna)
    echo "Processing $sample"
    mob_recon \
        --infile "$fna" \
        --outdir "$OUTDIR/${sample}" \
        --database_directory "$DBDIR" \
        --num_threads 6 --force
done

echo "Finished MOB-suite at $(date)"
