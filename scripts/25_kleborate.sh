#!/bin/bash
#SBATCH --job-name=kleborate_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/09_Kleborate/kleborate_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/09_Kleborate/kleborate_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate kleborate_env

FNADIR=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES
OUTDIR=/scratch/users/k22017808/KP_Research_Project/09_Kleborate

mkdir -p "$OUTDIR"

echo "Starting Kleborate at $(date)"

kleborate \
    --assemblies "$FNADIR"/*.fna \
    --outdir "$OUTDIR" \
    --preset kpsc

echo "Finished Kleborate at $(date)"
