#!/bin/bash
#SBATCH --job-name=scoary_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/12_Scoary/scoary_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/12_Scoary/scoary_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate mlst_env

scoary \
    --genes /scratch/users/k22017808/KP_Research_Project/gene_presence_absence_scoary.csv \
    --traits /scratch/users/k22017808/KP_Research_Project/scoary_traits_ST23.csv \
    --outdir /scratch/users/k22017808/KP_Research_Project/12_Scoary \
    --correction I \
    --no-time \
    --start_col 2

echo "Finished Scoary at $(date)"
