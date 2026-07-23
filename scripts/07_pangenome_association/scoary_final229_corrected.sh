#!/bin/bash
#SBATCH --job-name=scoary_kp229
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/12_Scoary_final229/scoary_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/12_Scoary_final229/scoary_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate mlst_env

echo "Starting Scoary (final229) at $(date)"

scoary --genes /scratch/users/k22017808/KP_Research_Project/06_Pangenome_final229/gene_presence_absence.Rtab \
       --traits /scratch/users/k22017808/KP_Research_Project/scoary_traits_ST23_final229.csv \
       --outdir /scratch/users/k22017808/KP_Research_Project/12_Scoary_final229 \
       --correction I --no-time --start_col 2

echo "Finished Scoary (final229) at $(date)"
