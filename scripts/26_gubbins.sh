#!/bin/bash
#SBATCH --job-name=gubbins_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --time=12:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/10_Gubbins/gubbins_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/10_Gubbins/gubbins_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate iqtree_env

OUTDIR=/scratch/users/k22017808/KP_Research_Project/10_Gubbins
ALN=/scratch/users/k22017808/KP_Research_Project/06_Pangenome_alignment/core_gene_alignment.aln

mkdir -p "$OUTDIR"
cd "$OUTDIR"

echo "Starting Gubbins at $(date)"

run_gubbins.py --prefix kp_gubbins --threads 6 "$ALN"

echo "Finished Gubbins at $(date)"
