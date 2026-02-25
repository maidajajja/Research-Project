#!/bin/bash
#SBATCH --job-name=gtdbtk_classify
#SBATCH --output=gtdbtk_classify_%j.out
#SBATCH --error=gtdbtk_classify_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=200G

source ~/.bashrc
conda activate gtdbtk_env

export GTDBTK_DATA_PATH="/scratch/users/k22017808/databases/gtdbtk"

mkdir -p /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/GTDBTK

gtdbtk classify_wf \
    --genome_dir /scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES \
    --extension fna \
    --out_dir /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/GTDBTK \
    --cpus 16

echo "GTDB-Tk classification complete!"
