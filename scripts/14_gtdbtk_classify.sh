#!/bin/bash
#SBATCH --job-name=gtdbtk_classify
#SBATCH --output=gtdbtk_classify_%j.out
#SBATCH --error=gtdbtk_classify_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=3
#SBATCH --mem=11G
#SBATCH --partition=msc_appbio

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh

conda activate gtdbtk_env

export GTDBTK_DATA_PATH="/scratch/users/k22017808/databases/gtdbtk"

mkdir -p /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/GTDBTK

gtdbtk classify_wf \
    --genome_dir /scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES \
    --extension fna \
    --out_dir /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/GTDBTK \
    --cpus 3

echo "GTDB-Tk classification complete!"
