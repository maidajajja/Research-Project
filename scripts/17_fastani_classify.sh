#!/bin/bash
#SBATCH --job-name=fastani_ref
#SBATCH --output=fastani_ref_%j.out
#SBATCH --error=fastani_ref_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=msc_appbio

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate gtdbtk_env

mkdir -p /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/FastANI

ls /scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES/*.fna > /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/FastANI/genome_list.txt

fastANI --ql /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/FastANI/genome_list.txt \
        --rl /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/FastANI/genome_list.txt \
        -o /scratch/users/k22017808/KP_Research_Project/05_Taxonomy/FastANI/fastani_results.txt \
        --minFraction 0.1 \
        -t 2

echo "FastANI complete!"
