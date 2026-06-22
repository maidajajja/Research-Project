#!/bin/bash
#SBATCH --job-name=gtdbtk_db2
#SBATCH --output=gtdbtk_db2_%j.out
#SBATCH --error=gtdbtk_db2_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --partition=msc_appbio

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate gtdbtk_env

mkdir -p /scratch/users/k22017808/databases/gtdbtk

wget https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz \
    -P /scratch/users/k22017808/databases/

tar -xvzf /scratch/users/k22017808/databases/gtdbtk_r226_data.tar.gz \
    -C /scratch/users/k22017808/databases/gtdbtk --strip 1

rm /scratch/users/k22017808/databases/gtdbtk_r226_data.tar.gz

echo "Database download and extraction complete!"
