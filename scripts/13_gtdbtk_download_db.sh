#!/bin/bash
#SBATCH --job-name=gtdbtk_db
#SBATCH --output=gtdbtk_db_%j.out
#SBATCH --error=gtdbtk_db_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

source ~/.bashrc
conda activate gtdbtk_env

mkdir -p /scratch/users/k22017808/databases/gtdbtk

wget https://data.gtdb.aau.ecogenomic.org/releases/release226/226.0/auxillary_files/gtdbtk_package/full_package/gtdbtk_r226_data.tar.gz \
    -P /scratch/users/k22017808/databases/

tar -xvzf /scratch/users/k22017808/databases/gtdbtk_r226_data.tar.gz \
    -C /scratch/users/k22017808/databases/gtdbtk --strip 1

rm /scratch/users/k22017808/databases/gtdbtk_r226_data.tar.gz

conda env config vars set GTDBTK_DATA_PATH="/scratch/users/k22017808/databases/gtdbtk"

echo "Database download and extraction complete!"
