#!/bin/bash
#SBATCH --job-name=kraken2_db
#SBATCH --output=kraken2_db_%j.out
#SBATCH --error=kraken2_db_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=msc_appbio

module load kraken2/2.1.2-gcc-13.2.0

mkdir -p /scratch/users/k22017808/databases/kraken2

cd /scratch/users/k22017808/databases/kraken2

wget -q https://genome-idx.s3.amazonaws.com/kraken/k2_standard_08gb_20240904.tar.gz \
    -O k2_standard_08gb.tar.gz

tar -xzf k2_standard_08gb.tar.gz

rm k2_standard_08gb.tar.gz

ls -lh /scratch/users/k22017808/databases/kraken2/

echo "Kraken2 database ready!"
