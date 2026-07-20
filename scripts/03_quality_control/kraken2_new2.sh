#!/bin/bash
#SBATCH --job-name=kraken2_new2
#SBATCH --output=logs/kraken2_new2_%j.out
#SBATCH --error=logs/kraken2_new2_%j.err
#SBATCH --time=2:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --partition=msc_appbio

module load kraken2/2.1.2-gcc-13.2.0

DB=/scratch/users/k22017808/databases/kraken2
OUTDIR=/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/Kraken2
GENOMES=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES

mkdir -p ${OUTDIR}

for sample in GCA_002268655.1_ASM226865v1_genomic GCA_002813595.1_ASM281359v1_genomic; do
    kraken2 --db ${DB} \
            --threads 6 \
            --output ${OUTDIR}/${sample}.kraken \
            --report ${OUTDIR}/${sample}.report \
            ${GENOMES}/${sample}.fna
done

echo "Kraken2 classification complete for new genomes (H5, SGH10)!"
