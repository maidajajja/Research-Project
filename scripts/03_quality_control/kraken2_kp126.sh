#!/bin/bash
#SBATCH --job-name=kraken2_kp126
#SBATCH --output=logs/kraken2_kp126_%j.out
#SBATCH --error=logs/kraken2_kp126_%j.err
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

sample=GCF_046599435.1_ASM4659943v1_genomic
kraken2 --db ${DB} \
        --threads 6 \
        --output ${OUTDIR}/${sample}.kraken \
        --report ${OUTDIR}/${sample}.report \
        ${GENOMES}/${sample}.fna

echo "Kraken2 classification complete for KP126"
