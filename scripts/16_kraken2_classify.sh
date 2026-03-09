#!/bin/bash
#SBATCH --job-name=kraken2_classify
#SBATCH --output=kraken2_classify_%j.out
#SBATCH --error=kraken2_classify_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --partition=msc_appbio

module load kraken2/2.1.2-gcc-13.2.0

DB=/scratch/users/k22017808/databases/kraken2
OUTDIR=/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/Kraken2
GENOMES=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES

mkdir -p ${OUTDIR}

for genome in ${GENOMES}/*.fna; do
    sample=$(basename ${genome} .fna)
    kraken2 --db ${DB} \
            --threads 6 \
            --output ${OUTDIR}/${sample}.kraken \
            --report ${OUTDIR}/${sample}.report \
            ${genome}
done

echo "Kraken2 classification complete!"
