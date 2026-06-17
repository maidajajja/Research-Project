#!/bin/bash
#SBATCH --job-name=bakta_kp126
#SBATCH --output=logs/bakta_kp126_%j.out
#SBATCH --error=logs/bakta_kp126_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --partition=interruptible_cpu

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate bakta_env
cd /scratch/users/k22017808/KP_Research_Project/06_Software/bakta-1.12.0

OUTDIR="/scratch/users/k22017808/KP_Research_Project/04_Annotations/kp_annotations/batch5"
mkdir -p $OUTDIR

genome="/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES/GCF_046599435.1_ASM4659943v1_genomic.fna"
basename=$(basename $genome .fna)

./bin/bakta --db /scratch/users/k22017808/KP_Research_Project/06_Software/db-light --output $OUTDIR/$basename --prefix $basename --skip-crispr --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-sorf --skip-ori --skip-gap --skip-plot --threads 8 $genome
BAKTA_EXIT=$?

if [ $BAKTA_EXIT -eq 0 ]; then
    echo "Bakta annotation complete for KP126"
else
    echo "Bakta annotation FAILED for KP126 (exit code $BAKTA_EXIT)"
    exit 1
fi
