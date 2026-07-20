#!/bin/bash
#SBATCH --job-name=bakta_batch4
#SBATCH --output=bakta_batch4_%j.out
#SBATCH --error=bakta_batch4_%j.err
#SBATCH --time=10:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate bakta_env
cd /scratch/users/k22017808/kp_liver_project/bakta-1.12.0
OUTDIR="/scratch/users/k22017808/kp_annotations/batch4"
mkdir -p $OUTDIR

count=0
for genome in /scratch/users/k22017808/kp_liver_project/genomes_fasta/all_genomes/*.fna; do
    if [ $count -lt 110 ]; then count=$((count+1)); continue; fi
    if [ $count -ge 160 ]; then break; fi
    basename=$(basename $genome .fna)
    ./bin/bakta --db /scratch/users/k22017808/kp_liver_project/db-light --output $OUTDIR/$basename --prefix $basename --skip-crispr --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-sorf --skip-ori --skip-gap --skip-plot --threads 8 $genome
    count=$((count+1))
done
