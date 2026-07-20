#!/bin/bash
#SBATCH --job-name=bakta_batch1
#SBATCH --output=bakta_batch1_%j.out
#SBATCH --error=bakta_batch1_%j.err
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G

source ~/.bashrc
conda activate bakta_env
cd /scratch/users/k22017808/kp_liver_project/bakta-1.12.0
OUTDIR="/scratch/users/k22017808/kp_annotations/batch1"
mkdir -p $OUTDIR

count=0
for genome in /scratch/users/k22017808/kp_liver_project/genomes_fasta/all_genomes/*.fna; do
    if [ $count -ge 10 ]; then break; fi
    basename=$(basename $genome .fna)
    ./bin/bakta --db /scratch/users/k22017808/kp_liver_project/db-light --output $OUTDIR/$basename --prefix $basename --skip-crispr --skip-trna --skip-tmrna --skip-rrna --skip-ncrna --skip-ncrna-region --skip-sorf --skip-ori --skip-gap --skip-plot --threads 8 $genome
    count=$((count+1))
done
