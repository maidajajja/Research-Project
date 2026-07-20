#!/bin/bash
#SBATCH --job-name=iqtree
#SBATCH --output=iqtree_final229_%j.out
#SBATCH --error=iqtree_final229_%j.err
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=20G
#SBATCH --partition=msc_appbio

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate iqtree_env

ALIGNMENT=/scratch/users/k22017808/KP_Research_Project/06_Pangenome_final229/core_gene_alignment.aln
OUTDIR=/scratch/users/k22017808/KP_Research_Project/07_Phylogeny_final229
mkdir -p ${OUTDIR}

iqtree -s ${ALIGNMENT} \
       -m GTR+G \
       -bb 1000 \
       -nt AUTO \
       -ntmax 4 \
       --prefix ${OUTDIR}/kp_core_final229

echo "IQ-TREE complete (final 229-genome core alignment, IQ-TREE v3.0.1, GTR+G, 1000 bootstrap replicates)"
