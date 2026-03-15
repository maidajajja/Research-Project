#!/bin/bash
#SBATCH --job-name=panaroo_aln
#SBATCH --output=panaroo_aln_%j.out
#SBATCH --error=panaroo_aln_%j.err
#SBATCH --time=2-00:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --partition=msc_appbio

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate panaroo_env

OUTPUT=/scratch/users/k22017808/KP_Research_Project/06_Pangenome_alignment

mkdir -p $OUTPUT

panaroo \
    -i /scratch/users/k22017808/KP_Research_Project/ALL_GFF3_FILES/*.gff3 \
    -o $OUTPUT \
    --mode strict \
    --alignment core \
    --aligner mafft \
    --core_threshold 0.98 \
    -t 6

echo "Panaroo alignment complete!"
