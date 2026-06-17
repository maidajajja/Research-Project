#!/bin/bash
#SBATCH --job-name=panaroo
#SBATCH --output=panaroo_clean233_%j.out
#SBATCH --error=panaroo_clean233_%j.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=20G
#SBATCH --partition=msc_appbio

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate panaroo_env

OUTDIR=/scratch/users/k22017808/KP_Research_Project/06_Pangenome_clean233
mkdir -p ${OUTDIR}

find /scratch/users/k22017808/KP_Research_Project/04_Annotations/kp_annotations/ -name "*.gff3" > /tmp/gff_list_clean233.txt
echo "Found $(wc -l < /tmp/gff_list_clean233.txt) GFF3 files"

panaroo \
    -i $(cat /tmp/gff_list_clean233.txt | tr '\n' ' ') \
    -o ${OUTDIR} \
    --mode strict \
    -t 6 \
    -a core \
    --aligner mafft

echo "Panaroo complete! (clean 233-genome run, contaminant Salmonella genome GCA_009079735.1 excluded)"
