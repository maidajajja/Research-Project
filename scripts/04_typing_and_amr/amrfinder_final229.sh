#!/bin/bash
#SBATCH --job-name=amrfinder_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=16G
#SBATCH --time=06:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder/amrfinder_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder/amrfinder_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate amrfinder_env

INDIR=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES
OUTDIR=/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder_final229

mkdir -p "$OUTDIR"

echo "Starting AMRFinderPlus at $(date)"

for fna in "$INDIR"/*.fna; do
    sample=$(basename "$fna" .fna)
    echo "Processing $sample"
    amrfinder \
        --nucleotide "$fna" \
        --organism Klebsiella_pneumoniae \
        --threads 6 \
        --output "$OUTDIR/${sample}_amrfinder.tsv"
done

# Combine all results
head -1 "$OUTDIR"/$(ls "$OUTDIR" | head -1) > "$OUTDIR"/amrfinder_all_final229.tsv
for f in "$OUTDIR"/*_amrfinder.tsv; do
    tail -n +2 "$f" >> "$OUTDIR"/amrfinder_all_final229.tsv
done

echo "Finished AMRFinderPlus at $(date)"
