#!/bin/bash
#SBATCH --job-name=mobsuite_remaining
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mobsuite_remaining_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mobsuite_remaining_%j.err

source ~/.bashrc
conda activate mobsuite_env

OUTDIR="/scratch/users/k22017808/KP_Research_Project/11_Plasmids/MOBsuite"
FNADIR="/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES"
DB="/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mob_db"

# Only process samples missing contig_report.txt
for dir in "$OUTDIR"/*/; do
  sample=$(basename "$dir")
  [ "$sample" == "__tmp" ] && continue
  if [ ! -f "$dir/contig_report.txt" ]; then
    # Find the FNA file
    fna=$(find "$FNADIR" -name "${sample}*.fna" -o -name "${sample}.fna" 2>/dev/null | head -1)
    if [ -z "$fna" ]; then
      echo "FNA not found for $sample - skipping"
      continue
    fi
    echo "Processing: $sample"
    mob_recon \
      --infile "$fna" \
      --outdir "$dir" \
      --database_directory "$DB" \
      --num_threads 4 \
      --force 2>&1
    echo "Done: $sample"
  fi
done

echo "MOBsuite remaining samples complete"
