#!/bin/bash
#SBATCH --job-name=mobsuite_fix
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=16:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mobsuite_fix_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mobsuite_fix_%j.err

source /users/k22017808/.bashrc
source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate mobsuite_env

OUTDIR="/scratch/users/k22017808/KP_Research_Project/11_Plasmids/MOBsuite"
FNADIR="/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES"
DB="/scratch/users/k22017808/KP_Research_Project/11_Plasmids/mob_db"

echo "mob_recon path: $(which mob_recon)"
echo "Starting MOBsuite fix run"

for dir in "$OUTDIR"/*/; do
  sample=$(basename "$dir")
  [ "$sample" == "__tmp" ] && continue
  if [ ! -f "$dir/contig_report.txt" ]; then
    fna="$FNADIR/${sample}.fna"
    if [ ! -f "$fna" ]; then
      echo "FNA not found for $sample"
      continue
    fi
    echo "Processing: $sample"
    mob_recon \
      --infile "$fna" \
      --outdir "$dir" \
      --database_directory "$DB" \
      --num_threads 4 \
      --force
    echo "Done: $sample"
  fi
done

echo "Completed. Total with contig_report.txt:"
find "$OUTDIR" -name "contig_report.txt" | wc -l
