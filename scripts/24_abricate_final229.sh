#!/bin/bash
#SBATCH --job-name=abricate_kp
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=04:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/08_AMR/ABRicate/abricate_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/08_AMR/ABRicate/abricate_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate mlst_env

FNADIR=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES
OUTDIR=/scratch/users/k22017808/KP_Research_Project/08_AMR/ABRicate_final229

mkdir -p "$OUTDIR"

echo "Starting ABRicate at $(date)"

for db in card resfinder vfdb; do
    echo "Running database: $db"
    > "$OUTDIR/abricate_${db}.tsv"
    for fna in "$FNADIR"/*.fna; do
        abricate --db "$db" --threads 4 --minid 80 --mincov 80 "$fna" | tail -n +2 >> "$OUTDIR/abricate_${db}.tsv"
    done
    abricate --summary "$OUTDIR/abricate_${db}.tsv" > "$OUTDIR/abricate_${db}_summary.tsv"
    echo "Finished $db at $(date)"
done

echo "Finished ABRicate at $(date)"
