#!/bin/bash
#SBATCH --job-name=amr_plasmid_new2
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --output=logs/amr_plasmid_new2_%j.out
#SBATCH --error=logs/amr_plasmid_new2_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh

FNADIR=/scratch/users/k22017808/KP_Research_Project/ALL_FNA_FILES
GENOMES="GCA_002813595.1_ASM281359v1_genomic GCF_046599435.1_ASM4659943v1_genomic"

# --- AMRFinderPlus ---
conda activate amrfinder_env
AMR_OUT=/scratch/users/k22017808/KP_Research_Project/08_AMR/AMRFinder
mkdir -p "$AMR_OUT"
echo "Starting AMRFinderPlus at $(date)"
for sample in $GENOMES; do
    amrfinder \
        --nucleotide "$FNADIR/${sample}.fna" \
        --organism Klebsiella_pneumoniae \
        --threads 4 \
        --output "$AMR_OUT/${sample}_amrfinder.tsv"
done

for sample in $GENOMES; do
    tail -n +2 "$AMR_OUT/${sample}_amrfinder.tsv" >> "$AMR_OUT/amrfinder_all_withsample.tsv"
done
echo "Finished AMRFinderPlus at $(date)"

# --- ABRicate ---
conda activate mlst_env
ABR_OUT=/scratch/users/k22017808/KP_Research_Project/08_AMR/ABRicate
echo "Starting ABRicate at $(date)"
for db in card resfinder vfdb; do
    for sample in $GENOMES; do
        abricate --db "$db" --threads 4 --minid 80 --mincov 80 "$FNADIR/${sample}.fna" \
            | tail -n +2 >> "$ABR_OUT/abricate_${db}.tsv"
    done
    abricate --summary "$ABR_OUT/abricate_${db}.tsv" > "$ABR_OUT/abricate_${db}_summary.tsv"
done
echo "Finished ABRicate at $(date)"

# --- PlasmidFinder ---
conda activate mobsuite_env
PF_OUT=/scratch/users/k22017808/KP_Research_Project/11_Plasmids/PlasmidFinder
mkdir -p "$PF_OUT"
echo "Starting PlasmidFinder at $(date)"
for sample in $GENOMES; do
    mkdir -p "$PF_OUT/${sample}"
    plasmidfinder.py \
        -i "$FNADIR/${sample}.fna" \
        -o "$PF_OUT/${sample}" \
        -p /users/k22017808/.conda/envs/mobsuite_env/share/plasmidfinder-2.1.6/database \
        -l 0.60 \
        -t 0.95
done
echo "Finished PlasmidFinder at $(date)"

echo "All AMR/plasmid analyses complete for SGH10 and KP126"
