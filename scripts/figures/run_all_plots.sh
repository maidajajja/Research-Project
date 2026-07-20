#!/bin/bash
#SBATCH --job-name=r_plots
#SBATCH --partition=msc_appbio
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --output=/scratch/users/k22017808/KP_Research_Project/plots/plots_%j.out
#SBATCH --error=/scratch/users/k22017808/KP_Research_Project/plots/plots_%j.err

source /software/spackages_v0_21_prod/apps/linux-ubuntu22.04-zen2/gcc-13.2.0/anaconda3-2022.10-5wy43yh5crcsmws4afls5thwoskzarhe/etc/profile.d/conda.sh
conda activate myRenv

echo "Running ST distribution plot..."
Rscript ~/kp_liver_project/scripts/30_plot_ST_distribution.R

echo "Running virulence and resistance plots..."
Rscript ~/kp_liver_project/scripts/31_plot_virulence_resistance.R

echo "Running AMR heatmap..."
Rscript ~/kp_liver_project/scripts/32_plot_amr_heatmap.R

echo "Running pan-genome plots..."
Rscript ~/kp_liver_project/scripts/33_plot_pangenome.R

echo "All plots complete at $(date)"
