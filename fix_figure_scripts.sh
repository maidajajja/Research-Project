#!/usr/bin/env bash
# Fix the mismatch between FINAL_FIGURES/ and scripts/figures/.
# Run this from the ROOT of your local clone (same place as reorganize_repo.sh ran).
#
# What it does:
#   - Moves Fig5b_table.R from archive into scripts/figures/ (it's the ACTUAL script
#     that generates FINAL_FIGURES/Fig5b_ST_country_table_FINAL - it was sitting in
#     archive/, so that figure was not reproducible from the "final" scripts folder)
#   - Renames Fig3a_virulence_heatmap_v2.R -> FigS3_virulence_heatmap_v2.R (it actually
#     generates FigS3 output, not Fig3a - confirmed from its own ggsave() calls)
#   - Archives 3 scripts whose output does NOT appear anywhere in FINAL_FIGURES/:
#     FigS1_FastANI_QC.R, FigS2_SNP_heatmap.R, publication_plots_v2.R
#   - Replaces the README's Final Figure Scripts table with one verified against each
#     script's actual ggsave()/output filenames, not assumed from script names
#
# Does NOT commit or push - review with git status / git diff --cached --stat first.

set -euo pipefail

if [ ! -f "README.md" ] || [ ! -d "scripts/figures" ]; then
  echo "Run this from the repo root, after reorganize_repo.sh has already been run and pushed." >&2
  exit 1
fi

echo "== Moving Fig5b's real script into scripts/figures/ =="
git mv scripts/archive/fig4_fig5_old_versions/Fig5b_table.R scripts/figures/Fig5b_table.R

echo "== Renaming mislabelled virulence heatmap script =="
git mv scripts/figures/Fig3a_virulence_heatmap_v2.R scripts/figures/FigS3_virulence_heatmap_v2.R

echo "== Archiving 3 scripts whose output is not in FINAL_FIGURES/ =="
git mv scripts/figures/FigS1_FastANI_QC.R scripts/archive/old_figure_versions/FigS1_FastANI_QC.R
git mv scripts/figures/FigS2_SNP_heatmap.R scripts/archive/old_figure_versions/FigS2_SNP_heatmap.R
git mv scripts/figures/publication_plots_v2.R scripts/archive/old_figure_versions/publication_plots_v2.R

echo "== Updating README =="
rm -f README.md
cat > README.md << 'README_EOF'
# Comparative Genomics and AMR Profiling of Klebsiella pneumoniae in Liver Disease Patients

## Project Overview
MSc Applied Bioinformatics | King's College London | 2026
Supervisor: Dr Ellis Paintsil
Institution: Roger Williams Institute of Liver Studies, KCL

This project characterises the population structure, antimicrobial resistance, virulence determinants,
and MDR-hypervirulent convergence of Klebsiella pneumoniae isolated from liver disease patients,
using a publicly available genomic dataset of 229 isolates.

---

## Dataset
- 229 genomes (final QC-passed dataset - see `docs/qc_summary_final229.md` for full exclusion/addition detail)
- Sources: NCBI Assembly and SRA-assembled via SPAdes
- Host: Liver disease patients (liver abscess, liver transplant, other liver disease)
- Geography: Asia (n=175), Europe (n=33), North America (n=20), Unknown (n=1)
- Years: 2007-2023
- Metadata: `data/genomes_master.csv` - definitive metadata file

---

## Repository Structure

```
README.md
environment/              Conda environment YAMLs for each tool
data/
  genomes_master.csv       Definitive metadata file
  lists/                   Accession lists used during genome acquisition
docs/
  qc_summary_final229.md   Final QC exclusion/addition summary (matches dissertation Methods/Results)
  archive/                 Superseded documentation, retained for audit history
scripts/
  01_data_acquisition/     Genome download (NCBI assemblies, GenBank, SRA) and SPAdes assembly
  02_annotation/           Bakta annotation batches
  03_quality_control/      QUAST, Kraken2, FastANI species/contamination/assembly checks
  04_typing_and_amr/       Kleborate MLST/virulence/resistance scoring, AMRFinderPlus, ABRicate
  05_phylogenetics/        Panaroo core-genome alignment, IQ-TREE, Gubbins recombination filtering
  06_plasmid_analysis/     MOB-suite and PlasmidFinder replicon typing
  07_pangenome_association/  Scoary gene-association analysis
  figures/                 Final figure-generation scripts (Fig1-Fig8, FigS1-FigS2)
  archive/                 Superseded script versions, retained for reproducibility
iTOL/
  final/                   Definitive iTOL annotation files for the phylogenetic tree
  archive/                 Superseded iTOL versions
FINAL_FIGURES/             Publication-quality figures (PNG and PDF) - see note below
```

Numeric prefixes on individual scripts were removed in favour of staged folders (01-07 reflect
pipeline order); within each folder, script names describe what they do rather than a step number,
since several tools (e.g. Kraken2, Bakta) were re-run at more than one pipeline stage.

---

## Bioinformatic Pipeline

1. **Data acquisition** (`scripts/01_data_acquisition/`) - NCBI Assembly and GenBank download, SRA raw-read retrieval and SPAdes v4.2.0 assembly
2. **Annotation** (`scripts/02_annotation/`) - Bakta v1.12.0
3. **Quality control** (`scripts/03_quality_control/`) - QUAST v5.2.0, Kraken2 v2.1.2, FastANI v1.34 (species ID, contamination, assembly fragmentation, CDS-density checks; see `docs/qc_summary_final229.md`)
4. **Typing, AMR and virulence** (`scripts/04_typing_and_amr/`) - Kleborate v3.2.4 (MLST/virulence/resistance scores), AMRFinderPlus v4.2.7 (db 2026-03-24.1), ABRicate v1.2.0
5. **Phylogenetics** (`scripts/05_phylogenetics/`) - Panaroo v1.1.2 core-genome alignment, IQ-TREE v3.0.1 (GTR+G, 1000 bootstrap), Gubbins v2.4.1 recombination filtering
6. **Plasmid analysis** (`scripts/06_plasmid_analysis/`) - MOB-suite v3.1.9, PlasmidFinder v2.1.6
7. **Pan-genome association** (`scripts/07_pangenome_association/`) - Scoary v1.6.16, ST23 vs non-ST23
8. **Figures** (`scripts/figures/`) - R (ComplexHeatmap, ggplot2)

---

## Software Versions
- SPAdes v4.2.0
- Bakta v1.12.0
- Kraken2 v2.1.2
- FastANI v1.34
- Kleborate v3.2.4
- AMRFinderPlus v4.2.7 (database 2026-03-24.1)
- ABRicate v1.2.0
- Panaroo v1.1.2
- IQ-TREE v3.0.1
- MOBsuite v3.1.9
- PlasmidFinder v2.1.6
- Scoary v1.6.16
- Gubbins v2.4.1
- R v4.x (ComplexHeatmap, ggplot2, dplyr, circlize)

---

## Key Results
- 59 distinct sequence types identified; ST23 dominant (n=76, 33%)
- ST11 carried the broadest acquired AMR profile, enriched for blaKPC-2
- 20/229 isolates (8.7%) met convergence criteria (virulence score >=4, resistance score >=2)
- Convergent isolates predominantly ST11 (n=17/20)
- Open pan-genome: 14,051 total genes; 3,769 core genes (26.8%)

---

## Final Figure Scripts

This table maps each file in `FINAL_FIGURES/` to the script that actually generates it, verified
against each script's `ggsave()`/output calls rather than assumed from filenames.

| Final figure(s) | Script | Description |
|-----------------|--------|--------------|
| Fig1 | `scripts/05_phylogenetics/iqtree_final229.sh` + iTOL | Core-genome phylogeny |
| Fig2_AMR_heatmap | `scripts/figures/Fig2_AMR_heatmap_v19.R` | AMR gene heatmap by ST |
| Fig3b_virulence_bubble | `scripts/figures/Fig3b_virulence_bubble_v5.R` | Virulence prevalence bubble plot |
| Fig4a_plasmid_replicons, Fig4b_IncHI1B_convergence | `scripts/figures/Fig4_v5.R` | Plasmid replicon prevalence by ST; IncHI1B in convergent vs non-convergent ST11 |
| Fig5a_ST_over_time | `scripts/figures/Fig5a_v4.R` | ST distribution over time |
| Fig5b_ST_country_table | `scripts/figures/Fig5b_table.R` | ST distribution by country |
| Fig6a_pangenome_composition, Fig6b_pangenome_frequency, Fig6d_pangenome_accumulation, Fig6ab_pangenome (combined panel) | `scripts/figures/Fig6_pangenome_v3.R` | Pan-genome composition, gene-frequency distribution, and accumulation curve |
| Fig6c_pangenome_heatmap | `scripts/figures/Fig6c_pangenome_heatmap_v4.R` | Pan-genome presence/absence heatmap |
| Fig7_convergence | `scripts/figures/Fig7_convergence_v3.R` | MDR-hypervirulent convergence bubble plot |
| Fig8_convergent_plasmid | `scripts/figures/Fig8_convergent_plasmid_v4.R` | Convergent isolate heatmap |
| Fig8b_convergent_integrated | `scripts/figures/Fig8b_convergent_integrated_v1.R` | Integrated convergence figure |
| FigS1_kraken2_QC | `scripts/figures/FigS1_Kraken2.R` | Kraken2 QC supplementary figure |
| FigS2_FastANI_heatmap | `scripts/figures/FigS2_FastANI_v3.R` | FastANI supplementary heatmap |
| FigS3_virulence_heatmap | `scripts/figures/FigS3_virulence_heatmap_v2.R` | Virulence loci heatmap (supplementary). Note: the script file was previously misnamed `Fig3a_virulence_heatmap_v2.R` despite generating FigS3 output; renamed for consistency. |

`scripts/figures/run_all_plots.sh` regenerates the full figure set in one pass.

`scripts/archive/old_figure_versions/FigS1_FastANI_QC.R`, `FigS2_SNP_heatmap.R`, and
`publication_plots_v2.R` were moved out of `scripts/figures/` because none of their outputs appear
in `FINAL_FIGURES/` - they were exploratory/superseded and are retained in the archive for history
only, not as part of the current figure set.

---

## FINAL_FIGURES/

This folder holds the publication-quality PNG/PDF exports of every figure above, generated on the
KCL CREATE HPC cluster (`/scratch/users/k22017808/KP_Research_Project/FINAL_FIGURES/`).

---

## HPC
Cluster: KCL CREATE HPC (hpc.create.kcl.ac.uk)
Partition: msc_appbio
Conda environments: myRenv (R packages), kleborate_env (Kleborate) - see `environment/` for full YAMLs
README_EOF

git add -A

echo ""
echo "Done. Review with:  git status   and   git diff --cached --stat"
echo "Then commit, e.g.:"
echo "  git commit -m \"Fix scripts/figures to match FINAL_FIGURES: recover Fig5b script, rename FigS3 script, archive 3 unused figure scripts, correct README table\""
echo "  git push"
