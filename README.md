# Genomic epidemiology of Klebsiella pneumoniae  in liver disease patients: lineage-specific convergence of antimicrobial resistance and hypervirulence

## Project Overview
MSc Applied Bioinformatics
Supervisor: Prof Debbie Shawcross and Dr Ellis Paintsil
Institution: Roger Williams Institute of Liver Studies, King's College London

This project characterises the population structure, antimicrobial resistance, virulence determinants,
and MDR-hypervirulent convergence of Klebsiella pneumoniae isolated from liver disease patients,
using a publicly available genomic dataset

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
- ST258 and ST11 carried the broadest acquired AMR profiles (mean 19.0 and 18.6 genes/isolate respectively), both enriched for carbapenemase and ESBL genes
- 20/229 isolates (8.7%) met convergence criteria (virulence score >=4, resistance score >=2)
- Convergent isolates predominantly ST11 (n=17/20)
- Open pan-genome: 14,051 total genes; 3,769 core genes (26.8%)

---

## Final Figure Scripts

This table maps each file in `FINAL_FIGURES/` to the script that actually generates it, verified
against each script's `ggsave()`/output calls rather than assumed from filenames.

| File in `FINAL_FIGURES/` | Dissertation figure | Script | Description |
| --- | --- | --- | --- |
| (not script-generated) | Figure 1 | (original schematic) | K. pneumoniae species complex vs. Enterobacterales, adapted from Wyres et al. (2020) |
| (not script-generated) | Figure 2 | (original schematic) | Chromosomal vs. plasmid-encoded virulence loci, adapted from Choby et al. (2020) |
| (not script-generated) | Figure 3 | (original schematic) | MDR/hypervirulence/convergence Venn diagram |
| Fig1_FINAL | Figure 4 | `scripts/05_phylogenetics/iqtree_final229.sh` + iTOL | Core-genome phylogeny |
| Fig2_AMR_heatmap_FINAL | Figure 5 | `scripts/figures/Fig2_AMR_heatmap_v19.R` | AMR gene heatmap by ST |
| Fig3b_virulence_bubble_FINAL | Figure 6 | `scripts/figures/Fig3b_virulence_bubble_v5.R` | Virulence prevalence bubble plot |
| Fig4a_plasmid_replicons_FINAL | Figure 7 | `scripts/figures/Fig4_v5.R` | Plasmid replicon prevalence by ST |
| Fig4b_IncHI1B_convergence_FINAL | Figure 8 | `scripts/figures/Fig4_v5.R` | IncHI1B in convergent vs non-convergent ST11 |
| Fig5a_ST_over_time_FINAL | Figure 9 | `scripts/figures/Fig5a_v4.R` | ST distribution over time |
| Fig5b_ST_country_table_FINAL | Figure 10 | `scripts/figures/Fig5b_table.R` | ST distribution by country |
| Fig6a_pangenome_composition_FINAL, Fig6b_pangenome_frequency_FINAL, Fig6ab_pangenome_FINAL | Figure 11a/11b | `scripts/figures/Fig6_pangenome_v3.R` | Pan-genome composition and gene-frequency distribution |
| Fig6c_pangenome_accumulation_FINAL | Figure 12 | `scripts/figures/Fig6_pangenome_v3.R` | Pan-genome accumulation curve |
| Fig7_convergence_FINAL | Figure 13 | `scripts/figures/Fig7_convergence_v3.R` | MDR-hypervirulent convergence bubble plot |
| Fig8_convergent_plasmid_FINAL, Fig8b_convergent_integrated_FINAL | Figure 14 | `scripts/figures/Fig8_convergent_plasmid_v4.R` | Convergent isolate heatmap |
| FigS1_kraken2_QC_FINAL | Figure S1 | `scripts/figures/FigS1_Kraken2.R` | Kraken2 QC (supplementary) |
| FigS2_FastANI_heatmap_FINAL | Figure S2 | `scripts/figures/FigS2_FastANI_v3.R` | FastANI heatmap (supplementary) |
| FigS3_virulence_heatmap_FINAL | Figure S3 | `scripts/figures/FigS3_virulence_heatmap_v2.R` | Virulence loci heatmap, per-isolate (supplementary) |
| FigS4_pangenome_heatmap_FINAL | Figure S4 | `scripts/figures/FigS4_pangenome_heatmap_v4.R` | Pan-genome presence/absence heatmap (supplementary) |
| FigS5_scoary_volcano_FINAL | Figure S5 | `scripts/figures/Fig_scoary_volcano_v2.R` | Scoary volcano plot, corrected final229 rerun (supplementary) |

Note: filenames in `FINAL_FIGURES/` retain their original working names and were not renamed when
three new conceptual schematics (Figures 1-3) were added to the dissertation Introduction, which
shifted every subsequent figure's dissertation number up by 3. The "Dissertation figure" column
reflects the final numbering used in the submitted dissertation; the "File in FINAL_FIGURES/" column
reflects the actual filenames in this repository.

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
