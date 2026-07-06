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
- 229 genomes (final QC-passed dataset)
- Sources: NCBI Assembly (152 genomes) and SRA-assembled via SPAdes (77 genomes)
- Host: Liver disease patients (liver abscess, liver transplant, other liver disease)
- Geography: Asia (n=175), Europe (n=33), North America (n=20), Unknown (n=1)
- Years: 2007-2023
- Metadata: genomes_master.csv — definitive metadata file

---

## Bioinformatic Pipeline

Step 1 - Genome download: NCBI datasets — 01_download_assemblies.sh
Step 2 - SRA assembly: SPAdes v4.2.0 — 06_sra_assembly_pipeline.sh
Step 3 - Quality control: QUAST, Kraken2 v2.1.2, FastANI v1.34 — 07_quast_all_genomes.sh, 16_kraken2_classify.sh, 17_fastani_classify.sh
Step 4 - Genome annotation: Bakta v1.12.0 — 01-12_annotate_*.sh
Step 5 - MLST + AMR + virulence: Kleborate v3.2.4 — 25_kleborate.sh
Step 6 - AMR profiling: AMRFinderPlus v4.2.7 (db 2026-03-24.1) — 23_amrfinder_final229.sh
Step 7 - Phylogenetic analysis: Panaroo v1.1.2 + IQ-TREE v3.0.1 — 26_panaroo_final229.sh, 27_iqtree_final229.sh
Step 8 - Pan-genome analysis: Panaroo v1.1.2 — 26_panaroo_final229.sh
Step 9 - Plasmid typing: MOBsuite v3.1.9, PlasmidFinder v2.1.6 — 27_mobsuite.sh, 29_plasmidfinder.sh
Step 10 - Figure generation: R (ComplexHeatmap, ggplot2) — Fig*.R scripts

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

## Repository Structure
scripts/                   Final pipeline and figure scripts
scripts/archive/           Intermediate versions retained for reproducibility
FINAL_FIGURES/             Publication-quality figures (PNG and PDF)
genomes_master.csv         Definitive metadata file

---

## Final Figure Scripts
Fig2  - Fig2_AMR_heatmap_v19.R            AMR gene heatmap by ST
Fig3a - Fig3a_virulence_heatmap_v2.R      Virulence loci heatmap
Fig3b - Fig3b_virulence_bubble_v5.R       Virulence prevalence bubble plot
Fig4  - Fig4_v5.R                         Plasmid replicon prevalence by ST
Fig5a - Fig5a_v4.R                        ST distribution over time
Fig6  - Fig6_pangenome_v3.R               Pan-genome composition and frequency
Fig6c - Fig6c_pangenome_heatmap_v4.R      Pan-genome presence/absence heatmap
Fig7  - Fig7_convergence_v3.R             MDR-hypervirulent convergence bubble plot
Fig8  - Fig8_convergent_plasmid_v4.R      Convergent isolate heatmap
Fig8b - Fig8b_convergent_integrated_v1.R  Integrated convergence figure

---

## HPC
Cluster: KCL CREATE HPC (hpc.create.kcl.ac.uk)
Partition: msc_appbio
Conda environments: myRenv (R packages), kleborate_env (Kleborate)
