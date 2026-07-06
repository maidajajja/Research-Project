# Comparative Genomics and AMR Profiling of Klebsiella pneumoniae from Liver Disease Patients

## Overview
MSc Applied Bioinformatics dissertation project (King's College London, 2026)
Supervisor: Dr Ellis Paintsil
Institution: Roger Williams Institute of Liver Studies, KCL

## Dataset
229 Klebsiella pneumoniae genomes from liver disease patients (final QC-passed dataset)
Sources: NCBI Assembly (152 genomes) and SRA-assembled via SPAdes (77 genomes)

## Pipeline Summary
1. Genome download and assembly (SRA genomes via SPAdes v4.2.0)
2. Quality control (QUAST, Kraken2 v2.1.2, FastANI v1.34)
3. Genome annotation (Bakta v1.12.0)
4. MLST and resistance/virulence profiling (Kleborate v3.2.4)
5. AMR profiling (AMRFinderPlus v4.2.7, ABRicate v1.2.0)
6. Phylogenetic analysis (Panaroo v1.1.2 core genome alignment, IQ-TREE v3.0.1, iTOL)
7. Pan-genome analysis (Panaroo v1.1.2)
8. Plasmid typing (MOBsuite v3.1.9, PlasmidFinder v2.1.6)
9. Pan-genome association (Scoary v1.6.16)
10. Figure generation (R scripts using ComplexHeatmap, ggplot2)

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
- R v4.x (ComplexHeatmap, ggplot2, dplyr)

## Key Results
- 229 genomes spanning 59 sequence types; ST23 dominant (n=76, 33%)
- ST11 carried the broadest acquired AMR profile including blaKPC-2
- 20/229 isolates (8.7%) met convergence criteria (virulence score >=4, resistance score >=2)
- Convergent isolates predominantly ST11 (n=17/20)
- Open pan-genome: 14,051 total genes, 3,769 core (26.8%)

## Final Figures
All final publication-quality figures are in FINAL_FIGURES/

## HPC
KCL CREATE HPC cluster (hpc.create.kcl.ac.uk)
Partition: msc_appbio
