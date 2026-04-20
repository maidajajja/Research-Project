# Comparative Genomics and AMR Profiling of Klebsiella pneumoniae from Liver Disease Patients

## Overview
MSc Applied Bioinformatics dissertation project (King's College London, 2026)
Supervisor: Dr Ellis Paintsil

## Dataset
234 Klebsiella pneumoniae genomes from liver disease patients
Sources: NCBI/SRA public databases

## Pipeline Summary
1. Genome download and assembly (SRA genomes via SPAdes)
2. Genome annotation (Bakta)
3. Quality control (QUAST, Kraken2, FastANI)
4. MLST typing (mlst)
5. Phylogenetic analysis (Panaroo core genome alignment, IQ-TREE, iTOL)
6. AMR profiling (AMRFinderPlus, Kleborate)
7. Virulence profiling (Kleborate)
8. Pan-genome analysis (Panaroo)
9. Plasmid typing (PlasmidFinder)
10. Figure generation (R scripts)

## Software Versions
- Bakta v1.x
- Panaroo v1.3.4
- AMRFinderPlus v3.x
- Kleborate v2.x
- PlasmidFinder v2.1.6
- R v4.x

## HPC
KCL CREATE HPC cluster
