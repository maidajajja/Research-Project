# Comparative Genomics and AMR Profiling of *Klebsiella pneumoniae* in Liver Disease Patients

**Author:** Maida Jajja  
**Institution:** King's College London  
**Programme:** MSc Applied Bioinformatics  
**Supervisor:** Dr Ellis Paintsil  
**Year:** 2024-2025  

---

## Project Overview

This project performs comparative genomics and antimicrobial resistance (AMR) profiling of 234 clinical *Klebsiella pneumoniae* genomes isolated from liver disease patients. The analysis spans quality control, genome assembly, annotation, phylogenetics, pan-genome analysis, AMR profiling, virulence characterisation, and plasmid detection.

---

## Repository Structure---

## Pipeline Overview

### 1. Quality Control
- **Tool:** FastQC, fastp, QUAST
- **Environment:** assembly_env
- **Input:** Raw genome assemblies
- **Output:** QC reports, filtered assemblies

### 2. Genome Assembly
- **Tool:** SPAdes v3.15.5
- **Environment:** assembly_env
- **Output:** Assembled contigs

### 3. Genome Annotation
- **Tool:** Bakta v1.12.0 (light database)
- **Environment:** bakta_env
- **Output:** GFF3, FFN, FAA, TSV annotation files

### 4. Taxonomic Classification
- **Tool:** FastANI v1.34, Kraken2 v2.1.2
- **Environment:** mlst_env
- **Output:** ANI values, taxonomic classifications

### 5. MLST
- **Tool:** MLST (PubMLST Klebsiella scheme)
- **Environment:** mlst_env
- **Output:** mlst_results.tsv — 53 distinct STs identified

### 6. Pan-genome Analysis
- **Tool:** Panaroo v1.3.4 (strict mode, MAFFT alignment)
- **Environment:** panaroo_env
- **Note:** Source code manually patched for Python 3.13 compatibility
- **Output:** Core genome alignment, gene_presence_absence.csv
- **Key results:** 23,699 total genes, 3,154 core genes (open pan-genome)

### 7. Phylogenetics
- **Tool:** IQ-TREE v2 (GTR+G model, 1000 bootstraps)
- **Environment:** iqtree_env
- **Input:** Panaroo core genome alignment (3.83 Mbp)
- **Recombination removal:** Gubbins v3
- **Output:** kp_core.treefile, visualised in iTOL

### 8. AMR Profiling
- **Tool:** AMRFinderPlus v4.2.7 (--organism Klebsiella_pneumoniae)
- **Tool:** ABRicate (CARD, ResFinder, VFDB databases)
- **Environment:** amrfinder_env
- **Key result:** blaSHV-11 in 131/234 isolates

### 9. Virulence Characterisation
- **Tool:** Kleborate v3.2.4 (kpsc preset)
- **Environment:** kleborate_env
- **Key result:** 61% isolates virulence score 4-5 (hypervirulent)
- **Key result:** 20 convergent strains (high virulence + resistance)

### 10. Plasmid Detection
- **Tool:** PlasmidFinder v2.1.6, MOB-suite
- **Environment:** mobsuite_env
- **Key result:** 214/234 carry plasmids; IncHI1B dominant in ST23

### 11. Pan-genome Association
- **Tool:** Scoary v1.6.16
- **Environment:** mlst_env

---

## Conda Environments

| Environment | Key Tools |
|-------------|-----------|
| assembly_env | SPAdes, FastQC, fastp, QUAST |
| bakta_env | Bakta v1.12.0 |
| mlst_env | MLST, ABRicate, Scoary |
| amrfinder_env | AMRFinderPlus v4.2.7 |
| kleborate_env | Kleborate v3.2.4 |
| mobsuite_env | MOB-suite, PlasmidFinder, BLAST |
| panaroo_env | Panaroo v1.3.4, MAFFT |
| iqtree_env | IQ-TREE v2, Gubbins |
| myRenv | R v4.3.1, ggplot2, ComplexHeatmap, pheatmap |

---

## HPC Configuration

- **Cluster:** KCL CREATE HPC
- **Scheduler:** SLURM
- **Partition:** msc_appbio
- **Resources:** 6 CPUs, 22GB RAM, max 2 concurrent jobs
- **Data location:** /scratch/users/k22017808/KP_Research_Project/

---

## Key Findings

1. ST23 is the dominant hypervirulent lineage (76/234 isolates)
2. Open pan-genome confirmed (23,699 genes total, 3,154 core)
3. blaSHV-11 most prevalent AMR gene (131/234 isolates)
4. 20 convergent strains carrying both high virulence and resistance
5. IncHI1B virulence plasmid predominantly in ST23 isolates

---

## Citation

If using this pipeline, please cite the individual tools listed above.
