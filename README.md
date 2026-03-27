# Comparative Genomics of *Klebsiella pneumoniae* in Liver Disease Patients

> MSc Applied Bioinformatics Dissertation | King's College London | 2025–2026
> Supervisor: Dr Ellis Paintsil

---

## Project Summary

This repository contains the bioinformatics pipeline and analysis scripts for a comparative genomic study of **234 *Klebsiella pneumoniae* genomes** derived from liver disease patients. Using whole genome sequencing data sourced from NCBI GenBank and the Sequence Read Archive, this project characterises the population structure, antimicrobial resistance gene repertoire, virulence factor profiles, pan-genome diversity, and phylogenetic relationships of *K. pneumoniae* circulating in this clinically vulnerable patient group.

---

## Repository Structure
```
Research-Project/
├── scripts/           # SLURM job submission scripts
├── lists/             # Genome accession lists
├── logs/              # Pipeline log files
└── quast_results/     # Assembly quality assessment metrics
```

---

## Dataset

| Source | Count | Format |
|--------|-------|--------|
| NCBI Assembly (GCA/GCF) | 152 | Pre-assembled FASTA |
| NCBI SRA | 82 | Raw Illumina reads → assembled |
| **Total** | **234** | *K. pneumoniae* genomes |

All genomes were confirmed as *Klebsiella pneumoniae* by FastANI (≥98.8% ANI) and Kraken2.

---

## Pipeline

### 1. Data Acquisition
| Step | Tool | Script |
|------|------|--------|
| Download NCBI assemblies | NCBI Datasets | `01_download_assemblies.sh` |
| Download GenBank isolates | Custom Python | `02_download_genbank_isolates.py` |
| Download additional assemblies | Custom Python | `05_download_genbank_genomes.py` |
| SRA assembly pipeline | SPAdes v3.15.5 | `06_sra_assembly_pipeline.sh` |

### 2. Assembly Quality Control
| Step | Tool | Script |
|------|------|--------|
| Assembly QC | QUAST | `07_quast_all_genomes.sh` |
| Consolidate genomes | Custom | `03_make_all_genomes.sh` |

### 3. Genome Annotation
| Step | Tool | Script |
|------|------|--------|
| Annotation (batched) | Bakta v1.12.0 (light db) | `08_annotate_batch1.sh` → `12_annotate_fasta.sh` |

### 4. Taxonomy & Sequence Typing
| Step | Tool | Script |
|------|------|--------|
| k-mer classification | Kraken2 v2.1.2 | `16_kraken2_classify.sh` |
| ANI-based classification | FastANI v1.34 | `17_fastani_classify.sh` |
| Multi-locus sequence typing | MLST (PubMLST Klebsiella) | `22_mlst.sh` |

### 5. Pan-genome & Phylogeny
| Step | Tool | Script |
|------|------|--------|
| Pan-genome analysis | Panaroo v1.3.4 | `18_panaroo_pangenome.sh` |
| Core genome alignment | Panaroo + MAFFT | `19_panaroo_alignment.sh` |
| Phylogenetic tree | IQ-TREE v2 (GTR+G, 1000 bootstraps) | — |

### 6. Resistome & Virulome
| Step | Tool | Script |
|------|------|--------|
| AMR gene detection | AMRFinderPlus v4.2.7 | `23_amrfinder.sh` |
| AMR & virulence screening | ABRicate (CARD, ResFinder, VFDB) | `24_abricate.sh` |
| Resistance + virulence scoring | Kleborate v3.2.4 | `25_kleborate.sh` |

### 7. Downstream Analyses (In Progress)
| Step | Tool |
|------|------|
| Recombination detection | Gubbins |
| Plasmid analysis | MOB-suite, PlasmidFinder |
| Population structure | PopPUNK |
| Pan-genome association | Scoary |
| Functional annotation | eggNOG-mapper |

---

## Key Results

| Metric | Value |
|--------|-------|
| Total genomes | 234 |
| Core genome size | 3,154 genes |
| Total pan-genome | 23,699 genes |
| Dominant lineage | ST23 (n=23) |
| Most prevalent AMR gene | fosA (205/234 genomes) |
| ESBL gene detected | blaSHV-11 (131/234 genomes) |
| Hypervirulence marker (iutA) | 396 hits across collection |

---

## HPC Environment

- **Cluster:** KCL CREATE HPC
- **Scheduler:** SLURM
- **Partition:** `msc_appbio` (6 CPUs, 22GB RAM, 2 jobs max)
- **Working directory:** `/scratch/users/k22017808/KP_Research_Project/`

### Conda Environments

| Environment | Tools |
|-------------|-------|
| `assembly_env` | SPAdes |
| `bakta_env` | Bakta |
| `mlst_env` | MLST, ABRicate |
| `amrfinder_env` | AMRFinderPlus |
| `kleborate_env` | Kleborate |
| `panaroo_env` | Panaroo, MAFFT |
| `iqtree_env` | IQ-TREE |
| `gubbins_env` | Gubbins (installing) |

---

## Notable Technical Details

- Panaroo source code manually patched for Python 3.13 compatibility (Bio.Alphabet removal, rU mode fix, MAFFT subprocess handling)
- Bakta source code patched to handle AMRFinder dependency crashes in HPC environment
- All scripts version controlled and submitted via SLURM

---

## Author

**Maida Jajja**
MSc Applied Bioinformatics, King's College London
GitHub: [@maidajajja](https://github.com/maidajajja)# Research-Project
