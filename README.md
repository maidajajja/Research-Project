# Comparative Genomics of *Klebsiella pneumoniae* in Liver Disease Patients

> MSc Applied Bioinformatics Dissertation | King's College London | 2026
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

All genomes confirmed as *Klebsiella pneumoniae* by FastANI (≥98.8% ANI) and Kraken2.

---

## Pipeline

### 1. Data Acquisition
| Step | Tool | Script |
|------|------|--------|
| Download NCBI assemblies | NCBI Datasets | `01_download_assemblies.sh` |
| Download GenBank isolates | Custom Python | `02_download_genbank_isolates.py` |
| Download additional assemblies | Custom Python | `05_download_genbank_genomes.py` |
| Download new assemblies | Custom | `04_download_new_assemblies.sh` |
| SRA assembly pipeline | SPAdes v3.15.5 | `06_sra_assembly_pipeline.sh` |

### 2. Assembly Quality Control
| Step | Tool | Script |
|------|------|--------|
| Assembly QC | QUAST | `07_quast_all_genomes.sh` |
| Consolidate genomes | Custom | `03_make_all_genomes.sh` |

### 3. Genome Annotation
| Step | Tool | Script |
|------|------|--------|
| Annotation batch 1 | Bakta v1.12.0 (light db) | `08_annotate_batch1_WORKED.sh` |
| Annotation batch 2 | Bakta v1.12.0 (light db) | `09_annotate_batch2_WORKED.sh` |
| Annotation batch 3 | Bakta v1.12.0 (light db) | `10_annotate_batch3_WORKED.sh` |
| Annotation batch 4 | Bakta v1.12.0 (light db) | `11_annotate_batch4_WORKED.sh` |
| Annotation FASTA genomes | Bakta v1.12.0 (light db) | `12_annotate_fasta_WORKED.sh` |

### 4. Taxonomy & Sequence Typing
| Step | Tool | Script |
|------|------|--------|
| Kraken2 database download | Kraken2 v2.1.2 | `15_kraken2_db_download.sh` |
| k-mer classification | Kraken2 v2.1.2 | `16_kraken2_classify.sh` |
| ANI-based classification | FastANI v1.34 | `17_fastani_classify.sh` |
| Multi-locus sequence typing | MLST (PubMLST Klebsiella) | `22_mlst.sh` |

### 5. Pan-genome & Phylogeny
| Step | Tool | Script |
|------|------|--------|
| Pan-genome analysis | Panaroo v1.3.4 | `18_panaroo_pangenome.sh` |
| Core genome alignment | Panaroo + MAFFT | `19_panaroo_alignment.sh` |
| Phylogenetic tree | IQ-TREE v2 (GTR+G, 1000 bootstraps) | — |
| Recombination analysis | Gubbins | `26_gubbins.sh` |

### 6. Resistome & Virulome
| Step | Tool | Script |
|------|------|--------|
| AMR gene detection | AMRFinderPlus v4.2.7 | `23_amrfinder.sh` |
| AMR screening (CARD + ResFinder) | ABRicate | `24_abricate.sh` |
| Virulence + resistance scoring | Kleborate v3.2.4 (kpsc) | `25_kleborate.sh` |

### 7. Plasmid Analysis
| Step | Tool | Script |
|------|------|--------|
| Plasmid classification | MOB-suite | `27_mobsuite.sh` |
| Plasmid replicon identification | PlasmidFinder v2.1.6 | `29_plasmidfinder.sh` |

### 8. Association Analysis
| Step | Tool | Script |
|------|------|--------|
| Pan-genome gene association | Scoary v1.6.16 | `28_scoary.sh` |

### 9. Analyses Blocked by HPC Constraints
| Step | Tool | Reason |
|------|------|--------|
| Population structure | PopPUNK | glibc incompatibility on HPC |
| Functional annotation | eggNOG-mapper | SSL network restriction on HPC |

---

## Key Results

| Metric | Value |
|--------|-------|
| Total genomes | 234 |
| Sequence types identified | 20 |
| Dominant lineage | ST23 (n=23, 9.8%) |
| Core genome size | 3,154 genes |
| Total pan-genome | 23,699 genes |
| Hypervirulent isolates (score ≥4) | 142/234 (61%) |
| Convergent strains (high virulence + resistance) | 20 |
| ESBL gene (blaSHV-11) | 131/234 (56%) |
| Plasmid carriage | 214/234 (91.5%) |

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
| `mlst_env` | MLST, ABRicate, Scoary, PopPUNK |
| `amrfinder_env` | AMRFinderPlus |
| `kleborate_env` | Kleborate |
| `mobsuite_env` | MOB-suite, PlasmidFinder |
| `panaroo_env` | Panaroo, MAFFT |
| `iqtree_env` | IQ-TREE, Gubbins |
| `eggnog_env` | eggNOG-mapper (db download blocked) |

---

## Author

**Maida Jajja**
MSc Applied Bioinformatics, King's College London
GitHub: [@maidajajja](https://github.com/maidajajja)
