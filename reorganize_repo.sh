#!/usr/bin/env bash
# Reorganise Research-Project repo into a distinction-level structure.
# Run this from the ROOT of your local clone of github.com/maidajajja/Research-Project
# (e.g. on the HPC or your laptop, wherever you have git push access).
#
# What it does:
#   - Creates staged pipeline folders under scripts/
#   - Moves data/env/doc files into data/, environment/, docs/
#   - Drops 5 files confirmed as exact byte-for-byte duplicates (annotate batch scripts
#     committed twice under different number prefixes)
#   - Archives the stale 234-genome QC summary and adds an accurate 229-genome one
#   - Replaces README.md with a version matching the new structure
#   - Creates FINAL_FIGURES/ (empty - see the separate figure-copy commands)
#
# It does NOT commit or push. Review with `git status` / `git diff --stat` first,
# then commit and push yourself.

set -euo pipefail

if [ ! -f "README.md" ] || [ ! -d "scripts" ]; then
  echo "Run this from the repo root (README.md and scripts/ must be present here)." >&2
  exit 1
fi

mkdir -p scripts/01_data_acquisition scripts/02_annotation scripts/03_quality_control \
         scripts/04_typing_and_amr scripts/05_phylogenetics scripts/06_plasmid_analysis \
         scripts/07_pangenome_association scripts/figures environment data docs docs/archive \
         FINAL_FIGURES

echo "== 1. Data acquisition =="
git mv scripts/01_download_assemblies.sh scripts/01_data_acquisition/download_assemblies.sh
git mv scripts/02_download_genbank_isolates.py scripts/01_data_acquisition/download_genbank_isolates.py
git mv scripts/03_make_all_genomes.sh scripts/01_data_acquisition/make_all_genomes.sh
git mv scripts/04_download_new_assemblies.sh scripts/01_data_acquisition/download_new_assemblies.sh
git mv scripts/05_download_genbank_genomes.py scripts/01_data_acquisition/download_genbank_genomes.py
git mv scripts/06_sra_assembly_pipeline.sh scripts/01_data_acquisition/sra_assembly_pipeline.sh
git mv scripts/lookup_assemblies.sh scripts/01_data_acquisition/lookup_assemblies.sh

echo "== 2. Annotation (dropping 5 confirmed exact duplicates: 08-12_*_WORKED.sh) =="
git mv scripts/01_annotate_batch1_WORKED.sh scripts/02_annotation/annotate_batch1.sh
git mv scripts/02_annotate_batch2_WORKED.sh scripts/02_annotation/annotate_batch2.sh
git mv scripts/03_annotate_batch3_WORKED.sh scripts/02_annotation/annotate_batch3.sh
git mv scripts/04_annotate_batch4_WORKED.sh scripts/02_annotation/annotate_batch4.sh
git mv scripts/05_annotate_fasta_WORKED.sh scripts/02_annotation/annotate_fasta.sh
git rm -q scripts/08_annotate_batch1_WORKED.sh scripts/09_annotate_batch2_WORKED.sh \
          scripts/10_annotate_batch3_WORKED.sh scripts/11_annotate_batch4_WORKED.sh \
          scripts/12_annotate_fasta_WORKED.sh
git mv scripts/23_bakta_sgh10.sh scripts/02_annotation/bakta_sgh10.sh
git mv scripts/25_bakta_kp126.sh scripts/02_annotation/bakta_kp126.sh
git mv scripts/ANNOTATION_README.txt scripts/02_annotation/ANNOTATION_README.txt

echo "== 3. Quality control =="
git mv scripts/07_quast_all_genomes.sh scripts/03_quality_control/quast_all_genomes.sh
git mv scripts/15_kraken2_db_download.sh scripts/03_quality_control/kraken2_db_download.sh
git mv scripts/16_kraken2_classify.sh scripts/03_quality_control/kraken2_classify.sh
git mv scripts/22_kraken2_new2.sh scripts/03_quality_control/kraken2_new2.sh
git mv scripts/24_kraken2_kp126.sh scripts/03_quality_control/kraken2_kp126.sh
git mv scripts/kraken2_summary.sh scripts/03_quality_control/kraken2_summary.sh
git mv scripts/17_fastani_classify.sh scripts/03_quality_control/fastani_classify.sh

echo "== 4. Typing, AMR and virulence =="
git mv scripts/25_kleborate.sh scripts/04_typing_and_amr/kleborate.sh
git mv scripts/23_amrfinder_final229.sh scripts/04_typing_and_amr/amrfinder_final229.sh
git mv scripts/24_abricate_final229.sh scripts/04_typing_and_amr/abricate_final229.sh

echo "== 5. Phylogenetics =="
git mv scripts/26_panaroo_final229.sh scripts/05_phylogenetics/panaroo_final229.sh
git mv scripts/27_iqtree_final229.sh scripts/05_phylogenetics/iqtree_final229.sh
git mv scripts/26_gubbins_v2.sh scripts/05_phylogenetics/gubbins_v2.sh

echo "== 6. Plasmid analysis =="
git mv scripts/27_mobsuite.sh scripts/06_plasmid_analysis/mobsuite.sh
git mv scripts/27_mobsuite_force.sh scripts/06_plasmid_analysis/mobsuite_force.sh
git mv scripts/27b_mobsuite_remaining.sh scripts/06_plasmid_analysis/mobsuite_remaining.sh
git mv scripts/27c_mobsuite_fix.sh scripts/06_plasmid_analysis/mobsuite_fix.sh
git mv scripts/29_plasmidfinder.sh scripts/06_plasmid_analysis/plasmidfinder.sh
git mv scripts/28_amr_plasmid_new2.sh scripts/06_plasmid_analysis/amr_plasmid_join.sh

echo "== 7. Pan-genome association =="
git mv scripts/28_scoary_final229.sh scripts/07_pangenome_association/scoary_final229.sh

echo "== 8. Figures =="
git mv scripts/Fig2_AMR_heatmap_v19.R scripts/figures/Fig2_AMR_heatmap_v19.R
git mv scripts/Fig3a_virulence_heatmap_v2.R scripts/figures/Fig3a_virulence_heatmap_v2.R
git mv scripts/Fig3b_virulence_bubble_v5.R scripts/figures/Fig3b_virulence_bubble_v5.R
git mv scripts/Fig4_v5.R scripts/figures/Fig4_v5.R
git mv scripts/Fig5a_v4.R scripts/figures/Fig5a_v4.R
git mv scripts/Fig6_pangenome_v3.R scripts/figures/Fig6_pangenome_v3.R
git mv scripts/Fig6c_pangenome_heatmap_v4.R scripts/figures/Fig6c_pangenome_heatmap_v4.R
git mv scripts/Fig7_convergence_v3.R scripts/figures/Fig7_convergence_v3.R
git mv scripts/Fig8_convergent_plasmid_v4.R scripts/figures/Fig8_convergent_plasmid_v4.R
git mv scripts/Fig8b_convergent_integrated_v1.R scripts/figures/Fig8b_convergent_integrated_v1.R
git mv scripts/FigS1_FastANI_QC.R scripts/figures/FigS1_FastANI_QC.R
git mv scripts/FigS1_Kraken2.R scripts/figures/FigS1_Kraken2.R
git mv scripts/FigS2_FastANI_v3.R scripts/figures/FigS2_FastANI_v3.R
git mv scripts/FigS2_SNP_heatmap.R scripts/figures/FigS2_SNP_heatmap.R
git mv scripts/34_run_all_plots.sh scripts/figures/run_all_plots.sh
git mv scripts/36_publication_plots_v2.R scripts/figures/publication_plots_v2.R

echo "== 9. Environment files =="
git mv environment.yml environment/environment.yml
git mv amrfinder_environment.yml environment/amrfinder_environment.yml
git mv iqtree_environment.yml environment/iqtree_environment.yml
git mv kleborate_environment.yml environment/kleborate_environment.yml
git mv myRenv_environment.yml environment/myRenv_environment.yml
git mv panaroo_environment.yml environment/panaroo_environment.yml

echo "== 10. Data files =="
git mv genomes_master.csv data/genomes_master.csv
git mv lists data/lists

echo "== 11. Archive the stale QC summary =="
git mv genome_quality_summary.md docs/archive/genome_quality_summary_PRE-FINAL-QC_234genomes.md

echo "== 12. Add new docs and README, and a FINAL_FIGURES placeholder =="
cat > docs/qc_summary_final229.md << 'QC_EOF'
# Genome Quality Control Summary - Final Dataset (229 genomes)

This supersedes `docs/archive/genome_quality_summary_PRE-FINAL-QC_234genomes.md`, which reflected an
earlier, provisional QC pass (234 genomes, contig/N50-bucket criteria) before the final exclusion
and addition decisions below. This file matches Methods 2.1–2.2 and Results 3.1 of the dissertation.

## Final dataset size: 229 genomes

## Exclusions (8 total)

| # | Genome(s) | Reason | Evidence |
|---|-----------|--------|----------|
| 1 | GCA_009079735.1 | Species misidentification (Salmonella enterica, not K. pneumoniae) | FastANI 80–82% vs cohort, Kraken2 species call = S. enterica, GC 52.00% vs cohort mean 57.2% |
| 2 | GCA_000338255.1 (hvKP1, 2013 draft) | Assembly fragmentation (>200 contigs) | QUAST, 99th-percentile cohort threshold |
| 3 | GCA_016805585.1 | Assembly fragmentation (>200 contigs) | QUAST |
| 4 | GCA_002268655.1 | Assembly fragmentation (>200 contigs) | QUAST |
| 5 | GCA_020526065.1 (strain BJ1) | Assembly fragmentation (>200 contigs) | QUAST |
| 6 | GCA_002027655.1 | CDS-density/collapsed-repeat artefact | 12,812 predicted CDS vs cohort mean 5,198; passed QUAST contiguity thresholds (11 contigs, N50 1.1 Mbp) |
| 7 | One of FK3038 pair (GCA_025762655.1 / GCA_041198405.1) | Duplicate biological isolate | FastANI 100% |
| 8 | One of KP6 pair (GCA_002831525.1 / GCA_026623445.1) | Duplicate biological isolate | FastANI 99.96% |

## Additions (2 total)

| Genome | Reason |
|--------|--------|
| SGH10 (GCA_002813595.1) | Present in project metadata but absent from initial download; retrieved and QC-verified |
| KP126 (GCF_046599435.1) | Present in project metadata but absent from initial download; retrieved and QC-verified |

## Note on the earlier (archived) QC summary

The archived `genome_quality_summary_PRE-FINAL-QC_234genomes.md` used a provisional contig/N50-bucket
classification on 234 genomes, before the contaminant, duplicate, and CDS-density issues above were
identified. It is retained for audit history only. **Do not cite it as the current QC methodology.**
The 8 exclusions and 2 additions above are the figures used throughout the dissertation.
QC_EOF

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

| Figure | Script | Description |
|--------|--------|--------------|
| Fig1 | `scripts/05_phylogenetics/iqtree_final229.sh` + iTOL | Core-genome phylogeny |
| Fig2 | `scripts/figures/Fig2_AMR_heatmap_v19.R` | AMR gene heatmap by ST |
| Fig3a | `scripts/figures/Fig3a_virulence_heatmap_v2.R` | Virulence loci heatmap (supplementary) |
| Fig3b | `scripts/figures/Fig3b_virulence_bubble_v5.R` | Virulence prevalence bubble plot |
| Fig4 | `scripts/figures/Fig4_v5.R` | Plasmid replicon prevalence by ST |
| Fig5a | `scripts/figures/Fig5a_v4.R` | ST distribution over time |
| Fig6 | `scripts/figures/Fig6_pangenome_v3.R` | Pan-genome composition and frequency |
| Fig6c | `scripts/figures/Fig6c_pangenome_heatmap_v4.R` | Pan-genome accumulation curve |
| Fig7 | `scripts/figures/Fig7_convergence_v3.R` | MDR-hypervirulent convergence bubble plot |
| Fig8 | `scripts/figures/Fig8_convergent_plasmid_v4.R` | Convergent isolate heatmap |
| Fig8b | `scripts/figures/Fig8b_convergent_integrated_v1.R` | Integrated convergence figure |
| FigS1 | `scripts/figures/FigS1_FastANI_QC.R`, `FigS1_Kraken2.R` | QC supplementary figures |
| FigS2 | `scripts/figures/FigS2_FastANI_v3.R`, `FigS2_SNP_heatmap.R` | FastANI / SNP supplementary heatmaps |

`scripts/figures/run_all_plots.sh` and `publication_plots_v2.R` regenerate the full figure set in one pass.

---

## FINAL_FIGURES/

This folder holds the publication-quality PNG/PDF exports of every figure above, generated on the
KCL CREATE HPC cluster (`/scratch/users/k22017808/KP_Research_Project/FINAL_FIGURES/`). Copy the
contents of that HPC folder here and commit them before submission.

---

## HPC
Cluster: KCL CREATE HPC (hpc.create.kcl.ac.uk)
Partition: msc_appbio
Conda environments: myRenv (R packages), kleborate_env (Kleborate) - see `environment/` for full YAMLs
README_EOF

touch FINAL_FIGURES/.gitkeep
echo "Copy the real PNG/PDF figures here from the HPC (see figure-copy commands) - this file only reserves the empty folder in git." > FINAL_FIGURES/README.md

git add -A

echo ""
echo "Done. Review with:  git status   and   git diff --cached --stat"
echo "Then commit, e.g.:"
echo "  git commit -m \"Reorganise repo: stage scripts by pipeline step, fix duplicate numbering, tidy data/docs/env, add accurate final QC summary\""
echo "  git push"
