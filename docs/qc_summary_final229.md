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
