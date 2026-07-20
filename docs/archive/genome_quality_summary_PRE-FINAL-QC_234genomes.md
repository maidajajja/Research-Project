# Genome Quality Assessment Summary
**Date:** 2026-02-19  
**Analyst:** Maida Jajja  
**Total Genomes:** 234

## Overview
Quality assessment performed using QUAST v5.2.0 on all downloaded K. pneumoniae genomes.

## Quality Categories

### Excellent Quality (73%, n=172)
- **Criteria:** ≤100 contigs, N50 >100kb, genome size 5-6 Mb
- **Status:** Ready for all downstream analyses

### Questionable Quality (13%, n=30)
- **Criteria:** 100-150 contigs OR N50 50-100kb
- **Status:** Acceptable for most analyses, may need review

### Poor Quality (14%, n=32)
- **Criteria:** >150 contigs OR N50 <50kb
- **Status:** Requires supervisor review before proceeding

## Poor Quality Genomes (>150 contigs OR N50 <50kb)

### Highly Fragmented (>250 contigs):
1. SRR16202829_assembled - 802 contigs, N50=206,112 bp
2. GCA_008630125.1_ASM863012v1_genomic - 456 contigs, N50=321,095 bp
3. SRR26938680_assembled - 287 contigs, N50=243,700 bp
4. GCA_016805585.1_ASM1680558v1_genomic - 283 contigs, N50=64,511 bp
5. GCA_000338255.1_KP1_contigs_1-272_genomic - 272 contigs, N50=43,309 bp

### Moderately Fragmented (150-250 contigs):
6. GCA_026623375.1_ASM2662337v1_genomic - 242 contigs
7. GCA_026624205.1_ASM2662420v1_genomic - 232 contigs
8. GCA_026624025.1_ASM2662402v1_genomic - 207 contigs
9. GCA_026624045.1_ASM2662404v1_genomic - 202 contigs
10. GCA_024136685.1_ASM2413668v1_genomic - 195 contigs
11. SRR16202849_assembled - 194 contigs
12. GCA_024136665.1_ASM2413666v1_genomic - 193 contigs
13. GCA_024136865.1_ASM2413686v1_genomic - 193 contigs
14. GCA_024136585.1_ASM2413658v1_genomic - 189 contigs
15. GCA_024136745.1_ASM2413674v1_genomic - 189 contigs
16. GCA_024136625.1_ASM2413662v1_genomic - 187 contigs
17. SRR16202839_assembled - 186 contigs
18. GCA_024136605.1_ASM2413660v1_genomic - 182 contigs
19. GCA_024136765.1_ASM2413676v1_genomic - 182 contigs
20. SRR26896786_assembled - 182 contigs
21. GCA_024136565.1_ASM2413656v1_genomic - 180 contigs
22. GCA_024136695.1_ASM2413669v1_genomic - 178 contigs
23. GCA_024136805.1_ASM2413680v1_genomic - 177 contigs
24. GCA_019249995.1_ASM1924999v1_genomic - 174 contigs
25. GCA_024136645.1_ASM2413664v1_genomic - 171 contigs
26. GCA_024136905.1_ASM2413690v1_genomic - 171 contigs
27. GCA_025762435.1_ASM2576243v1_genomic - 171 contigs
28. GCA_025762675.1_ASM2576267v1_genomic - 169 contigs
29. GCA_002831525.1_ASM283152v1_genomic - 168 contigs
30. GCA_025762465.1_ASM2576246v1_genomic - 164 contigs
31. GCA_025762455.1_ASM2576245v1_genomic - 162 contigs
32. GCA_026156485.1_ASM2615648v1_genomic - 153 contigs

## SRA Assembly Performance
Of the 10 SRA-assembled genomes:
- **Good quality:** 5 genomes (50%)
- **Poor quality:** 5 genomes (50%)
  - SRR16202829 (802 contigs)
  - SRR26938680 (287 contigs)
  - SRR16202849 (194 contigs)
  - SRR16202839 (186 contigs)
  - SRR26896786 (182 contigs)

## Overall Statistics
- **Mean genome size:** 5.6 Mb (range: 4.8-6.2 Mb)
- **Mean GC content:** 57.1% (correct for K. pneumoniae)
- **Mean N50:** 294 kb (median: 301 kb)
- **Mean contig count:** 76 (median: 45)

## Recommendations
1. **Proceed with all 234 genomes** for Prokka annotation
2. **Discuss exclusion criteria** with supervisor before comparative genomics
3. **Likely threshold:** Exclude genomes with >200 contigs OR N50 <50kb
4. **Flagged for review:** 32 poor quality genomes listed above

## Next Steps
1. Prokka annotation (all 234 genomes)
2. Taxonomic confirmation (GTDB-Tk + Kraken2/FastANI)
3. Quality threshold discussion with supervisor
4. Proceed with comparative genomics on filtered dataset

## Files Generated
- `quast_results/report.html` - Interactive quality report
- `quast_results/report.txt` - Text summary
- `quast_results/report.tsv` - Spreadsheet format for further analysis
