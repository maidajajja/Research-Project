#!/bin/bash

KRAKEN_DIR="/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/Kraken2"
OUTPUT="/scratch/users/k22017808/KP_Research_Project/05_Taxonomy/kraken2_summary.tsv"

echo -e "Sample\tKlebsiella_%\tKlebsiella_pneumoniae_%\tUnclassified_%" > $OUTPUT

for report in $KRAKEN_DIR/*.report; do
    sample=$(basename $report .report)
    klebsiella=$(awk -F'\t' '$NF ~ /Klebsiella$/ {print $1}' $report | head -1)
    kp=$(awk -F'\t' '$NF ~ /Klebsiella pneumoniae$/ {print $1}' $report | head -1)
    unclass=$(awk -F'\t' '$NF ~ /unclassified$/ {print $1}' $report | head -1)
    echo -e "${sample}\t${klebsiella:-0}\t${kp:-0}\t${unclass:-0}" >> $OUTPUT
done

echo "Done!"
