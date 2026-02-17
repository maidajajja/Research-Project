#!/bin/bash
# Look up Assembly accessions from BioSample IDs

BIOSAMPLES=(
  "SAMN41061928"
  "SAMN41061926"
  "SAMN41061929"
  "SAMN41061925"
  "SAMN41061924"
  "SAMN41061923"
  "SAMN48077203"
)

GENOME_IDS=(
  "573.80604"
  "573.80605"
  "573.80606"
  "573.80607"
  "573.80608"
  "573.80610"
  "573.80611"
  "573.81747"
)

echo "Looking up Assembly accessions for failed genomes..."
echo ""

for i in "${!BIOSAMPLES[@]}"; do
  biosample="${BIOSAMPLES[$i]}"
  genome_id="${GENOME_IDS[$i]}"
  
  echo "[$((i+1))/8] $genome_id (BioSample: $biosample)"
  
  # Use esearch to find assembly from biosample
  assembly=$(esearch -db biosample -query "$biosample" 2>/dev/null | \
             elink -target assembly 2>/dev/null | \
             esummary 2>/dev/null | \
             xtract -pattern DocumentSummary -element AssemblyAccession 2>/dev/null)
  
  if [ -n "$assembly" ]; then
    echo "  ✓ Found: $assembly"
  else
    echo "  ✗ No assembly found"
  fi
  echo ""
done
