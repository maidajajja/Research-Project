#!/bin/bash
#=============================================================================
# Download genome assemblies from NCBI using assembly accessions
# Author: Maida Jajja
# Date: 2025-02-17
# Purpose: Download 142 K. pneumoniae genome assemblies for cirrhosis project
#=============================================================================

set -e  # Exit on any error

# Paths
PROJECT_DIR="$HOME/kp_liver_project"
ACCESSION_FILE="${PROJECT_DIR}/lists/new_assembly_accessions.txt"
OUTPUT_DIR="${PROJECT_DIR}/genomes_fasta/all_genomes"
TEMP_DIR="${PROJECT_DIR}/genomes_fasta/temp_download"
LOG_FILE="${PROJECT_DIR}/logs/download_new_assemblies_$(date +%Y%m%d_%H%M%S).log"

# Logging function
log() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1" | tee -a "$LOG_FILE"
}

#=============================================================================
# Main workflow
#=============================================================================

log "========================================="
log "Starting genome download"
log "========================================="
log "Accessions file: $ACCESSION_FILE"
log "Output directory: $OUTPUT_DIR"
log ""

# Check dependencies
if ! command -v datasets &> /dev/null; then
    log "ERROR: NCBI 'datasets' CLI not found"
    log "Install with: conda install -c conda-forge ncbi-datasets-cli"
    exit 1
fi

# Check accession file exists
if [ ! -f "$ACCESSION_FILE" ]; then
    log "ERROR: Accession file not found: $ACCESSION_FILE"
    exit 1
fi

TOTAL=$(grep -c . "$ACCESSION_FILE")
log "Total accessions to download: $TOTAL"
log ""

# Create temp directory
mkdir -p "$TEMP_DIR"
cd "$TEMP_DIR"

#=============================================================================
# Step 1: Download genomes
#=============================================================================

log "[1/3] Downloading genomes from NCBI..."
datasets download genome accession \
    --inputfile "$ACCESSION_FILE" \
    --include genome \
    --filename genomes.zip \
    2>&1 | tee -a "$LOG_FILE"

if [ $? -ne 0 ]; then
    log "ERROR: Download failed"
    exit 1
fi

log "Download complete: genomes.zip"
log ""

#=============================================================================
# Step 2: Extract genomes
#=============================================================================

log "[2/3] Extracting genomes..."
unzip -q genomes.zip
log "Extraction complete"
log ""

#=============================================================================
# Step 3: Move FASTA files to final location
#=============================================================================

log "[3/3] Moving FASTA files to $OUTPUT_DIR..."

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Find and move all .fna files
FASTA_COUNT=0
for fasta in $(find ncbi_dataset/data -name "*.fna"); do
    # Get just the filename
    filename=$(basename "$fasta")
    
    # Move to output directory
    cp "$fasta" "$OUTPUT_DIR/"
    FASTA_COUNT=$((FASTA_COUNT + 1))
done

log "Moved $FASTA_COUNT FASTA files to $OUTPUT_DIR"
log ""

#=============================================================================
# Step 4: Cleanup
#=============================================================================

log "Cleaning up temporary files..."
cd "$PROJECT_DIR"
rm -rf "$TEMP_DIR"
log "Cleanup complete"
log ""

#=============================================================================
# Summary
#=============================================================================

FINAL_COUNT=$(ls -1 "$OUTPUT_DIR"/*.fna 2>/dev/null | wc -l)

log "========================================="
log "Download Complete!"
log "========================================="
log "Requested: $TOTAL genomes"
log "Downloaded: $FASTA_COUNT genomes"
log "Total in directory: $FINAL_COUNT genomes"
log "Location: $OUTPUT_DIR"
log "Log file: $LOG_FILE"
log "========================================="
