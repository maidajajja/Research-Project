#!/usr/bin/env python3
"""
Download GenBank sequences and combine by genome ID
Author: Maida Jajja
Date: 2025-02-17
Purpose: Download 80 K. pneumoniae genomes using GenBank accessions (CP/JB numbers)
"""

import csv
import subprocess
import sys
import os
from pathlib import Path
from datetime import datetime

# Paths
PROJECT_DIR = Path.home() / "kp_liver_project"
INPUT_CSV = PROJECT_DIR / "lists" / "genbank_map.csv"
OUTPUT_DIR = PROJECT_DIR / "genomes_fasta" / "all_genomes"
TEMP_DIR = PROJECT_DIR / "genomes_fasta" / "temp_genbank_download"
LOG_FILE = PROJECT_DIR / "logs" / f"download_genbank_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

def log(message):
    """Write message to both console and log file"""
    timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
    log_msg = f"[{timestamp}] {message}"
    print(log_msg)
    with open(LOG_FILE, 'a') as f:
        f.write(log_msg + '\n')

def download_sequence(accession, output_file):
    """Download a single GenBank sequence using efetch"""
    try:
        cmd = [
            "efetch",
            "-db", "nuccore",
            "-id", accession,
            "-format", "fasta"
        ]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        
        if result.stdout and len(result.stdout) > 0:
            with open(output_file, 'w') as f:
                f.write(result.stdout)
            return True
        else:
            log(f"  WARNING: Empty response for {accession}")
            return False
    except subprocess.CalledProcessError as e:
        log(f"  ERROR downloading {accession}: {e}")
        return False

def main():
    log("=" * 60)
    log("Starting GenBank genome download")
    log("=" * 60)
    log(f"Input CSV: {INPUT_CSV}")
    log(f"Output directory: {OUTPUT_DIR}")
    log("")
    
    # Check dependencies
    try:
        subprocess.run(["efetch", "-version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        log("ERROR: EDirect 'efetch' not found")
        log("Install with: conda install -c bioconda entrez-direct")
        sys.exit(1)
    
    # Check input file
    if not INPUT_CSV.exists():
        log(f"ERROR: Input CSV not found: {INPUT_CSV}")
        sys.exit(1)
    
    # Create directories
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    TEMP_DIR.mkdir(parents=True, exist_ok=True)
    
    # Read CSV
    genomes = []
    with open(INPUT_CSV, 'r', encoding='utf-8-sig') as f:
        reader = csv.DictReader(f)
        for row in reader:
            genome_id = row['Genome ID']
            accessions = row['GenBank Accessions'].replace('"', '').split(',')
            genomes.append((genome_id, accessions))
    
    log(f"Found {len(genomes)} genomes to download")
    log("")
    
    # Download each genome
    success_count = 0
    fail_count = 0
    
    for i, (genome_id, accessions) in enumerate(genomes, 1):
        log(f"[{i}/{len(genomes)}] Processing {genome_id} ({len(accessions)} sequences)")
        
        # Download all sequences for this genome
        temp_files = []
        download_success = True
        
        for acc in accessions:
            acc = acc.strip()
            temp_file = TEMP_DIR / f"{acc}.fasta"
            
            if download_sequence(acc, temp_file):
                temp_files.append(temp_file)
                log(f"  ✓ Downloaded {acc}")
            else:
                log(f"  ✗ Failed to download {acc}")
                download_success = False
        
        # Combine sequences into one file
        if temp_files:
            combined_file = OUTPUT_DIR / f"{genome_id}.fasta"
            with open(combined_file, 'w') as outf:
                for temp_file in temp_files:
                    with open(temp_file, 'r') as inf:
                        outf.write(inf.read())
            
            log(f"  → Combined into {genome_id}.fasta")
            
            if download_success:
                success_count += 1
            else:
                fail_count += 1
                log(f"  ⚠ Warning: Some sequences missing for {genome_id}")
        else:
            log(f"  ✗ No sequences downloaded for {genome_id}")
            fail_count += 1
        
        log("")
    
    # Cleanup temp directory
    log("Cleaning up temporary files...")
    for file in TEMP_DIR.glob("*.fasta"):
        file.unlink()
    TEMP_DIR.rmdir()
    
    # Summary
    log("=" * 60)
    log("Download Complete!")
    log("=" * 60)
    log(f"Successfully downloaded: {success_count} genomes")
    log(f"Failed or incomplete: {fail_count} genomes")
    log(f"Output location: {OUTPUT_DIR}")
    log(f"Log file: {LOG_FILE}")
    log("=" * 60)

if __name__ == "__main__":
    main()
