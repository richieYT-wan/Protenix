#!/bin/bash
# Sets up a local SAbDab-nano (nanobody) MMseqs2 database from a prebuilt FASTA.
# FASTA is fetched from GCS.
#
# Prerequisites:
#   - mmseqs2 installed and on PATH
#   - gcloud (or gsutil) for GCS download
#
# Usage:
#   bash scripts/setup_sabdab_nano_db.sh

set -euo pipefail

DB_DIR="$HOME/search_database/sabdab_nano"
FASTA="$DB_DIR/sabdab_nano.fasta"
GCS_FASTA="gs://em52-ab-develop-analytics-prod-f684/data/denovo-design/snakemake/data/01_raw/protenix-inputs/search_database/sabdab_nano/sabdab_nano.fasta"

mkdir -p "$DB_DIR"
cd "$DB_DIR"

##############################################################################
# Step 1: Fetch FASTA from GCS                                               #
##############################################################################
echo "=== Step 1/2: Fetching FASTA from GCS ==="

if [[ -f "$FASTA" ]]; then
    echo "[info] FASTA already present at $FASTA"
else
    if command -v gcloud >/dev/null 2>&1; then
        gcloud storage cp "$GCS_FASTA" "$FASTA"
    elif command -v gsutil >/dev/null 2>&1; then
        gsutil cp "$GCS_FASTA" "$FASTA"
    else
        echo "[error] Neither gcloud nor gsutil found. Install one, or manually place the FASTA at:"
        echo "        $FASTA"
        exit 1
    fi
fi

n_seqs=$(grep -c "^>" "$FASTA")
echo "[info] FASTA: $FASTA ($n_seqs sequences, $(du -h "$FASTA" | cut -f1))"

##############################################################################
# Step 2: Build MMseqs2 database                                             #
##############################################################################
echo ""
echo "=== Step 2/2: Building MMseqs2 database ==="
mmseqs createdb "$FASTA" "$DB_DIR/sabdab_nano_db"
mkdir -p "$DB_DIR/tmp"
mmseqs createindex "$DB_DIR/sabdab_nano_db" "$DB_DIR/tmp" --search-type 1
rm -rf "$DB_DIR/tmp"

echo ""
echo "=========================================================================="
echo "[SUCCESS] SAbDab-nano database ready at: $DB_DIR/sabdab_nano_db"
echo ""
echo "Test it with:"
echo "  mmseqs easy-search my_vhh.fasta $DB_DIR/sabdab_nano_db hits.m8 /tmp/mmseqs_tmp"
echo "=========================================================================="