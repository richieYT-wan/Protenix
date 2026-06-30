#!/bin/bash
# Sets up a local SAbDab-nano (nanobody) MMseqs2 database
# SAbDab-nano: https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/nano/
#
# Prerequisites:
#   - mmseqs2 installed and on PATH
#   - python3 with biopython (for parsing PDB/FASTA)
#
# Usage:
#   bash scripts/setup_sabdab_nano_db.sh [output_dir]

set -euo pipefail

DB_DIR="$HOME/search_database/sabdab_nano/"
mkdir -p "$DB_DIR"
cd "$DB_DIR"

##############################################################################
# Step 1: Download SAbDab-nano summary + sequences                           #
##############################################################################
echo "=== Step 1/3: Downloading SAbDab-nano ==="


# Download the all-nano structures archive (contains FASTA-derivable sequences)
# Option A: download the full PDB archive (large but easy)

ARCHIVE_URL="https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/archive/all_nano/"
wget -nc -O sabdab_nano_archive.zip "$ARCHIVE_URL" || {
    echo "[warn] Archive download failed. Falling back to per-PDB download from summary."
}
##############################################################################
# Step 2: Extract sequences into a FASTA                                     #
##############################################################################
echo ""
echo "=== Step 2/3: Extracting sequences ==="
cd $DB_DIR

python3 - <<'PYEOF'
import csv
import gzip
import re
import zipfile
from pathlib import Path

out_fasta = Path("sabdab_nano.fasta")
seqs_seen = set()
count = 0

# Try to extract from the archive if present
archive = Path("sabdab_nano_archive.zip")
if archive.exists() and archive.stat().st_size > 0:
    print(f"[info] Extracting from {archive}")
    with zipfile.ZipFile(archive) as zf, open(out_fasta, "w") as fout:
        for name in zf.namelist():
            if not name.endswith((".pdb", ".pdb.gz", ".cif", ".cif.gz")):
                continue
            try:
                with zf.open(name) as fh:
                    data = fh.read()
                    if name.endswith(".gz"):
                        data = gzip.decompress(data)
                    text = data.decode("utf-8", errors="ignore")

                # Pull SEQRES records out of PDB (simple parser)
                chains = {}
                for line in text.splitlines():
                    if line.startswith("SEQRES"):
                        chain_id = line[11]
                        residues = line[19:].split()
                        chains.setdefault(chain_id, []).extend(residues)

                # Convert 3-letter to 1-letter and filter by length
                three_to_one = {
                    'ALA':'A','CYS':'C','ASP':'D','GLU':'E','PHE':'F','GLY':'G',
                    'HIS':'H','ILE':'I','LYS':'K','LEU':'L','MET':'M','ASN':'N',
                    'PRO':'P','GLN':'Q','ARG':'R','SER':'S','THR':'T','VAL':'V',
                    'TRP':'W','TYR':'Y'
                }
                for chain_id, residues in chains.items():
                    seq = ''.join(three_to_one.get(r, 'X') for r in residues)
                    seq = seq.replace('X', '')  # drop unknowns
                    # VHH typical length 100-150
                    if 90 < len(seq) < 160 and seq not in seqs_seen:
                        seqs_seen.add(seq)
                        pdb_id = Path(name).stem.split('.')[0]
                        fout.write(f">{pdb_id}_{chain_id}\n{seq}\n")
                        count += 1
            except Exception as e:
                print(f"[warn] Skipping {name}: {e}")
else:
    print("[error] No archive found and no fallback configured.")
    print("[error] You may need to manually download from:")
    print("        https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabdab/nano/")
    exit(1)

print(f"[done] Wrote {count} unique VHH sequences to {out_fasta}")
PYEOF

##############################################################################
# Step 3: Build MMseqs2 database                                             #
##############################################################################
echo ""
echo "=== Step 3/3: Building MMseqs2 database ==="
mmseqs createdb ${DB_DIR}sabdab_nano.fasta sabdab_nano_db
mkdir -p tmp
mmseqs createindex sabdab_nano_db tmp --search-type 1
rm -rf tmp

echo ""
echo "=========================================================================="
echo "[SUCCESS] SAbDab-nano database ready at: ${DB_DIR}sabdab_nano_db"
echo ""
echo "Test it with:"
echo "  mmseqs easy-search my_vhh.fasta ${DB_DIR}sabdab_nano_db hits.m8 /tmp/mmseqs_tmp"
echo "=========================================================================="