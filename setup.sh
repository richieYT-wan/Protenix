#!/usr/bin/env bash
# Run-once setup script for our fork of Protenix.

set -euxo pipefail

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

ENV_NAME="protenix_test"
GCS_BASE="gs://em52-ab-develop-analytics-prod-f684/data/denovo-design/snakemake/data/01_raw/protenix-inputs"
CHECKPOINT_DIR="$HOME/checkpoint"
SEARCH_DB_DIR="$HOME/search_database"

# ---- Conda env ----
source "$(conda info --base)/etc/profile.d/conda.sh"
conda env create -f protenix_env.yaml
conda activate "$ENV_NAME"

# ---- Python install ----
pip install -e .
CUDA_VERSION="$(python -c 'import torch; print(torch.version.cuda)')"
pip install "cuda-toolkit==${CUDA_VERSION}"
protenix pred --help

# ---- Weights ----
echo "Downloading weights to ${CHECKPOINT_DIR}"
mkdir -p "$CHECKPOINT_DIR"
gsutil -m cp -r "${GCS_BASE}/checkpoint/*" "${CHECKPOINT_DIR}/"
ls "$CHECKPOINT_DIR"

# ---- Search databases ----
# Contains:
#   - pdb_seqres_2022_09_28.fasta          (template search)
#   - sabdab_nano/sabdab_nano_db*          (VHH MMseqs2 DB, prebuilt)
echo "Downloading search databases to ${SEARCH_DB_DIR}"
mkdir -p "$SEARCH_DB_DIR"
gsutil -m cp -r "${GCS_BASE}/search_database/*" "${SEARCH_DB_DIR}/"
ls -R "$SEARCH_DB_DIR"

# Sanity: check the paths make_json_from_csv.py expects
test -f "${SEARCH_DB_DIR}/sabdab_nano/sabdab_nano_db.dbtype"
test -f "${SEARCH_DB_DIR}/pdb_seqres_2022_09_28.fasta"

# ---- MMseqs2 + ColabFold DBs ----
MMSEQS_VERSION="16-747c6"
MMSEQS_DIR="$HOME/mmseqs2"
COLABFOLD_DB_DIR="$HOME/colabfold_db"

mkdir -p "$MMSEQS_DIR"
cd "$MMSEQS_DIR"
wget -q "https://github.com/soedinglab/MMseqs2/archive/refs/tags/${MMSEQS_VERSION}.tar.gz"
tar xzf "${MMSEQS_VERSION}.tar.gz"
cd "MMseqs2-${MMSEQS_VERSION}"
mkdir -p build && cd build
cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
make -j"$(nproc)"
make install

# Make mmseqs available on PATH for downstream scripts (make_json_from_csv.py calls it)
ln -sf "${MMSEQS_DIR}/MMseqs2-${MMSEQS_VERSION}/build/bin/mmseqs" \
       "${CONDA_PREFIX}/bin/mmseqs"
mmseqs version

cd "$MMSEQS_DIR"
git clone https://github.com/sokrypton/ColabFold.git
cd ColabFold
mkdir -p "$COLABFOLD_DB_DIR"
MMSEQS_NO_INDEX=1 ./setup_databases.sh "$COLABFOLD_DB_DIR"

echo "Setup complete."
echo "  Checkpoints:   ${CHECKPOINT_DIR}"
echo "  Search DBs:    ${SEARCH_DB_DIR}"
echo "  ColabFold DBs: ${COLABFOLD_DB_DIR}"
echo "  MMseqs2:       $(command -v mmseqs)"