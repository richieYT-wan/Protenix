#!/usr/bin/env bash
#
# Unified installation script for the forked Protenix repository.
# This script sets up the conda environment, installs the package,
# downloads required binaries and data, and sets up the necessary tools.
#

set -euo pipefail

# --- Configuration Variables ---
# (Edit these to match your local setup)

# Name for the Conda environment
ENV_NAME="protenix"

# Base directory for installing external binaries like MMseqs2
INSTALL_PREFIX="$HOME/protenix_bin"

# Directory to store model weights/checkpoints
WEIGHTS_DIR="$HOME/checkpoint"

# Directory to store search databases (PDB, VHH, etc.)
DB_DIR="$HOME/protenix_dbs"

# GCS bucket for downloading weights and databases.
# This should be the base path containing 'checkpoint' and 'search_database' folders.
GCS_BUCKET="gs://em52-ab-develop-analytics-prod-f684/data/denovo-design/snakemake/data/01_raw/protenix-inputs"

# Version for MMseqs2 to be built from source
MMSEQS_VERSION="16-747c6"

# --- Script Start ---

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

echo "=== Protenix Unified Installer ==="
echo "  Repo Dir:       $REPO_DIR"
echo "  Conda Env:      $ENV_NAME"
echo "  Install Prefix: $INSTALL_PREFIX"
echo "  Weights Dir:    $WEIGHTS_DIR"
echo "  Database Dir:   $DB_DIR"
echo "  GCS Bucket:     $GCS_BUCKET"
echo "=================================="

# --- Step 1: Conda Environment Setup ---
echo "=== Step 1: Setting up Conda environment: '$ENV_NAME' ==="

# Use mamba if available, otherwise fall back to conda
CONDA_CMD="conda"
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Mamba detected, using it for faster environment creation."
fi

# Check if environment already exists
if "$CONDA_CMD" env list | grep -q "^${ENV_NAME}\s"; then
    echo "Conda environment '$ENV_NAME' already exists. Skipping creation."
else
    echo "Creating Conda environment from 'protenix_env.yaml'..."
    # Logic from: protenix_env.yaml, setup.sh
    "$CONDA_CMD" env create -f protenix_env.yaml
    echo "Environment created."
fi

# Activate the environment for the rest of the script
# shellcheck disable=SC1091
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "Conda environment activated."
echo "Python version: $(python --version)"
echo "Conda prefix: $CONDA_PREFIX"

# --- Step 2: Install Protenix Package ---
echo "=== Step 2: Installing Protenix (non-editable) ==="
# Logic from: setup.sh, run_setup.sh
pip install .
echo "Protenix package installed."

# --- Step 3: Install CUDA Toolkit ---
echo "=== Step 3: Installing CUDA Toolkit for PyTorch ==="
# Logic from: setup.sh, run_setup.sh
if ! python -c "import torch; print(torch.version.cuda)" &> /dev/null; then
    echo "PyTorch with CUDA support not found. Skipping CUDA Toolkit installation."
    echo "Please ensure you have the correct PyTorch version installed for your GPU."
else
    CUDA_VERSION="$(python -c 'import torch; print(torch.version.cuda)')"
    echo "Detected CUDA version $CUDA_VERSION from PyTorch."
    pip install "cuda-toolkit==${CUDA_VERSION}"
    echo "CUDA Toolkit installed."
fi

# --- Step 4: Install HMMER ---
echo "=== Step 4: Verifying HMMER installation ==="
# Logic from: protenix_env.yaml
if command -v hmmsearch &> /dev/null; then
    echo "HMMER (hmmsearch) found at: $(command -v hmmsearch)"
else
    echo "WARNING: HMMER not found. It should have been installed by conda."
    echo "You may need to install it manually: conda install -c bioconda hmmer"
fi

# --- Step 5: Install MMseqs2 from Source ---
echo "=== Step 5: Installing MMseqs2 v${MMSEQS_VERSION} from source ==="
# Logic from: setup.sh
MMSEQS_BIN_PATH="${INSTALL_PREFIX}/mmseqs2/bin/mmseqs"
if [[ -f "$MMSEQS_BIN_PATH" ]]; then
    echo "MMseqs2 already found at $MMSEQS_BIN_PATH. Skipping installation."
else
    echo "MMseqs2 not found. Proceeding with installation..."
    mkdir -p "$INSTALL_PREFIX"
    cd "$INSTALL_PREFIX"

    echo "Downloading MMseqs2 source..."
    wget -q "https://github.com/soedinglab/MMseqs2/archive/refs/tags/${MMSEQS_VERSION}.tar.gz"
    tar xzf "${MMSEQS_VERSION}.tar.gz"

    MMSEQS_SRC_DIR="MMseqs2-${MMSEQS_VERSION}"
    cd "$MMSEQS_SRC_DIR"

    echo "Building MMseqs2... (this may take a few minutes)"
    mkdir -p build && cd build
    cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX}/mmseqs2" ..
    make -j"$(nproc)"
    make install
    echo "MMseqs2 build complete."

    cd "$REPO_DIR" # Return to repo root
fi

# Symlink the binary into the conda environment to make it available on the PATH
if [[ -f "${CONDA_PREFIX}/bin/mmseqs" ]]; then
    echo "MMseqs2 symlink already exists in conda env."
else
    echo "Creating symlink for mmseqs in conda environment..."
    ln -sf "$MMSEQS_BIN_PATH" "${CONDA_PREFIX}/bin/mmseqs"
fi

echo "MMseqs2 is available at: $(command -v mmseqs)"
"$(command -v mmseqs)" version

# --- Step 6: Download Model Weights and Databases ---
echo "=== Step 6: Downloading model weights and search databases ==="

# Download Protenix-V2 weights
echo "Downloading weights to ${WEIGHTS_DIR}..."
# Logic from: setup.sh
mkdir -p "$WEIGHTS_DIR"
gsutil -m rsync -r "${GCS_BUCKET}/checkpoint/" "${WEIGHTS_DIR}/"
echo "Weights download complete. Contents:"
ls -1 "$WEIGHTS_DIR"

# Download search databases (VHH, PDB seqres)
echo "Downloading search databases to ${DB_DIR}..."
# Logic from: setup.sh
mkdir -p "$DB_DIR"
gsutil -m rsync -r "${GCS_BUCKET}/search_database/" "${DB_DIR}/"
echo "Database download complete. Contents:"
ls -1R "$DB_DIR"

# --- Finalization ---
echo ""
echo "✅ Protenix setup complete!"
echo ""
echo "To activate the environment, run:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Paths:"
echo "  - Model Weights:   $WEIGHTS_DIR"
echo "  - Databases:       $DB_DIR"
echo "  - MMseqs2 Binary:  $MMSEQS_BIN_PATH"
echo ""
echo "To verify the installation, you can run:"
echo "  protenix pred --help"
echo ""