#!/usr/bin/bash
#
# Unified installation script for the forked Protenix repository.
# This script sets up the conda environment, installs the package,
# downloads required binaries and data, and sets up the necessary tools.
#
# Usage:
#   ./install.sh              # Full install (env + package + binaries + weights + DBs)
#   ./install.sh --no-db      # Skip database download (weights still downloaded)
#   ./install.sh -h | --help  # Show usage
#

set -euo pipefail

# --- Argument Parsing ---
SKIP_DB=0

usage() {
    cat <<EOF
Usage: $(basename "$0") [OPTIONS]

Options:
  --no-db       Skip downloading the large ColabFold database (colabfold_db.tar).
                The VHH search database is still downloaded.
  -h, --help    Show this help message and exit
EOF
}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --no-db)
            SKIP_DB=1
            shift
            ;;
        -h|--help)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option: $1" >&2
            usage
            exit 1
            ;;
    esac
done

# --- Configuration Variables ---
ENV_NAME="protenix"
INSTALL_PREFIX="$HOME/protenix_bin"
WEIGHTS_DIR="$HOME/checkpoint"
DB_DIR="$HOME/protenix_dbs"
GCS_BUCKET="gs://em52-ab-develop-analytics-prod-f684/data/denovo-design/snakemake/data/01_raw/protenix-inputs"
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
echo "  Skip DB:        $([[ $SKIP_DB -eq 1 ]] && echo yes || echo no)"
echo "=================================="

# --- Step 1: Conda Environment Setup ---
echo "=== Step 1: Setting up Conda environment: '$ENV_NAME' ==="

CONDA_CMD="conda"
if command -v mamba &> /dev/null; then
    CONDA_CMD="mamba"
    echo "Mamba detected, using it for faster environment creation."
fi

if "$CONDA_CMD" env list | grep -q "^${ENV_NAME}\s"; then
    echo "Conda environment '$ENV_NAME' already exists. Skipping creation."
else
    echo "Creating Conda environment from 'protenix_env.yaml'..."
    "$CONDA_CMD" env create --quiet -f protenix_env.yaml
    echo "Environment created."
fi

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$ENV_NAME"

echo "Conda environment activated."
echo "Python version: $(python --version)"
echo "Conda prefix: $CONDA_PREFIX"

# --- Step 2: Install Protenix Package ---
echo "=== Step 2: Installing Protenix (non-editable) ==="
pip install .
echo "Protenix package installed."

# --- Step 3: Install real CUDA Toolkit (nvcc + headers) into conda env ---
echo "=== Step 3: Installing CUDA Toolkit (nvcc + headers) ==="

if ! python -c "import torch" &> /dev/null; then
    echo "ERROR: PyTorch not importable. Aborting."
    exit 1
fi

CUDA_VERSION="$(python -c 'import torch; print(torch.version.cuda)')"
echo "PyTorch was built against CUDA ${CUDA_VERSION}."

# Remove the useless stub pip package if present
pip uninstall -y cuda-toolkit 2>/dev/null || true

if command -v nvcc &> /dev/null && find "$CONDA_PREFIX" -name "cuda_runtime_api.h" -print -quit | grep -q .; then
    echo "nvcc and CUDA headers already present. Skipping toolkit install."
else
    echo "Installing NVIDIA CUDA Toolkit ${CUDA_VERSION} (with headers) into conda env..."
    "$CONDA_CMD" install -y -n "$ENV_NAME" \
        -c "nvidia/label/cuda-${CUDA_VERSION}.0" \
        cuda-toolkit cuda-cudart-dev cuda-nvcc cuda-cccl
fi

# Set env vars for the current session
export CUDA_HOME="$CONDA_PREFIX"
export PATH="$CUDA_HOME/bin:$PATH"
export CPATH="$CONDA_PREFIX/targets/x86_64-linux/include:${CPATH:-}"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"

# Persist via conda activate hook so future 'conda activate protenix' sets them automatically
ACTIVATE_DIR="${CONDA_PREFIX}/etc/conda/activate.d"
mkdir -p "$ACTIVATE_DIR"
cat > "${ACTIVATE_DIR}/cuda_home.sh" <<'EOF'
export CUDA_HOME="$CONDA_PREFIX"
export PATH="$CUDA_HOME/bin:$PATH"
export CPATH="$CONDA_PREFIX/targets/x86_64-linux/include:${CPATH:-}"
export LD_LIBRARY_PATH="$CUDA_HOME/lib64:${LD_LIBRARY_PATH:-}"
EOF

echo "CUDA_HOME: $CUDA_HOME"
echo "CPATH:     $CPATH"
nvcc --version | tail -1

# --- Step 4: Install HMMER ---
echo "=== Step 4: Verifying HMMER installation ==="
if command -v hmmsearch &> /dev/null; then
    echo "HMMER (hmmsearch) found at: $(command -v hmmsearch)"
else
    echo "WARNING: HMMER not found. It should have been installed by conda."
    echo "You may need to install it manually: conda install -c bioconda hmmer"
fi

# --- Step 5: Install MMseqs2 from Source ---
echo "=== Step 5: Installing MMseqs2 v${MMSEQS_VERSION} from source ==="
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

    cd "$REPO_DIR"
fi

if [[ -L "${CONDA_PREFIX}/bin/mmseqs" ]]; then
    echo "MMseqs2 symlink already exists in conda env. Refreshing..."
    ln -sf "$MMSEQS_BIN_PATH" "${CONDA_PREFIX}/bin/mmseqs"
elif [[ -e "${CONDA_PREFIX}/bin/mmseqs" ]]; then
    echo "Existing (non-symlink) mmseqs found in conda env. Replacing with symlink to built binary..."
    rm -f "${CONDA_PREFIX}/bin/mmseqs"
    ln -sf "$MMSEQS_BIN_PATH" "${CONDA_PREFIX}/bin/mmseqs"
else
    echo "Creating symlink for mmseqs in conda environment..."
    ln -sf "$MMSEQS_BIN_PATH" "${CONDA_PREFIX}/bin/mmseqs"
fi

echo "MMseqs2 is available at: $(command -v mmseqs)"
"$(command -v mmseqs)" version

# --- Step 6: Download Model Weights and Databases ---
echo "=== Step 6: Downloading model weights and search databases ==="

echo "Downloading weights to ${WEIGHTS_DIR}..."
mkdir -p "$WEIGHTS_DIR"
gsutil -m rsync -r "${GCS_BUCKET}/checkpoint/" "${WEIGHTS_DIR}/"
echo "Weights download complete. Contents:"
ls -1 "$WEIGHTS_DIR"

echo "Downloading VHH search database to ${DB_DIR}..."
mkdir -p "$DB_DIR"
gsutil -m rsync -r "${GCS_BUCKET}/search_database/" "${DB_DIR}/"
echo "VHH database download complete. Contents:"
ls -1R "$DB_DIR"

COLABFOLD_TAR_REMOTE="${GCS_BUCKET}/colabfold_db.tar"
COLABFOLD_TAR_LOCAL="${DB_DIR}/colabfold_db.tar"

if [[ $SKIP_DB -eq 1 ]]; then
    echo "--no-db flag set: skipping ColabFold database download."
    echo "  (Would have downloaded: ${COLABFOLD_TAR_REMOTE})"
else
    if [[ -f "$COLABFOLD_TAR_LOCAL" ]]; then
        echo "ColabFold database already exists at ${COLABFOLD_TAR_LOCAL}. Skipping download."
    else
        echo "Downloading ColabFold database (this may take a while)..."
        gsutil -m cp "$COLABFOLD_TAR_REMOTE" "$COLABFOLD_TAR_LOCAL"
        echo "Extracting ColabFold database..."
        tar -xvf "$COLABFOLD_TAR_LOCAL" -C "$DB_DIR"
        echo "ColabFold database download and extraction complete."
    fi
fi

# --- Finalization ---
# Need to run protenix <command> a first time to build the fast_layer_norm_cuda_v2 module 
# If this works, then everything passes the test
protenix pred --help
echo ""
echo "Protenix setup complete!"
echo ""
echo "To activate the environment, run:"
echo "  conda activate $ENV_NAME"
echo ""
echo "Paths:"
echo "  - Model Weights:   $WEIGHTS_DIR"
if [[ $SKIP_DB -eq 1 ]]; then
    echo "  - Databases:       (skipped, re-run without --no-db to fetch)"
else
    echo "  - Databases:       $DB_DIR"
fi
echo "  - MMseqs2 Binary:  $MMSEQS_BIN_PATH"
echo "  - CUDA_HOME:       $CUDA_HOME"
echo ""
echo "To verify the installation, you can run:"
echo "  protenix pred --help"
echo ""
protenix pred --help