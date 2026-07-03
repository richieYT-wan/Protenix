set -euxo pipefail
conda env create -f protenix_env.yaml -y

# Manually use pip to install various things because the environment didn't work before
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate protenix_test
cd "$(dirname ${0})"
pip install -e .
# TO BE CHANGED TO GIVE A PROPER PATH
pip install ~/ab-develop/projects/uc_denovo_vhh/foldlab
CUDA_VERSION="$(python -c "import torch; print(torch.version.cuda)")"
pip install cuda-toolkit==$CUDA_VERSION
protenix pred --help

echo "Downloading weights to ~/checkpoint"
mkdir ~/checkpoint
gsutil cp -r gs://em52-ab-develop-analytics-prod-f684/data/denovo-design/snakemake/data/01_raw/protenix-inputs/checkpoint/* ~/checkpoint/
ls ~/checkpoint

# Set up MMSEQS and databases if needed here
#wget https://github.com/soedinglab/MMseqs2/archive/refs/tags/16-747c6.tar.gz
#tar xzf 16-747c6.tar.gz
#cd MMseqs2-16-747c6/
#mkdir build && cd build
#cmake -DCMAKE_BUILD_TYPE=RELEASE -DCMAKE_INSTALL_PREFIX=. ..
#make -j8
#make install
#git clone https://github.com/sokrypton/ColabFold.git
#cd ColabFold
## Configure database:
#MMSEQS_NO_INDEX=1 ./setup_databases.sh <path/to/colabfold_db>

#export CUDA_HOME=$CONDA_PREFIX
## Symlink headers
#ln -sf $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/
#
## Symlink libraries too (needed for linking step)
#ln -sf $CONDA_PREFIX/targets/x86_64-linux/lib/* $CONDA_PREFIX/lib/
#
#protenix pred -i examples/input.json -o ./output -n protenix_base_20250630_v1.0.0
