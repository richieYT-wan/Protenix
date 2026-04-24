conda create -c conda-forge -c pytorch -c bioconda -n protenix python=3.12 pip torch -y
conda activate protenix
pip install --upgrade protenix --index-url https://pypi.org/simple
CUDA_VERSION="$(python -c "import torch; print(torch.version.cuda)")"
conda install -c nvidia cuda-toolkit=$CUDA_VERSION
export CUDA_HOME=$CONDA_PREFIX

# Symlink headers
ln -sf $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/

# Symlink libraries too (needed for linking step)
ln -sf $CONDA_PREFIX/targets/x86_64-linux/lib/* $CONDA_PREFIX/lib/

protenix pred -i examples/input.json -o ./output -n protenix_base_20250630_v1.0.0
