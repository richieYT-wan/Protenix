set -euxo pipefail
conda env create -f protenix_env.yaml -y

# Manually use pip to install various things because the environment didn't work before
sleep 1
conda activate protenix
CUDA_VERSION="$(python -c "import torch; print(torch.version.cuda)")"
pip install cuda-toolkit==$CUDA_VERSION
export CUDA_HOME=$CONDA_PREFIX

# Symlink headers
ln -sf $CONDA_PREFIX/targets/x86_64-linux/include/* $CONDA_PREFIX/include/

# Symlink libraries too (needed for linking step)
ln -sf $CONDA_PREFIX/targets/x86_64-linux/lib/* $CONDA_PREFIX/lib/

protenix pred -i examples/input.json -o ./output -n protenix_base_20250630_v1.0.0
