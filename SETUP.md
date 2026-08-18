# Protenix Fork Setup Guide

This document provides instructions for setting up the development and execution environment for this forked Protenix repository.

The unified `install.sh` script automates the setup process, ensuring a reproducible environment.

## Prerequisites

Before running the installation script, please ensure your system meets the following requirements:

1.  **Operating System**: A Debian-based Linux distribution (e.g., Ubuntu 20.04+). The script is designed for environments like a Google Vertex AI Workbench instance or a Google Cloud Workstation.
2.  **GPU**: A CUDA-capable NVIDIA GPU with the appropriate drivers installed.
3.  **Conda/Mamba**: `conda` or `mamba` must be installed and available in your shell. `mamba` is recommended for faster environment setup.
4.  **Google Cloud SDK**: The `gcloud` CLI must be installed and authenticated. You need to have `gsutil` available, which is included with the SDK.
5.  **Build Tools**: Standard build-essential tools (`make`, `g++`, `cmake`, etc.) are required for compiling `MMseqs2`. These are typically pre-installed on Vertex AI instances.

## Installation

The installation process is handled by a single script.

1.  **Clone the repository**:
    ```bash
    git clone https://github.com/richieYT-wan/Protenix.git
    cd Protenix
    ```

2.  **Run the installer**:
    ```bash
    bash install.sh
    ```

The script will perform the following steps:
- Create a Conda environment named `protenix`.
- Install the Protenix package mode.
- Install the correct CUDA toolkit version matching your PyTorch installation.
- Verify the `HMMER` installation from Conda.
- Download and build a specific version of `MMseqs2` from source.
- Download the required model weights and search databases from Google Cloud Storage.

## Post-Installation

### Activating the Environment

To use the Protenix tools and run the code, you must first activate the Conda environment:

```bash
conda activate protenix
```

### Installed Components

The script places files in the following locations in your home directory:

-   **Model Weights**: `~/checkpoint/`
-   **Search Databases**: `~/protenix_dbs/`
-   **External Binaries**: `~/protenix_bin/` (contains the `mmseqs2` installation)

These paths can be configured by editing the variables at the top of the `install.sh` script before running it.

### Verifying the Installation

To confirm that the installation was successful and the command-line interface is working, activate the environment and run the following command:

```bash
protenix pred --help
```

This should display the help menu for the prediction command, indicating that the `protenix` entry point is correctly installed in your PATH.

You can also verify that `mmseqs2` is correctly linked and available:
```bash
mmseqs version
```
This should output the version `16-747c6`.