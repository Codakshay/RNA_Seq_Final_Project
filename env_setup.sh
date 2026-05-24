#!/usr/bin/env bash
set -euo pipefail

# 1. Initialize Mamba shell environment hooks
eval "$(mamba shell hook --shell bash)"

# 2. Create the base environment with robust bioconda & python data packages
echo "[setup] Creating core mamba environment..."
mamba create --name rna_seq -c conda-forge -c bioconda \
    python=3.10 pip samtools=1.19 subread=2.0.6 r-base=4.3 -y

# Activate the environment context immediately
mamba activate rna_seq

# 3. Safely install python pipeline managers using pip
echo "[setup] Installing python workflow packages..."
pip install --upgrade pip
pip install snakemake==7.32.4 matplotlib pandas

# 4. Configure Apptainer staging environments on high-capacity scratch storage
export APPTAINER_CACHEDIR="/home/dzk5572/scratch/apptainer_cache"
export APPTAINER_TMPDIR="/home/dzk5572/scratch/apptainer_tmp"
mkdir -p "$APPTAINER_CACHEDIR" "$APPTAINER_TMPDIR"

# 5. Pull and build the Parabricks SIF binary
echo "[setup] Building Clara Parabricks Apptainer SIF (this takes a few minutes)..."
IMAGE_NAME="clara-parabricks_4.7.0-1.sif"
apptainer pull --no-xattrs "$IMAGE_NAME" docker://nvcr.io/nvidia/clara/clara-parabricks:4.7.0-1

# 6. Verify image build and hardware link
echo "[setup] Verifying GPU/Container interface topology..."
apptainer exec --nv "$IMAGE_NAME" pbrun version

# 7. Create a permanent run alias inside this mamba environment
# This writes a tiny activation script so whenever you run 'mamba activate rna_seq', 
# 'pbrun' automatically points to this local Apptainer SIF container with GPU mapping.
ACTIVATE_DIR="${CONDA_PREFIX}/etc/conda/activate.d"
mkdir -p "$ACTIVATE_DIR"
cat > "${ACTIVATE_DIR}/parabricks_alias.sh" <<EOF
#!/bin/sh
alias pbrun="apptainer exec --nv --bind \$(pwd):/data $(pwd)/$IMAGE_NAME pbrun"
echo "[rna_seq] Parabricks alias 'pbrun' registered with GPU execution hooks."
EOF

echo "[setup] Environment complete. Run 'mamba activate rna_seq' to begin tracking pipelines."