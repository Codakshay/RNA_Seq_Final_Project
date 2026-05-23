#!/usr/bin/env bash
#SBATCH --job-name=rna_seq_analysis
#SBATCH --output=logs/rna_seq_analysis_%j.out
#SBATCH --error=logs/rna_seq_analysis_%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G

set -euo pipefail

echo "[pipeline] Loading modules..."
module load sra-toolkit/3.0.9
module load hisat2/2.2.1
module load python/3.10
module load r/4.2

echo "[pipeline] Setting up Python virtual environment..."
ENVDIR=$(mktemp -d)
virtualenv --no-download "$ENVDIR"
source "$ENVDIR/bin/activate"

pip install --no-index --upgrade pip
pip install -r requirements.txt

echo "[pipeline] Starting Snakemake workflow..."
snakemake     --snakefile pipeline.smk     --cores 32     --rerun-incomplete     --printshellcmds

echo "[pipeline] Workflow complete. Cleaning up virtual environment..."
deactivate
rm -rf "$ENVDIR"
