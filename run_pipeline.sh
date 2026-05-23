#!/usr/bin/env bash
# Run the RNA-Seq DGE pipeline on a SLURM cluster.
#
# Usage:
#   sbatch run_pipeline.sh <PRJNA_ID> <HH:MM:SS> <GSE_ACCESSION> [CONDITION_FIELD]
#
# Arguments:
#   PRJNA_ID        NCBI BioProject accession (e.g. PRJNA316873)
#   HH:MM:SS        SLURM time limit for the job
#   GSE_ACCESSION   GEO Series accession for sample metadata (e.g. GSE80336)
#   CONDITION_FIELD Column in GEO phenoData used to derive condition labels
#                   (default: title)
#
# Example:
#   sbatch run_pipeline.sh PRJNA316873 24:00:00 GSE80336 title

PROJECT_ID="${1:?Error: PRJNA_ID is required (e.g. PRJNA316873)}"
TIME_LIMIT="${2:?Error: time limit is required (e.g. 24:00:00)}"
GSE_ACCESSION="${3:?Error: GSE_ACCESSION is required (e.g. GSE80336)}"
CONDITION_FIELD="${4:-title}"

# Create logs/ on the host before sbatch so SLURM can write its own log files
mkdir -p logs

# Generate config.yaml from the supplied arguments
cat > config.yaml <<YAML
project: ${PROJECT_ID}
gse_accession: ${GSE_ACCESSION}
condition_field: ${CONDITION_FIELD}
YAML

# Create and submit the SLURM job script
cat > run_job.sbatch <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=rna_seq_analysis
#SBATCH --output=logs/rna_seq_analysis_%j.out
#SBATCH --error=logs/rna_seq_analysis_%j.err
#SBATCH --time=${TIME_LIMIT}
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G

set -euo pipefail

echo "[pipeline] Loading modules..."
module load sra-toolkit/3.0.9
module load hisat2/2.2.1
module load python/3.10
module load r/4.2

echo "[pipeline] Setting up Python virtual environment..."
ENVDIR=\$(mktemp -d)
virtualenv --no-download "\$ENVDIR"
source "\$ENVDIR/bin/activate"

pip install --no-index --upgrade pip
pip install -r requirements.txt

echo "[pipeline] Starting Snakemake workflow..."
snakemake \
    --snakefile pipeline.smk \
    --cores 32 \
    --rerun-incomplete \
    --printshellcmds

echo "[pipeline] Workflow complete. Cleaning up virtual environment..."
deactivate
rm -rf "\$ENVDIR"
EOF

sbatch run_job.sbatch
echo "Job submitted. Monitor with: sq"
echo "Logs: logs/rna_seq_analysis_<jobid>.out / .err"
