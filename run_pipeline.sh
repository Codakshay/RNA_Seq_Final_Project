#!/usr/bin/env bash
# Run the RNA-Seq DGE pipeline on a SLURM cluster.
#
# Usage:
#   sbatch run_pipeline.sh <PRJNA_ID> <HH:MM:SS> <GSE_ACCESSION> [CONDITION_FIELD] [--gpu]
#
# Arguments:
#   PRJNA_ID        NCBI BioProject accession (e.g. PRJNA316873)
#   HH:MM:SS        SLURM time limit for the job
#   GSE_ACCESSION   GEO Series accession for sample metadata (e.g. GSE80336)
#   CONDITION_FIELD Column in GEO phenoData used to derive condition labels
#                   (default: title)
#   --gpu           Optional flag. When set, the job is submitted to a GPU node
#                   and pipeline.smk uses NVIDIA Parabricks STAR for alignment.
#
# Examples:
#   sbatch run_pipeline.sh PRJNA316873 24:00:00 GSE80336 title          # CPU (HISAT2)
#   sbatch run_pipeline.sh PRJNA316873 04:00:00 GSE80336 title --gpu    # GPU (Parabricks STAR)

PROJECT_ID="${1:?Error: PRJNA_ID is required (e.g. PRJNA316873)}"
TIME_LIMIT="${2:?Error: time limit is required (e.g. 24:00:00)}"
GSE_ACCESSION="${3:?Error: GSE_ACCESSION is required (e.g. GSE80336)}"

# CONDITION_FIELD is optional, --gpu is optional; either may appear in slot 4 or 5
CONDITION_FIELD="title"
USE_GPU=0
for arg in "${@:4}"; do
    case "$arg" in
        --gpu) USE_GPU=1 ;;
        *)     CONDITION_FIELD="$arg" ;;
    esac
done

# Create logs/ on the host before sbatch so SLURM can write its own log files
mkdir -p logs

# Compose the aligner-specific bits
if [ "$USE_GPU" -eq 1 ]; then
    ALIGNER="parabricks_star"
    SBATCH_GPU_LINE="#SBATCH --gres=gpu:ampere:1"
    echo "[run_pipeline.sh] GPU mode: aligner=parabricks_star, --gres=gpu:ampere:1"
else
    ALIGNER="hisat2"
    SBATCH_GPU_LINE="# (CPU mode — no GPU allocation)"
    echo "[run_pipeline.sh] CPU mode: aligner=hisat2"
fi

# Generate config.yaml from the supplied arguments
cat > config.yaml <<YAML
project: ${PROJECT_ID}
gse_accession: ${GSE_ACCESSION}
condition_field: ${CONDITION_FIELD}
aligner: ${ALIGNER}
download_parallel: 4
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
${SBATCH_GPU_LINE}

set -euo pipefail

eval "$(mamba shell hook --shell bash)"
mamba activate rna_seq

echo "[pipeline] Starting Snakemake workflow (aligner=${ALIGNER})..."

snakemake \
    --snakefile pipeline.smk \
    --cores "${SLURM_CPUS_PER_TASK}" \
    --rerun-incomplete \
    --printshellcmds

echo "[pipeline] Workflow complete."
EOF

sbatch run_job.sbatch
echo "Job submitted. Monitor with: sq"
echo "Logs: logs/rna_seq_analysis_<jobid>.out / .err"
