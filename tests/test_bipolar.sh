#!/usr/bin/env bash
#
# bipolar case study
# We deliberately use the bipolar BioProject
# Usage:
#   sbatch tests/test_minimal.sh           # CPU mode
#   sbatch tests/test_minimal.sh --gpu     # GPU mode (Parabricks STAR)
#
# To repurpose for a different BioProject, edit the PROJECT_ID/GSE_ACCESSION/
# CONDITION_FIELD constants below.

set -euo pipefail

USE_GPU=0
[[ "${1:-}" == "--gpu" ]] && USE_GPU=1

PROJECT_ID=PRJNA318642
GSE_ACCESSION=GSE80336
CONDITION_FIELD=title
MAX_SAMPLES=0
SUBSAMPLE_READS=0

if [[ -n "${SLURM_JOB_ID:-}" ]]; then
    REPO_ROOT="${SLURM_SUBMIT_DIR}"
else
    REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
fi

if [[ ! -f "$REPO_ROOT/pipeline.smk" ]]; then
    echo "ERROR: REPO_ROOT=$REPO_ROOT does not contain pipeline.smk" >&2
    echo "Run: sbatch tests/test_minimal.sh from the repository root." >&2
    exit 1
fi

WORKSPACE="$REPO_ROOT/tests/output/bipolar"
mkdir -p "$WORKSPACE/logs"

for f in pipeline.smk merge_transcripts.py metadata.py run_deg_analysis.R requirements.txt; do
    ln -sf "$REPO_ROOT/$f" "$WORKSPACE/$f"
done

if [ "$USE_GPU" -eq 1 ]; then
    ALIGNER=parabricks_star
    GPU_SBATCH="#SBATCH --gres=gpu:ampere:1"
    WALLTIME="24:30:00"
else
    ALIGNER=hisat2
    GPU_SBATCH="# CPU mode — no GPU allocation"
    WALLTIME="30:00:00"
fi

cat > "$WORKSPACE/config.yaml" <<YAML
project: $PROJECT_ID
gse_accession: $GSE_ACCESSION
condition_field: $CONDITION_FIELD
aligner: $ALIGNER
download_parallel: 4
subsample_reads: $SUBSAMPLE_READS
max_samples: $MAX_SAMPLES
YAML

cat > "$WORKSPACE/run_test.sbatch" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=rna_test_min
#SBATCH --output=$WORKSPACE/logs/test_minimal_%j.out
#SBATCH --error=$WORKSPACE/logs/test_minimal_%j.err
#SBATCH --time=$WALLTIME
#SBATCH --cpus-per-task=32
#SBATCH --mem=80G
$GPU_SBATCH
set -euo pipefail
cd "$WORKSPACE"
echo "[test_minimal] aligner=$ALIGNER  walltime=$WALLTIME  samples=$MAX_SAMPLES"
set +eu
if [ -f "\$HOME/.bashrc" ]; then
    source "\$HOME/.bashrc"
fi
set -eu
source "\$(mamba info --base)/etc/profile.d/mamba.sh"
mamba activate parabricks_env
snakemake --snakefile pipeline.smk --cores 32 --rerun-incomplete --printshellcmds deg_results.csv
echo "[test_minimal] Pipeline finished execution. Evaluating output footprints..."
EXIT_OK=1
if [ ! -f deg_results.csv ]; then
    echo "FAIL: deg_results.csv missing"
    EXIT_OK=0
fi
for p in PCAPlot MAPlot resMAPlot VolcanoPlot DispersionPlot HeatmapPairwisePlot HeatmapDEGPlot; do
    if [ ! -f "data/plots/\${p}.png" ]; then
        echo "FAIL: data/plots/\${p}.png missing"
        EXIT_OK=0
    fi
done
set -e
if [ "\$EXIT_OK" -eq 1 ]; then
    echo "[test_minimal] ALL PLOTS AND DATA VERIFIED: PASS"
else
    echo "[test_minimal] PIPELINE EVALUATION OUTCOME: FAIL"
    exit 1
fi
echo "[pipeline] Workflow complete."
EOF

sbatch "$WORKSPACE/run_test.sbatch"
echo "[test_minimal] Submitted. Workspace: $WORKSPACE"
echo "[test_minimal] Tail logs: tail -f $WORKSPACE/logs/test_minimal_<jobid>.out"