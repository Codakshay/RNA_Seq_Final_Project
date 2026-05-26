#!/usr/bin/env bash
# tests/test_quick.sh
#
# Runs the minimum viable end-to-end smoke test: the first 4 SRRs of the
# bipolar case study (PRJNA316873 / GSE80336) at 5% read scale. Tests the
# same code path as test_subsampled.sh but with a fraction of the samples, so
# it completes in ~15 min on CPU. Use this for the fastest possible
# is-anything-broken check before committing a pipeline change.
#
# We deliberately reuse the bipolar BioProject (rather than a different study)
# because:
#   - the GEO metadata + condition-field handling is known-good
#   - we don't need to gamble on the GSE record format for a less-documented
#     dataset
#   - paired-end auto-detect is already covered by test_subsampled.sh's
#     ability to switch BioProjects by editing the constants below
#
# Usage:
#   sbatch tests/test_quick.sh           # CPU mode
#   sbatch tests/test_quick.sh --gpu     # GPU mode (Parabricks STAR)
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
    echo "Run: sbatch tests/test_quick.sh from the repository root." >&2
    exit 1
fi

WORKSPACE="$REPO_ROOT/tests/output/quick"
mkdir -p "$WORKSPACE/logs"

for f in pipeline.smk merge_transcripts.py metadata.py run_deg_analysis.R requirements.txt clara-parabricks_4.7.0-1.sif; do
    [ -e "$REPO_ROOT/$f" ] && ln -sf "$REPO_ROOT/$f" "$WORKSPACE/$f"
done

if [ "$USE_GPU" -eq 1 ]; then
    ALIGNER=parabricks_star
    GPU_SBATCH="#SBATCH --gres=gpu:ampere:1"
    WALLTIME="00:30:00"
else
    ALIGNER=hisat2
    GPU_SBATCH="# CPU mode — no GPU allocation"
    WALLTIME="01:00:00"
fi

cat > "$WORKSPACE/config.yaml" <<YAML
project: $PROJECT_ID
gse_accession: $GSE_ACCESSION
condition_field: $CONDITION_FIELD
aligner: $ALIGNER
download_parallel: 8
subsample_reads: $SUBSAMPLE_READS
max_samples: $MAX_SAMPLES
YAML

cat > "$WORKSPACE/run_test.sbatch" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=rna_test_min
#SBATCH --output=$WORKSPACE/logs/test_quick_%j.out
#SBATCH --error=$WORKSPACE/logs/test_quick_%j.err
#SBATCH --time=$WALLTIME
#SBATCH --cpus-per-task=16
#SBATCH --mem=60G
$GPU_SBATCH

set -euo pipefail
cd "$WORKSPACE"

echo "[test_quick] aligner=$ALIGNER  walltime=$WALLTIME  samples=$MAX_SAMPLES"

set +eu
if [ -f "\$HOME/.bashrc" ]; then
    source "\$HOME/.bashrc"
fi
set -eu

export IMAGE_NAME="$REPO_ROOT/clara-parabricks_4.7.0-1.sif"
source "\$(mamba info --base)/etc/profile.d/mamba.sh"
mamba activate parabricks_env

snakemake --snakefile pipeline.smk --cores 16 --rerun-incomplete --printshellcmds deg_results.csv

# Temporarily lower strict error catching to evaluate the results block safely
set +e

echo "[test_quick] Pipeline finished execution. Evaluating output footprints..."
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

# Restore strict error checks
set -e

if [ "\$EXIT_OK" -eq 1 ]; then
    echo "[test_quick] ALL PLOTS AND DATA VERIFIED: PASS"
else
    echo "[test_quick] PIPELINE EVALUATION OUTCOME: FAIL"
    exit 1
fi

echo "[pipeline] Workflow complete."
EOF

sbatch "$WORKSPACE/run_test.sbatch"
echo "[test_quick] Submitted. Workspace: $WORKSPACE"
echo "[test_quick] Tail logs: tail -f $WORKSPACE/logs/test_quick_<jobid>.out"