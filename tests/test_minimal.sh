#!/usr/bin/env bash
# tests/test_minimal.sh
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
#   sbatch tests/test_minimal.sh           # CPU mode
#   sbatch tests/test_minimal.sh --gpu     # GPU mode (Parabricks STAR)
#
# To repurpose for a different BioProject, edit the PROJECT_ID/GSE_ACCESSION/
# CONDITION_FIELD constants below.

set -euo pipefail

USE_GPU=0
[[ "${1:-}" == "--gpu" ]] && USE_GPU=1

PROJECT_ID=PRJNA316873
GSE_ACCESSION=GSE80336
CONDITION_FIELD=title
MAX_SAMPLES=4
SUBSAMPLE_READS=50000

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKSPACE="$REPO_ROOT/tests/output/minimal"
mkdir -p "$WORKSPACE/logs"

for f in pipeline.smk merge_transcripts.py metadata.py run_deg_analysis.R requirements.txt; do
    ln -sf "$REPO_ROOT/$f" "$WORKSPACE/$f"
done

if [ "$USE_GPU" -eq 1 ]; then
    ALIGNER=parabricks_star
    GPU_SBATCH="#SBATCH --gres=gpu:v100:1"
    EXTRA_MOD="module load parabricks cuda star/2.7.10b subread/2.0.6 samtools/1.17"
    WALLTIME="00:30:00"
else
    ALIGNER=hisat2
    GPU_SBATCH="# CPU mode — no GPU allocation"
    EXTRA_MOD="module load hisat2/2.2.1 subread/2.0.6 samtools/1.17"
    WALLTIME="00:45:00"
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
#SBATCH --mem=64G
$GPU_SBATCH

set -euo pipefail
cd "$WORKSPACE"

echo "[test_minimal] aligner=$ALIGNER  walltime=$WALLTIME  samples=$MAX_SAMPLES"
module load sra-toolkit/3.0.9 python/3.10 r/4.2
$EXTRA_MOD

ENVDIR=\$(mktemp -d)
virtualenv --no-download "\$ENVDIR"
source "\$ENVDIR/bin/activate"
pip install --no-index --upgrade pip
pip install -r requirements.txt

snakemake --snakefile pipeline.smk --cores 32 --rerun-incomplete --printshellcmds

# Pass/fail gate
EXIT_OK=1
[ -f deg_results.csv ] || { echo "FAIL: deg_results.csv missing"; EXIT_OK=0; }
for p in PCAPlot MAPlot resMAPlot VolcanoPlot DispersionPlot HeatmapPairwisePlot HeatmapDEGPlot; do
    [ -f "data/plots/\$p.png" ] || { echo "FAIL: data/plots/\$p.png missing"; EXIT_OK=0; }
done
if [ "\$EXIT_OK" -eq 1 ]; then
    echo "[test_minimal] PASS"
else
    echo "[test_minimal] FAIL"
    exit 1
fi

deactivate
rm -rf "\$ENVDIR"
EOF

sbatch "$WORKSPACE/run_test.sbatch"
echo "[test_minimal] Submitted. Workspace: $WORKSPACE"
echo "[test_minimal] Tail logs: tail -f $WORKSPACE/logs/test_minimal_<jobid>.out"
