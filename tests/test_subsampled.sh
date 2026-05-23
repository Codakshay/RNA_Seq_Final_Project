#!/usr/bin/env bash
# tests/test_subsampled.sh
#
# Runs the published bipolar case study (PRJNA316873 / GSE80336) end-to-end at
# 1% read scale: all 36 samples, but each FASTQ is capped at 100k reads via
# `fasterq-dump -X 100000`. Whole pipeline finishes in ~30 min on CPU, ~10 min
# on GPU. Use this to validate that a pipeline change still produces the right
# shape (counts.csv columns, deg_results.csv, all 7 plot PNGs) without paying
# the full 12-hour case-study cost.
#
# Usage:
#   sbatch tests/test_subsampled.sh           # CPU mode (HISAT2)
#   sbatch tests/test_subsampled.sh --gpu     # GPU mode (Parabricks STAR)

set -euo pipefail

USE_GPU=0
[[ "${1:-}" == "--gpu" ]] && USE_GPU=1

PROJECT_ID=PRJNA316873
GSE_ACCESSION=GSE80336
CONDITION_FIELD=title
SUBSAMPLE_READS=100000

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
WORKSPACE="$REPO_ROOT/tests/output/subsampled"
mkdir -p "$WORKSPACE/logs"

# Symlink pipeline files into the workspace (so we never clobber a real run's
# fastq_files/, alignment/, data/, etc. in the repo root)
for f in pipeline.smk merge_transcripts.py metadata.py run_deg_analysis.R requirements.txt; do
    ln -sf "$REPO_ROOT/$f" "$WORKSPACE/$f"
done

if [ "$USE_GPU" -eq 1 ]; then
    ALIGNER=parabricks_star
    GPU_SBATCH="#SBATCH --gres=gpu:v100:1"
    EXTRA_MOD="module load parabricks cuda star/2.7.10b subread/2.0.6 samtools/1.17"
    WALLTIME="00:45:00"
else
    ALIGNER=hisat2
    GPU_SBATCH="# CPU mode — no GPU allocation"
    EXTRA_MOD="module load hisat2/2.2.1 subread/2.0.6 samtools/1.17"
    WALLTIME="01:30:00"
fi

cat > "$WORKSPACE/config.yaml" <<YAML
project: $PROJECT_ID
gse_accession: $GSE_ACCESSION
condition_field: $CONDITION_FIELD
aligner: $ALIGNER
download_parallel: 4
subsample_reads: $SUBSAMPLE_READS
YAML

cat > "$WORKSPACE/run_test.sbatch" <<EOF
#!/usr/bin/env bash
#SBATCH --job-name=rna_test_sub
#SBATCH --output=$WORKSPACE/logs/test_subsampled_%j.out
#SBATCH --error=$WORKSPACE/logs/test_subsampled_%j.err
#SBATCH --time=$WALLTIME
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
$GPU_SBATCH

set -euo pipefail
cd "$WORKSPACE"

echo "[test_subsampled] aligner=$ALIGNER  walltime=$WALLTIME"
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
    echo "[test_subsampled] PASS"
else
    echo "[test_subsampled] FAIL"
    exit 1
fi

deactivate
rm -rf "\$ENVDIR"
EOF

sbatch "$WORKSPACE/run_test.sbatch"
echo "[test_subsampled] Submitted. Workspace: $WORKSPACE"
echo "[test_subsampled] Tail logs: tail -f $WORKSPACE/logs/test_subsampled_<jobid>.out"
