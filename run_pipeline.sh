#!/usr/bin/env bash
# Usage: ./run_pipeline.sh PRJNAxxxxxx 40:30:00

PROJECT_ID=$1 # BioProject ID
TIME_LIMIT=$2 # SLURM time limit

# Create a temporary Slurm job script
cat <<EOF > run_job.sbatch
#!/usr/bin/env bash
#SBATCH --job-name=rna_seq_analysis
#SBATCH --output=logs/rna_seq_analysis.out
#SBATCH --error=logs/rna_seq_analysis.err
#SBATCH --time=${TIME_LIMIT}
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G

echo "Import modules and configure virtual environment"
module load sra-toolkit/3.0.9
module load hisat2/2.2.1
module load python/3.10

echo "Creating and activating virtual environment"
ENVDIR=\$(mktemp -d)
virtualenv --no-download \$ENVDIR
source \$ENVDIR/bin/activate

# Upgrade pip and install dependencies
pip install --no-index --upgrade pip
pip install --no-index numpy
pip install msgpack
pip install packaging
pip install blosc2
pip install tables
pip install HTSeq
pip install biomart
pip install --no-index loguru
pip install --no-index pandas
pip install --no-index biopython

echo "project: ${PROJECT_ID}" > config.yaml

snakemake make_directories --cores 32
snakemake download_fastq --cores 32

snakemake --cores 32

EOF

# Submit the job to Slurm
sbatch run_job.sbatch