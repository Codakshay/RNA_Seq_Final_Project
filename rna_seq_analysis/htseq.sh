#!/bin/bash
#SBATCH --job-name=rna_seq_analysis
#SBATCH --mem=100G
#SBATCH --time=03:30:00
#SBATCH --ntasks=1                # One task per job
#SBATCH --cpus-per-task=32        # Number of threads per job
#SBATCH --array=1-1              # Dynamically set array size
#SBATCH --output=logs/%A_%a.out   # Log output for each array job
#SBATCH --error=logs/%A_%a.err    # Error output for each array job

echo "Import modules and configure virtual environment"
module load sra-toolkit/3.0.9
module load hisat2/2.2.1
module load python/3.10

echo "Generating virtual environment"
ENVDIR=/tmp/$RANDOM
virtualenv --no-download $ENVDIR
source $ENVDIR/bin/activate
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
pip freeze --local > requirements.txt
deactivate

echo "activating temporary environment"
module load python/3.10
VENV=/tmp/$RANDOM
virtualenv --no-download $VENV
source $VENV/bin/activate
pip install --no-index --upgrade pip
pip install -r requirements.txt

# Check if the environment was set up correctly
if [ ! -d $VENV ]; then
    echo "Error setting up virtual environment."
    exit 1
fi

hisat2-build -p 32 "data/GCF_000001405.40_GRCh38.p14_genomic.fna" "ref_genome/ref_gen"

# Activate the virtual environment
source $VENV/bin/activate
echo "Successfully activated virtual environment"

# Generate ht-seq raw read counts
echo "Generating ht-seq raw read counts"
hisat2 -p 32 -x ref_genome/ref_gen -U "fastq_files/SRR3393492.fastq" -S "alignment/SRR3393492.sam"
htseq-count "alignment/SRR3393492.sam" "data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz" -f sam -r pos -i gene_id -t exon -n 32 > "transcripts/SRR3393492.csv"
echo "ht-seq raw read counts for SRR3393523.sam already exist. Skipping read count generation."

# Deactivate and remove the virtual environment
deactivate
echo "deactivated"
rm -rf $VENV

