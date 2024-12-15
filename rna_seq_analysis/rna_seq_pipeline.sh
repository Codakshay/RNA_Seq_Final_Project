#!/bin/bash
#SBATCH --job-name=rna_seq_analysis
#SBATCH --mem=128G
#SBATCH --time=6:00:00
#SBATCH --ntasks=1                # One task per job
#SBATCH --cpus-per-task=16        # Number of threads per job
#SBATCH --array=1-36		  # Dynamically set array size
#SBATCH --output=logs/%A_%a.out   # Log output for each array job
#SBATCH --error=logs/%A_%a.err    # Error output for each array job

echo "Import modules and configure virtual environment"
module load sra-toolkit/3.0.9
module load minimap2/2.28
module load python/3.10
echo "generating virtual environment"
ENVDIR=/tmp/$RANDOM
virtualenv --no-download $ENVDIR
source $ENVDIR/bin/activate
pip install --no-index --upgrade pip
pip install --no-index HTSeq
pip freeze --local > requirements.txt
deactivate
rm -rf $ENVDIR

echo "Ensure directories exist"
mkdir -p fastq_files alignment logs transcripts results

echo "Fetch SRR numbers if PRJNA is provided"

echo "Get the SRR ID for the current job"
SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SRR.numbers)

echo "Downloading and converting SRA to FASTQ"                               

echo "Completed SRA -> FASTQ conversion"

echo "Download reference genome and GTF file from NCBI Homo Sapiens database"
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    echo "Downloading reference genome"
    curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    echo "Downloading GTF file from NCBI database"
    curl -O https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
    echo "Successfully downloaded reference genome and GTF file"

    # Prepare reference genome for alignment
    echo "Indexing reference genome"
    minimap2 -d ref.mmi GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    echo "Completed indexing reference genome"
fi

echo "Ensuring that the reference genome is indexed..."
while [ ! -f ref.mmi ]; do
    sleep 10
done

echo "Check if aligned SAM file already exists"
if [ ! -f "alignment/${SRR}.sam" ]; then
    echo "Assembling RNA-seq data with SRR ID: $SRR"
    minimap2 -ax sr -t 16 ref.mmi "fastq_files/${SRR}.fastq" > "alignment/${SRR}.sam"
    echo "Completed RNA-seq assembly for SRR ID: $SRR"
else
    echo "Minimap2 alignment for ${SRR}.fastq already exists. Skipping alignment."
fi

echo "Activating virtual environment"

srun --ntasks $SLURM_NNODES --tasks-per-node=1 bash << EOF

module load python/3.12

virtualenv --no-download $SLURM_TMPDIR/test_env
source $SLURM_TMPDIR/env/bin/activate

pip install --no-index --upgrade pip
pip install --no-index -r requirements.txt

EOF



echo "Generate ht-seq raw read counts"
if [ ! -f "transcripts/${SRR}.csv" ]; then
    echo "Generating ht-seq raw read counts for SRR ID: $SRR"
    htseq-count "alignment/${SRR}.sam" GCF_000001405.40_GRCh38.p14_genomic.gtf.gz -f sam -r pos -i gene_id -t exon -n 16 -o "transcripts/${SRR}.csv"
    echo "Generated ht-seq raw read counts for SRR ID: $SRR"
else
    echo "ht-seq raw read counts for ${SRR}.sam already exist. Skipping read count generation."
fi

