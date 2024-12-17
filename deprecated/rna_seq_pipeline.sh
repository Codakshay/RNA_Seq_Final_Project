#!/bin/bash
#SBATCH --job-name=rna_seq_analysis
#SBATCH --mem=120G
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --array=1-36
#SBATCH --output=logs/%A_%a.out
#SBATCH --error=logs/%A_%a.err

echo "Loading required modules..."
module load sra-toolkit/3.0.9
module load hisat2/2.2.1
module load python/3.10

# Ensure required directories exist
echo "Creating required directories..."
mkdir -p fastq_files alignment logs transcripts results data ref_genome

# Get the SRR ID for the current job
if [ ! -f "SRR.numbers" ]; then
    echo "Error: SRR.numbers file not found!"
    exit 1
fi

SRR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" SRR.numbers)
if [ -z "$SRR" ]; then
    echo "Error: No SRR ID found for task ID ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

echo "Processing SRR ID: $SRR"

# Download and convert SRA to FASTQ
if [ ! -f "fastq_files/${SRR}.fastq" ]; then
    echo "Downloading and converting SRA to FASTQ for $SRR..."
    fasterq-dump -O fastq_files --split-files "$SRR"
    if [ $? -ne 0 ]; then
        echo "Error: Failed to download or convert $SRR to FASTQ."
        exit 1
    fi
else
    echo "FASTQ file for $SRR already exists. Skipping download."
fi

# Download reference genome and GTF file (only once for task 1)
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    if [ ! -f "data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz" ]; then
        echo "Downloading reference genome..."
        curl -o data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
    fi

    if [ ! -f "data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz" ]; then
        echo "Downloading GTF file..."
        curl -o data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
    fi

    # Prepare reference genome for alignment
    if [ ! -f "ref_genome/ref_gen.1.ht2" ]; then
        echo "Indexing reference genome..."
        gunzip -c data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz > data/GCF_000001405.40_GRCh38.p14_genomic.fna
        hisat2-build -p 32 data/GCF_000001405.40_GRCh38.p14_genomic.fna ref_genome/ref_gen
    fi
fi

# Wait for the reference genome to be ready (race condition prevention)
while [ ! -f "ref_genome/ref_gen.1.ht2" ]; do
    echo "Waiting for reference genome indexing to complete..."
    sleep 60
done

# Align reads with HISAT2
if [ ! -f "alignment/${SRR}.sam" ]; then
    echo "Aligning reads for $SRR..."
    hisat2 -p 32 -x ref_genome/ref_gen -U fastq_files/${SRR}_1.fastq -S alignment/${SRR}.sam
    if [ $? -ne 0 ]; then
        echo "Error: Alignment failed for $SRR."
        exit 1
    fi
else
    echo "Alignment file for $SRR already exists. Skipping alignment."
fi

# Generate read counts with HTSeq
if [ ! -f "transcripts/${SRR}.csv" ]; then
    echo "Generating read counts for $SRR..."
    htseq-count -f sam -r pos -i gene_id -t exon -n 32 alignment/${SRR}.sam data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz > transcripts/${SRR}.csv
    if [ $? -ne 0 ]; then
        echo "Error: HTSeq count failed for $SRR."
        exit 1
    fi
else
    echo "Read counts for $SRR already exist. Skipping HTSeq count generation."
fi

# Merge transcripts and run DEG analysis (only on the last task)
if [ "${SLURM_ARRAY_TASK_ID}" -eq "${SLURM_ARRAY_TASK_MAX}" ]; then
    echo "Merging transcripts and running differential expression analysis..."
    python merge_transcripts.py
    Rscript run_deg_analysis.r
fi

echo "Task completed for $SRR."
