#!/bin/bash
#SBATCH --job-name=minimap2_assembly
#SBATCH --mem=100G
#SBATCH --time=02:30:00
#SBATCH --ntasks=1                # One task per job
#SBATCH --cpus-per-task=32        # Number of threads per job

module load hisat2/2.2.1

# perform genome assembly for the given fastq files

# Step 2: Get the SRR ID for the current job
hisat2-build -p 32 GCF_000001405.40_GRCh38.p14_genomic.fna ref_genome


