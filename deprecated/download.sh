#!/bin/bash
#SBATCH --mem-per-cpu=1G
#SBATCH --time=05:30:00
#SBATCH --ntasks-per-node 1
#SBATCH --cpus-per-task=16

module load sra-toolkit/3.0.9
xargs < SRR.numbers -n 1 fasterq-dump
