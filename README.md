# RNA_Seq_Final_Project
This repository presents a comprehensive analysis of RNA-Seq data to investigate the molecular underpinnings of bipolar disorder (BD), alongside the implementation of a computational pipeline optimized for high-performance computing (HPC) environments. The pipeline performs a robust and systematic analysis of raw Illumina high-throughput RNA-Seq reads, ensuring efficient and accurate data processing. The study outlines the rationale, methodology, and implementation details of each step in the RNA-Seq data analysis workflow. RNA-Seq data from the human dorsal striatum were utilized to profile transcriptomes, comparing cohorts of bipolar disorder patients to healthy controls. Differential gene expression (DGE) analysis was conducted using the \texttt{DESeq2} library to identify genes exhibiting significant changes in expression associated with bipolar disorder status. The results were further analyzed to uncover sets of co-expressed and correlated genes, providing insights into their potential association with bipolar disorder.

## Design and Specifications

The computational pipeline, which streamlines the complete RNA-Seq differential gene expression (DGE) analysis process, is publicly available on GitHub at [https://github.com/CebolaLab/RNA-seq](https://github.com/CebolaLab/RNA-seq). Designed as a command-line application, the pipeline offers an efficient, automated, and reproducible solution for high-performance computing (HPC) environments. It integrates multiple tools and processes seamlessly using the `slurm` workload manager and runs on Linux operating systems such as AlmaLinux 9 and CentOS 7.

The pipeline is implemented as a modular framework, combining Python, R, shell scripts, and a `Snakemake` workflow. By leveraging `Snakemake`, the pipeline allows task parallelization and efficient resource management across multiple compute nodes. This design significantly reduces processing time, making it suitable for large-scale RNA-Seq datasets.

| Feature | Specification |
|---|---|
| Operating System | Linux (AlmaLinux 9, CentOS 7) |
| Workload Management | slurm |
| Cores | minimum of 32 cores |
| Memory | minimum of 50 gigabytes + total size of FASTQ files |
| External Dependencies | python (≥ 3.9.0), R (≥ 4.0.5) |
| Estimated Processing Time | 24-48 hours depending on the FASTQ data size |

**Table 1: System Requirements**

To ensure optimal performance and reproducibility, the pipeline adheres to the specifications outlined in Table 1. These requirements include the minimum number of CPU cores, memory size, and software dependencies necessary for the successful execution of the analysis. Depending on the size of the input `FASTQ` files, the estimated runtime for a complete analysis ranges between 24 and 48 hours.

By consolidating all steps of RNA-Seq data analysis into a single pipeline, this implementation provides an effective one-click solution for large-scale genomic studies in high performance computing (HPC) environments.

## Installation

We recommend installing the pipeline directly from the GitHub repository to ensure access to the most recent version, including any updates or fixes. To install the pipeline on your computing cluster, execute the following command:

```sh
git clone [https://github.com/Codakshay/RNA_Seq_Final_Project.git](https://github.com/Codakshay/RNA_Seq_Final_Project.git)
```

For all subsequent steps, ensure you navigate to the root directory of the cloned repository.

Before running the pipeline, it is essential to verify that your computing cluster supports both R and Python as available modules. Most high-performance computing (HPC) environments use module management systems (e.g., module load) to provide software packages. You should confirm that R (version ≥ 4.0.5) and Python (version ≥ 3.9.0) are accessible. If these are not pre-installed, contact your system administrator to install or enable them.

Additionally, ensure that all Python dependencies listed in the requirements.txt file located in the root directory of the repository are installed on the cluster. The required libraries can be installed locally in a virtual environment or within a Conda environment. To set up the environment and install the necessary libraries, execute the following commands:

```sh
python3 -m venv env
source env/bin/activate
pip install -r requirements.txt
```
This will create an isolated virtual environment containing all the required Python dependencies for running the pipeline. If you encounter any issues, verify that the correct Python version is being used and that the required libraries are available.

## Running the Pipeline

This subsection provides a concise manual for executing the computational pipeline on a high-performance computing (HPC) environment. Before proceeding, please verify that your computing system meets the design specifications outlined in the previous subsection.

The pipeline requires only a BioProject ID as input. At its current stage, the pipeline exclusively supports RNA-Seq differential gene expression analysis for Homo sapiens data. The necessary datasets, including raw sequencing reads and reference annotations, are automatically retrieved from the NCBI FTP server by the pipeline script.

In the NCBI Sequence Read Archive (SRA) database, BioProject IDs provide unique identifiers for retrieving all raw sequencing datasets associated with a particular project. These identifiers are typically formatted as \texttt{PRJNA} followed by a series of digits. BioProject IDs can be located either via the Gene Expression Omnibus (GEO) database or by directly searching the NCBI database. As shown in Figure~\ref{fig:ncbi-navigator}, users can find the BioProject ID by selecting “BioProject” in the search criteria. Once the appropriate BioProject ID is retrieved, users must navigate to the root directory of the cloned pipeline repository and execute the following command on their computing cluster:
```sh
sbatch run_pipeline.sh PRJNA123456 24:00:00
```
In this command, PRJNA123456 represents the BioProject ID, and 24\:00\:00 specifies the maximum expected execution time in HH\:MM\:SS format. While the total processing time depends on the number of \texttt{FASTQ} files and their sizes, typical runtimes range between 24 and 48 hours. For the case study conducted in this work, the pipeline completed execution in approximately 12 hours.

It is important to note that if the specified time limit is insufficient, the job will be automatically terminated by the SLURM workload manager. In such cases, the pipeline must be re-executed with an increased time allocation. Users are advised to monitor the execution status of their job by running the `sq` command, which provides real-time updates on job progress within the computing cluster.

While the execution process itself does not require user intervention, the status and potential errors can be monitored by inspecting the `.out` and `.err` files generated in the \texttt{RNA\_Seq\_Final\_Project/logs/} directory. These log files provide detailed information about the ongoing pipeline execution, including standard outputs and any error messages encountered during runtime.
