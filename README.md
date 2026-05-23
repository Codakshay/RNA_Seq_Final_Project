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

The pipeline supports RNA-Seq differential gene expression analysis for *Homo sapiens* **single-end** sequencing data. All datasets (raw reads and reference annotations) are downloaded automatically.

### Arguments

```
sbatch run_pipeline.sh <PRJNA_ID> <HH:MM:SS> <GSE_ACCESSION> [CONDITION_FIELD]
```

| Argument | Description | Example |
|---|---|---|
| `PRJNA_ID` | NCBI BioProject accession | `PRJNA316873` |
| `HH:MM:SS` | SLURM wall-clock time limit | `24:00:00` |
| `GSE_ACCESSION` | GEO Series accession for sample metadata | `GSE80336` |
| `CONDITION_FIELD` | Column in GEO phenoData used to derive group labels (default: `title`) | `title` |

### Example — reproducing the bipolar disorder case study

```sh
sbatch run_pipeline.sh PRJNA316873 24:00:00 GSE80336 title
```

The script generates `config.yaml` automatically from the arguments above, so no manual file editing is required. To inspect or adjust the configuration beforehand, copy the template:

```sh
cp config.yaml.example config.yaml
# edit config.yaml as needed, then:
sbatch run_pipeline.sh PRJNA316873 24:00:00 GSE80336 title
```

While the job runs, monitor its status with `sq`. Logs are written to `logs/rna_seq_analysis_<jobid>.out` and `.err`.

### Output files

Upon successful completion the following files are produced:

| Path | Contents |
|---|---|
| `fastq_files/SRR.numbers` | List of SRR run accessions downloaded |
| `transcripts/<SRR>.csv` | Per-sample raw read counts from htseq-count |
| `data/counts.csv` | Merged annotated count matrix (Ensembl ID + gene info + all samples) |
| `deg_results.csv` | DESeq2 results ordered by adjusted p-value |
| `data/plots/PCAPlot.png` | Principal Component Analysis |
| `data/plots/MAPlot.png` | MA plot (raw LFC) |
| `data/plots/resMAPlot.png` | MA plot (LFC-shrunken) |
| `data/plots/VolcanoPlot.png` | Volcano plot of DEGs |
| `data/plots/DispersionPlot.png` | DESeq2 dispersion plot |
| `data/plots/HeatmapPairwisePlot.png` | Sample-to-sample distance heatmap |
| `data/plots/HeatmapDEGPlot.png` | Heatmap of top 2 000 DEGs |

> **Note**: if the specified time limit is insufficient the SLURM job will be terminated automatically. Re-submit with a longer allocation. Snakemake's `--rerun-incomplete` flag (used by default) ensures it picks up where it left off.
