configfile: "config.yaml"
PROJECT         = config["project"]
GSE_ACCESSION   = config["gse_accession"]
CONDITION_FIELD = config["condition_field"]

# Plot filenames that run_deg_analysis.R will produce
PLOT_NAMES = [
    "PCAPlot", "MAPlot", "resMAPlot",
    "VolcanoPlot", "DispersionPlot",
    "HeatmapPairwisePlot", "HeatmapDEGPlot",
]

########################################
# Helper: read SAMPLES after checkpoint
########################################
def get_samples(wildcards):
    """Input function that reads SRR.numbers only after download_fastq completes."""
    checkpoints.download_fastq.get(**wildcards)
    srr_file = "fastq_files/SRR.numbers"
    with open(srr_file) as f:
        return [line.strip() for line in f if line.strip()]


def all_transcript_csvs(wildcards):
    return expand("transcripts/{sample}.csv", sample=get_samples(wildcards))


########################################
# Setup Directories
########################################
rule make_directories:
    output:
        directory("data"),
        directory("logs"),
        directory("transcripts"),
        directory("fastq_files"),
        directory("ref_genome"),
        directory("alignment"),
    shell:
        # data/plots is created here so deg_analysis can write PNGs into it
        "mkdir -p data/plots logs transcripts fastq_files ref_genome alignment"


########################################
# Download FASTQ files (checkpoint)
########################################
checkpoint download_fastq:
    input:
        rules.make_directories.output,
    output:
        srr_numbers = "fastq_files/SRR.numbers",
    params:
        project = PROJECT,
    shell:
        """
        # Install EDirect if not already on PATH
        if ! command -v esearch &>/dev/null; then
            sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
            export PATH="${{HOME}}/edirect:${{PATH}}"
        fi

        # Fetch the list of SRR run accessions for this BioProject
        esearch -db sra -query {params.project} \
            | efetch -format runinfo \
            | cut -d ',' -f 1 \
            | tail -n +2 \
            > {output.srr_numbers}

        # Download all runs as single-end FASTQ files
        xargs -a {output.srr_numbers} -n 1 fasterq-dump -O fastq_files

        # Reject paired-end data: fasterq-dump produces _1/_2 suffixes for PE runs
        if ls fastq_files/*_2.fastq 2>/dev/null | grep -q .; then
            echo "ERROR: paired-end FASTQ files detected (*_2.fastq)." >&2
            echo "This pipeline only supports single-end RNA-Seq data." >&2
            exit 1
        fi
        """


########################################
# Download and Prepare Reference Genome
########################################
rule download_ref_fna:
    input:
        rules.make_directories.output,
    output:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
    shell:
        """
        curl -o {output} \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
        """

rule unzip_ref_fna:
    input:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz",
    output:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna",
    shell:
        "gunzip -c {input} > {output}"

rule download_gtf:
    input:
        rules.make_directories.output,
    output:
        "data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
    shell:
        """
        curl -o {output} \
            https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
        """

rule build_hisat2_index:
    input:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna",
    output:
        multiext("ref_genome/ref_gen",
                 ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                 ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"),
    threads: 32
    shell:
        "hisat2-build -p {threads} {input} ref_genome/ref_gen"


########################################
# Alignment and Counting (per sample)
########################################
rule align_reads:
    input:
        index = multiext("ref_genome/ref_gen",
                         ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                         ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"),
        fastq = "fastq_files/{sample}.fastq",
    output:
        "alignment/{sample}.sam",
    threads: 32
    shell:
        "hisat2 -p {threads} -x ref_genome/ref_gen -U {input.fastq} -S {output}"

rule count_reads:
    input:
        sam = "alignment/{sample}.sam",
        gtf = "data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz",
    output:
        "transcripts/{sample}.csv",
    threads: 32
    shell:
        """
        htseq-count -f sam -r pos -i gene_id -t exon -n {threads} \
            {input.sam} {input.gtf} > {output}
        """


########################################
# Post-processing: merge + DEG analysis
########################################
rule merge_results:
    input:
        all_transcript_csvs,
    output:
        "data/counts.csv",
    shell:
        """
        python3 merge_transcripts.py \
            --transcripts-dir transcripts \
            --out {output} \
            --ensembl-cache data/ensemble_df.pkl \
            --srr-map data/srr_to_gsm.tsv
        """

rule deg_analysis:
    input:
        counts = "data/counts.csv",
    output:
        results = "deg_results.csv",
        plots   = expand("data/plots/{name}.png", name=PLOT_NAMES),
    params:
        gse             = GSE_ACCESSION,
        condition_field = CONDITION_FIELD,
        plots_dir       = "data/plots",
    shell:
        """
        Rscript run_deg_analysis.R \
            --counts {input.counts} \
            --gse {params.gse} \
            --condition-field "{params.condition_field}" \
            --out {output.results} \
            --plots-dir {params.plots_dir}
        """


########################################
# Default target
########################################
rule all:
    input:
        all_transcript_csvs,
        "deg_results.csv",
        expand("data/plots/{name}.png", name=PLOT_NAMES),
