configfile: "config.yaml"
PROJECT = config["project"]

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
        directory("rna_seq_analysis/data")
    shell:
        """
        mkdir -p data logs transcripts fastq_files ref_genome alignment rna_seq_analysis/data
        """

rule download_fastq:
    input:
        PROJECT
    output:
        "fastq_files/SRR.numbers"
    params:
        sra_toolkit_module="sra-toolkit/3.0.9"
    shell:
        """
        module load {params.sra_toolkit_module}
        
        sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
        export PATH=${{HOME}}/edirect:${{PATH}}
        
        esearch -db sra -query {input} | efetch -format runinfo | cut -d "," -f 1 | tail -n +2 > {output}
        
        xargs < {output} -n 1 fasterq-dump -O fastq_files
        """


########################################
# Download and Prepare Reference Genome and Annotation
########################################
rule download_ref_fna:
    input:
        rules.make_directories.output
    output:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    shell:
        """
        curl -o {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
        """

rule unzip_ref_fna:
    input:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
    output:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna"
    shell:
        "gunzip -c {input} > {output}"

rule download_gtf:
    input:
        rules.make_directories.output
    output:
        "data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
    shell:
        """
        curl -o {output} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
        """

SAMPLES = [line.strip() for line in open("fastq_files/SRR.numbers") if line.strip()]

rule all:
    input:
        expand("transcripts/{sample}.csv", sample=SAMPLES),
        "deg_results.csv"

rule build_hisat2_index:
    input:
        "data/GCF_000001405.40_GRCh38.p14_genomic.fna"
    output:
        "ref_genome/ref_gen.1.ht2",
        "ref_genome/ref_gen.2.ht2",
        "ref_genome/ref_gen.3.ht2",
        "ref_genome/ref_gen.4.ht2",
        "ref_genome/ref_gen.5.ht2",
        "ref_genome/ref_gen.6.ht2",
        "ref_genome/ref_gen.7.ht2",
        "ref_genome/ref_gen.8.ht2"
    threads: 32
    shell:
        """
        hisat2-build -p {threads} {input} ref_genome/ref_gen
        """

########################################
# Alignment and Counting
########################################
rule align_reads:
    input:
        index_prefix = "ref_genome/ref_gen",
        fastq = "fastq_files/{sample}.fastq"
    output:
        "alignment/{sample}.sam"
    threads: 32
    shell:
        """
        hisat2 -p {threads} -x {input.index_prefix} -U {input.fastq} -S {output}
        """

rule count_reads:
    input:
        sam = "alignment/{sample}.sam",
        gtf = "data/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz"
    output:
        "transcripts/{sample}.csv"
    threads: 32
    shell:
        """
        htseq-count {input.sam} {input.gtf} -f sam -r pos -i gene_id -t exon -n {threads} > {output}
        """

########################################
# Post-processing and DEG Analysis
########################################
rule merge_results:
    input:
        expand("transcripts/{sample}.csv", sample=SAMPLES)
    output:
        "rna_seq_analysis/data/counts.csv"
    shell:
        """
        python3 merge_script.py > {output}
        """

rule deg_analysis:
    input:
        "rna_seq_analysis/data/counts.csv"
    output:
        "deg_results.csv"
    shell:
        """
        Rscript run_deg_analysis.R
        """