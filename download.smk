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