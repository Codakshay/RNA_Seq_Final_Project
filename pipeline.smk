import os

configfile: "config.yaml"
PROJECT         = config["project"]
GSE_ACCESSION   = config["gse_accession"]
CONDITION_FIELD = config["condition_field"]

# Optional execution-mode keys (added in Stage 8 — safe defaults preserve old behaviour)
LOCAL_FA          = config.get("local_reference_fa", "").strip()
LOCAL_GTF         = config.get("local_reference_gtf", "").strip()
DOWNLOAD_PARALLEL = int(config.get("download_parallel", 4))
ALIGNER           = config.get("aligner", "hisat2").strip()
assert ALIGNER in ("hisat2", "parabricks_star"), \
    f"config[aligner] must be 'hisat2' or 'parabricks_star', got {ALIGNER!r}"
DESEQ_WORKERS     = int(config.get("deseq_workers", 1))

# Plot filenames that run_deg_analysis.R will produce
PLOT_NAMES = [
    "PCAPlot", "MAPlot", "resMAPlot",
    "VolcanoPlot", "DispersionPlot",
    "HeatmapPairwisePlot", "HeatmapDEGPlot",
]

########################################
# Helpers: read SAMPLES + layouts after the download checkpoint
########################################
def get_samples(wildcards):
    """Input function that reads SRR.numbers only after download_fastq completes."""
    checkpoints.download_fastq.get(**wildcards)
    srr_file = "fastq_files/SRR.numbers"
    with open(srr_file) as f:
        return [line.strip() for line in f if line.strip()]


def all_transcript_csvs(wildcards):
    return expand("transcripts/{sample}.csv", sample=get_samples(wildcards))


def _read_layouts():
    """Parse fastq_files/layouts.tsv (SRR<TAB>single|paired). Returns {} if missing."""
    path = "fastq_files/layouts.tsv"
    if not os.path.exists(path):
        return {}
    out = {}
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            srr, layout = line.split("\t", 1)
            out[srr] = layout
    return out


def is_paired(sample):
    """True iff layouts.tsv marks this sample as paired-end."""
    return _read_layouts().get(sample, "single") == "paired"


def fastq_inputs(wildcards):
    """Return the FASTQ file(s) for {sample}, gated on the download checkpoint."""
    checkpoints.download_fastq.get()
    if is_paired(wildcards.sample):
        return [f"fastq_files/{wildcards.sample}_1.fastq",
                f"fastq_files/{wildcards.sample}_2.fastq"]
    return [f"fastq_files/{wildcards.sample}.fastq"]


def hisat2_read_args(wildcards):
    """Return the HISAT2 read-input fragment (-U single OR -1 r1 -2 r2)."""
    if is_paired(wildcards.sample):
        return (f"-1 fastq_files/{wildcards.sample}_1.fastq "
                f"-2 fastq_files/{wildcards.sample}_2.fastq")
    return f"-U fastq_files/{wildcards.sample}.fastq"


def featurecounts_paired_flag(wildcards):
    """featureCounts paired-end flag, empty for single-end."""
    return "-p --countReadPairs" if is_paired(wildcards.sample) else ""


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
# Download FASTQ files (checkpoint, parallel, paired-end aware)
########################################
checkpoint download_fastq:
    input:
        rules.make_directories.output,
    output:
        srr_numbers = "fastq_files/SRR.numbers",
        layouts     = "fastq_files/layouts.tsv",
    params:
        project    = PROJECT,
        parallel   = DOWNLOAD_PARALLEL,
    threads: 32
    shell:
        r"""
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

        # Parallel download. With {params.parallel} concurrent processes each
        # using 8 worker threads, we keep all 32 cores busy AND cut wall time
        # to roughly 1/{params.parallel} of the serial baseline (network/NCBI-
        # throttling permitting).
        xargs -a {output.srr_numbers} -n 1 -P {params.parallel} \
            fasterq-dump -e 8 --split-files -O fastq_files

        # Build layouts.tsv (SRR<TAB>single|paired) by inspecting the produced files.
        : > {output.layouts}
        while read SRR; do
            if [ -f "fastq_files/${{SRR}}_2.fastq" ]; then
                echo -e "${{SRR}}\tpaired" >> {output.layouts}
            else
                echo -e "${{SRR}}\tsingle" >> {output.layouts}
            fi
        done < {output.srr_numbers}
        """


########################################
# Reference genome + annotation
#
# When config[local_reference_fa] / [local_reference_gtf] are set, the
# corresponding download rule is replaced by a symlink-or-decompress rule.
# Otherwise the original NCBI GRCh38 GCF files are fetched.
########################################
if LOCAL_FA:
    rule provide_reference_fa:
        input:
            rules.make_directories.output,
        output:
            "data/reference.fna",
        params:
            src = LOCAL_FA,
        shell:
            r"""
            if [[ "{params.src}" == *.gz ]]; then
                gunzip -c "{params.src}" > {output}
            else
                ln -sf "$(readlink -f '{params.src}')" {output}
            fi
            """
else:
    rule download_reference_fa:
        input:
            rules.make_directories.output,
        output:
            "data/reference.fna",
        shell:
            """
            curl -o {output}.gz \
                https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz
            gunzip -f {output}.gz
            """

if LOCAL_GTF:
    rule provide_reference_gtf:
        input:
            rules.make_directories.output,
        output:
            "data/reference.gtf.gz",
        params:
            src = LOCAL_GTF,
        shell:
            r"""
            if [[ "{params.src}" == *.gz ]]; then
                ln -sf "$(readlink -f '{params.src}')" {output}
            else
                gzip -c "{params.src}" > {output}
            fi
            """
else:
    rule download_reference_gtf:
        input:
            rules.make_directories.output,
        output:
            "data/reference.gtf.gz",
        shell:
            """
            curl -o {output} \
                https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.gtf.gz
            """

########################################
# Aligner index (only the one matching config[aligner] is materialized)
########################################
if ALIGNER == "hisat2":
    rule build_hisat2_index:
        input:
            "data/reference.fna",
        output:
            multiext("ref_genome/ref_gen",
                     ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                     ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"),
        threads: 32
        shell:
            "hisat2-build -p {threads} {input} ref_genome/ref_gen"

else:  # parabricks_star
    rule build_star_index:
        # STAR genomeGenerate is CPU-only (no Parabricks step), runs once, ~30 min on
        # 32 cores. The index ends up in ref_genome/star_index/ and is reused for all
        # alignment jobs.
        input:
            fa     = "data/reference.fna",
            gtf_gz = "data/reference.gtf.gz",
        output:
            directory("ref_genome/star_index"),
        threads: 32
        shell:
            r"""
            mkdir -p {output}
            # STAR's --sjdbGTFfile wants an uncompressed GTF
            TMPGTF=$(mktemp --suffix=.gtf)
            gunzip -c {input.gtf_gz} > "$TMPGTF"
            STAR --runMode genomeGenerate --runThreadN {threads} \
                 --genomeDir {output} \
                 --genomeFastaFiles {input.fa} \
                 --sjdbGTFfile "$TMPGTF" \
                 --sjdbOverhang 99
            rm -f "$TMPGTF"
            """


########################################
# Alignment (the rule actually instantiated depends on ALIGNER)
########################################
if ALIGNER == "hisat2":
    rule align_reads:
        # HISAT2 pipes SAM into samtools sort so downstream counting always sees
        # sorted BAM, regardless of which aligner produced it.
        input:
            index = multiext("ref_genome/ref_gen",
                             ".1.ht2", ".2.ht2", ".3.ht2", ".4.ht2",
                             ".5.ht2", ".6.ht2", ".7.ht2", ".8.ht2"),
            fastq = fastq_inputs,
        output:
            "alignment/{sample}.bam",
        threads: 32
        params:
            reads = hisat2_read_args,  # "-U …" for single, "-1 … -2 …" for paired
        shell:
            "hisat2 -p {threads} -x ref_genome/ref_gen {params.reads} "
            "| samtools sort -@ {threads} -o {output} -"

else:  # parabricks_star
    rule align_parabricks:
        # GPU alignment via NVIDIA Parabricks (rna_fq2bam uses STAR under the hood).
        # Requires `module load parabricks` and a SLURM allocation with --gres=gpu:1
        # (set up by run_pipeline.sh when invoked with --gpu).
        input:
            fastq = fastq_inputs,
            idx   = "ref_genome/star_index",
            fa    = "data/reference.fna",
        output:
            "alignment/{sample}.bam",
        threads: 16
        resources:
            gpu    = 1,
            mem_mb = 80000,
        shell:
            # pbrun rna_fq2bam takes 1 or 2 paths after --in-fq; {input.fastq}
            # expands to a space-separated list, so both layouts Just Work.
            "pbrun rna_fq2bam "
            "--ref {input.fa} "
            "--genome-lib-dir {input.idx} "
            "--in-fq {input.fastq} "
            "--out-bam {output}"


########################################
# Counting (input is the canonical sorted BAM from whichever aligner ran)
########################################
rule count_reads:
    # featureCounts is ~10-20x faster than htseq-count on the same input;
    # output is reduced to the (gene_id<TAB>count) two-column TSV that
    # merge_transcripts.py already knows how to parse.
    input:
        bam = "alignment/{sample}.bam",
        gtf = "data/reference.gtf.gz",
    output:
        counts  = "transcripts/{sample}.csv",
    threads: 8
    params:
        paired_flag = featurecounts_paired_flag,  # "-p --countReadPairs" iff paired
    shell:
        """
        featureCounts -T {threads} -a {input.gtf} \
                      -t exon -g gene_id {params.paired_flag} \
                      -o {output.counts}.raw {input.bam}
        # Strip featureCounts header lines (1: comment, 2: column names)
        # and keep only Geneid (col 1) and the per-sample count (last col).
        awk 'NR>2 {{print $1"\\t"$NF}}' {output.counts}.raw > {output.counts}
        rm {output.counts}.raw {output.counts}.raw.summary
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
        workers         = DESEQ_WORKERS,
    shell:
        """
        Rscript run_deg_analysis.R \
            --counts {input.counts} \
            --gse {params.gse} \
            --condition-field "{params.condition_field}" \
            --out {output.results} \
            --plots-dir {params.plots_dir} \
            --workers {params.workers}
        """


########################################
# Default target
########################################
rule all:
    input:
        all_transcript_csvs,
        "deg_results.csv",
        expand("data/plots/{name}.png", name=PLOT_NAMES),
