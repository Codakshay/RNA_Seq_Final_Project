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

# Test-mode knobs (default 0 = no effect — production runs ignore them)
SUBSAMPLE_READS   = int(config.get("subsample_reads", 0))  # fasterq-dump -X N
MAX_SAMPLES       = int(config.get("max_samples", 0))      # head -n N of SRR list

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
        ".dirs_ready"
    shell:
        "mkdir -p data/plots logs transcripts fastq_files ref_genome alignment && touch .dirs_ready"


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
        project         = PROJECT,
        parallel        = DOWNLOAD_PARALLEL,
        max_samples     = MAX_SAMPLES,
        subsample_reads = SUBSAMPLE_READS,
    threads: 32
    shell:
        r"""
        # Install EDirect if not already on PATH
        if ! command -v esearch &>/dev/null; then
            sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"
            export PATH="${{HOME}}/edirect:${{PATH}}"
        fi

        # Fetch SRR accessions — retry up to 5 times to handle transient SSL/network errors
        for attempt in 1 2 3 4 5; do
            esearch -db sra -query {params.project} \
                | efetch -format runinfo \
                | cut -d ',' -f 1 \
                | tail -n +2 \
                > {output.srr_numbers}.full
            [ -s {output.srr_numbers}.full ] && break
            echo "esearch attempt $attempt returned empty results, retrying in 10s..." >&2
            sleep 10
        done

        # Optional cap (test mode): take half from head, half from tail of the SRR
        # list. GEO submissions are almost always grouped by condition, so this
        # heuristic ensures both groups are represented in the test subset.
        if [ "{params.max_samples}" -gt 0 ]; then
            HALF=$(( {params.max_samples} / 2 ))
            OTHER=$(( {params.max_samples} - HALF ))
            head -n $HALF {output.srr_numbers}.full > {output.srr_numbers}
            tail -n $OTHER {output.srr_numbers}.full >> {output.srr_numbers}
            # Dedup in case head and tail overlap on a small full list
            awk '!seen[$0]++' {output.srr_numbers} > {output.srr_numbers}.dedup
            mv {output.srr_numbers}.dedup {output.srr_numbers}
        else
            mv {output.srr_numbers}.full {output.srr_numbers}
        fi
        rm -f {output.srr_numbers}.full

        if [ ! -s {output.srr_numbers} ]; then
            echo "ERROR: esearch returned no SRR accessions for {params.project}" >&2
            exit 1
        fi

        # Download FASTQs from ENA's pre-converted .fastq.gz mirror (much faster
        # than fasterq-dump, which downloads .sra and locally converts). The
        # ENA filereport API returns the exact FTP URLs and tells us paired vs
        # single (1 URL = single, 2 URLs = paired). Fallback to fastq-dump only
        # if ENA has no record for a given SRR.
        download_srr() {{
            local srr=$1
            local subsample=$2
            local urls
            urls=$(curl -sf --retry 3 \
                "https://www.ebi.ac.uk/ena/portal/api/filereport?accession=${{srr}}&result=read_run&fields=fastq_ftp&format=tsv" \
                | tail -n +2 | cut -f1 | tr ';' '\n' | grep -v '^$' || true)

            if [ -z "$urls" ]; then
                echo "[$srr] ENA has no FASTQ URLs; falling back to fastq-dump" >&2
                if [ "$subsample" -gt 0 ]; then
                    fastq-dump -X "$subsample" --split-3 -O fastq_files "$srr"
                else
                    fasterq-dump -e 4 --split-files -O fastq_files "$srr"
                fi
                return
            fi

            for url in $urls; do
                # ENA may return URLs with or without protocol — strip any leading
                # protocol so we control it ourselves and avoid 'https://https://...'.
                local clean_url dl_url fname out
                clean_url=$(echo "$url" | sed -E 's|^[a-z]+://||')
                dl_url="https://${{clean_url}}"
                fname=$(basename "$clean_url" .gz)
                out="fastq_files/$fname"

                # Sanity check: ENA filenames must end in .fastq
                case "$fname" in
                    *.fastq) ;;
                    *) echo "[$srr] WARNING: unexpected ENA filename '$fname' (url=$url); skipping" >&2
                       continue ;;
                esac

                echo "[$srr] Downloading $dl_url -> $out" >&2

                if [ "$subsample" -gt 0 ]; then
                    # Stream + truncate. SIGPIPE upstream (curl, gunzip) when head
                    # closes the pipe is expected — disable pipefail/errexit locally.
                    set +eo pipefail
                    curl -sf --retry 3 -L "$dl_url" 2>/dev/null \
                        | gunzip 2>/dev/null \
                        | head -n $(( subsample * 4 )) > "$out"
                    set -eo pipefail

                    if [ ! -s "$out" ]; then
                        echo "[$srr] ERROR: empty output for $dl_url" >&2
                        rm -f "$out"
                        return 1
                    fi
                else
                    if ! curl -sf --retry 3 -L "$dl_url" -o "${{out}}.gz"; then
                        echo "[$srr] ERROR: curl failed for $dl_url" >&2
                        rm -f "${{out}}.gz"
                        return 1
                    fi
                    gunzip -f "${{out}}.gz"
                    if [ ! -s "$out" ]; then
                        echo "[$srr] ERROR: gunzip produced empty $out" >&2
                        rm -f "$out"
                        return 1
                    fi
                fi
            done
        }}
        export -f download_srr

        xargs -r -a {output.srr_numbers} -n 1 -P {params.parallel} -I '{{SRR}}' \
            bash -c 'download_srr "$0" {params.subsample_reads}' '{{SRR}}'

        # Build layouts.tsv (SRR<TAB>single|paired) by inspecting the produced
        # files. Fail loudly if a SRR has no FASTQ at all (e.g. download failed)
        # rather than letting it slip through as a misclassified 'single'.
        : > {output.layouts}
        missing=0
        while read SRR; do
            if [ -f "fastq_files/${{SRR}}_2.fastq" ]; then
                echo -e "${{SRR}}\tpaired" >> {output.layouts}
            elif [ -f "fastq_files/${{SRR}}.fastq" ]; then
                echo -e "${{SRR}}\tsingle" >> {output.layouts}
            else
                echo "ERROR: no FASTQ produced for $SRR" >&2
                missing=$((missing + 1))
            fi
        done < {output.srr_numbers}
        if [ "$missing" -gt 0 ]; then
            echo "ERROR: $missing SRR(s) failed to download — aborting" >&2
            exit 1
        fi
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
        params:
             image = os.environ.get("IMAGE_NAME", "clara-parabricks_4.7.0-1.sif"),
        shell:
            # pbrun rna_fq2bam takes 1 or 2 paths after --in-fq; {input.fastq}
            # expands to a space-separated list, so both layouts Just Work.
            "apptainer exec --nv {params.image} "
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
