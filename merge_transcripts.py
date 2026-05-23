"""Merge per-sample htseq-count CSV files into an annotated count matrix.

Usage (called by Snakemake rule merge_results):
    python3 merge_transcripts.py \
        --transcripts-dir transcripts \
        --out data/counts.csv \
        --ensembl-cache data/ensemble_df.pkl \
        --srr-map data/srr_to_gsm.tsv
"""
import argparse
import os
import pandas as pd
from biomart import BiomartServer
from loguru import logger
import time
import metadata

# HTSeq-count appends these special summary lines at the end of each output file.
_HTSEQ_SPECIAL = frozenset({
    '__no_feature', '__ambiguous', '__too_low_aQual',
    '__not_aligned', '__alignment_not_unique',
})


def _parse_args():
    p = argparse.ArgumentParser(
        description='Merge per-sample transcript counts into an annotated count matrix.'
    )
    p.add_argument('--transcripts-dir', default='transcripts',
                   help='Directory containing per-sample count CSV files (default: transcripts)')
    p.add_argument('--out', default='data/counts.csv',
                   help='Output merged counts CSV (default: data/counts.csv)')
    p.add_argument('--ensembl-cache', default='data/ensemble_df.pkl',
                   help='Pickle cache path for BioMart Ensembl lookup (default: data/ensemble_df.pkl)')
    p.add_argument('--srr-map', default='data/srr_to_gsm.tsv',
                   help='TSV cache path for SRR→GSM mapping (default: data/srr_to_gsm.tsv)')
    return p.parse_args()


# ---------------------------------------------------------------------------
# BioMart helpers
# ---------------------------------------------------------------------------

def _fetch_chunk(chunk, dataset, sleep_time=5):
    """Query BioMart for a list of HGNC symbols.  Returns a list of row-lists."""
    logger.info(f"Sleeping {sleep_time}s before BioMart query ({len(chunk)} symbols)...")
    time.sleep(sleep_time)
    response = dataset.search({
        'filters': {'hgnc_symbol': chunk},
        'attributes': ['ensembl_gene_id', 'hgnc_symbol', 'gene_biotype', 'chromosome_name'],
    })
    logger.success("BioMart query successful.")
    lines = response.text.strip().split('\n')
    return [line.split('\t') for line in lines if line.strip()]


def _fetch_data_in_chunks(gene_symbols, chunk_size=500):
    """Fetch BioMart annotation for all gene_symbols, batching to avoid URL limits."""
    all_results = []
    total = len(gene_symbols)

    # First pass: split into chunks of chunk_size
    raw_chunks = [gene_symbols[i:i + chunk_size] for i in range(0, total, chunk_size)]

    # Second pass: further sub-divide any chunk whose joined length exceeds ~4 000 chars
    chunks = []
    for chunk in raw_chunks:
        joined_len = len(' '.join(chunk))
        if joined_len >= 4000:
            n = joined_len // 4000 + 1
            k, m = divmod(len(chunk), n)
            for i in range(n):
                sub = chunk[i * k + min(i, m):(i + 1) * k + min(i + 1, m)]
                if sub:
                    chunks.append(sub)
        else:
            chunks.append(chunk)

    processed = 0
    for chunk in chunks:
        try:
            logger.info("Connecting to BioMart server...")
            server = BiomartServer("http://www.ensembl.org/biomart")
            hsapiens = server.datasets['hsapiens_gene_ensembl']
            result = _fetch_chunk(chunk, hsapiens)
            all_results.extend(result)
        except Exception as exc:
            logger.error(
                f"BioMart fetch failed for chunk starting with '{chunk[0]}': {exc}. "
                "Skipping chunk — those genes will be absent from the count matrix."
            )
        processed += len(chunk)
        logger.info(f"BioMart: {processed}/{total} gene symbols processed.")

    return all_results


# ---------------------------------------------------------------------------
# Core logic
# ---------------------------------------------------------------------------

def build_ensemble_df(transcripts_dir, cache_path):
    """Return a DataFrame (Ensembl_ID, GeneSymbol, Biotype, Chromosome) for all
    gene symbols found in transcripts_dir.  Result is pickled to cache_path.
    """
    if os.path.exists(cache_path):
        logger.info(f"Loading Ensembl annotation cache from {cache_path}.")
        return pd.read_pickle(cache_path)

    csv_files = [f for f in os.listdir(transcripts_dir) if f.endswith('.csv')]
    if not csv_files:
        logger.critical(f"{transcripts_dir} contains no CSV files.")
        raise FileNotFoundError(f"No CSV files in {transcripts_dir}")

    # Read one representative file to obtain the gene-symbol list
    sample_path = os.path.join(transcripts_dir, csv_files[0])
    logger.info(f"Building gene list from {sample_path}.")
    df = pd.read_csv(sample_path, sep='\t', header=None, names=['Gene', 'Expression'])
    df = df[~df['Gene'].isin(_HTSEQ_SPECIAL)]
    gene_symbols = df['Gene'].dropna().unique().tolist()
    logger.info(f"Found {len(gene_symbols)} unique gene symbols; querying BioMart...")

    results = _fetch_data_in_chunks(gene_symbols)
    ensemble_df = pd.DataFrame(
        results,
        columns=['Ensembl_ID', 'GeneSymbol', 'Biotype', 'Chromosome'],
    )
    # Keep the first Ensembl ID when a symbol maps to multiple genes
    ensemble_df = ensemble_df.drop_duplicates(subset='GeneSymbol', keep='first')

    os.makedirs(os.path.dirname(os.path.abspath(cache_path)), exist_ok=True)
    ensemble_df.to_pickle(cache_path)
    logger.success(f"Ensembl cache written to {cache_path} ({len(ensemble_df)} genes).")
    return ensemble_df


def build_count_matrix(transcripts_dir, ensemble_df, srr_to_gsm):
    """Inner-join each sample's counts onto the ensemble annotation DataFrame.

    Columns in the result: Ensembl_ID, GeneSymbol, Biotype, Chromosome, <GSM_ids…>
    """
    merged = ensemble_df.copy()

    for filename in sorted(os.listdir(transcripts_dir)):
        if not filename.endswith('.csv'):
            continue
        srr = filename[:-4]  # strip .csv
        gsm = srr_to_gsm.get(srr)
        if gsm is None:
            logger.warning(f"No GSM mapping for {srr}; skipping sample.")
            continue

        full_path = os.path.join(transcripts_dir, filename)
        sample_df = pd.read_csv(full_path, sep='\t', header=None,
                                names=['GeneSymbol', gsm])
        sample_df = sample_df[~sample_df['GeneSymbol'].isin(_HTSEQ_SPECIAL)]
        merged = merged.merge(sample_df, how='inner', on='GeneSymbol')
        logger.info(f"Merged {srr} ({gsm}): {len(merged)} genes remaining.")

    return merged


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

if __name__ == '__main__':
    args = _parse_args()

    # 1. Resolve SRR → GSM for every sample in the transcripts directory
    srr_list = sorted(
        f[:-4] for f in os.listdir(args.transcripts_dir) if f.endswith('.csv')
    )
    logger.info(f"Found {len(srr_list)} samples in {args.transcripts_dir}.")
    srr_to_gsm = metadata.build_srr_to_gsm_map(srr_list, args.srr_map)

    # 2. Fetch (or load cached) BioMart annotation for all gene symbols
    ensemble_df = build_ensemble_df(args.transcripts_dir, args.ensembl_cache)

    # 3. Merge per-sample counts into one matrix
    counts = build_count_matrix(args.transcripts_dir, ensemble_df, srr_to_gsm)
    n_samples = counts.shape[1] - 4  # subtract the 4 annotation columns
    logger.info(f"Final count matrix: {counts.shape[0]} genes × {n_samples} samples.")

    # 4. Write output (no pandas index — the Ensembl_ID column serves as the key)
    out_dir = os.path.dirname(os.path.abspath(args.out))
    os.makedirs(out_dir, exist_ok=True)
    counts.to_csv(args.out, index=False)
    logger.success(f"Count matrix written to {args.out}.")
