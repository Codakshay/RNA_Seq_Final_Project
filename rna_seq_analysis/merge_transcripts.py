import os
import pandas as pd
from biomart import BiomartServer
from loguru import logger
import time
import metadata

# Function to fetch data for a chunk of gene symbols
def fetch_chunk(chunk, dataset, sleep_time=5):
    try:
        logger.info("Going to sleep...")
        time.sleep(sleep_time)
        logger.info("Awake! Overwhelming the server...")
        response = dataset.search({
            'filters': {'hgnc_symbol': chunk},
            'attributes': ['ensembl_gene_id', 'hgnc_symbol', 'gene_biotype', 'chromosome_name']
        })
        logger.success("Successful retrieval from the server!")
        data = response.text.strip().split('\n')
        lines = [line.split('\t') for line in data if line.strip()]  # Skip empty lines
        return lines
    except Exception as e:
        logger.error(f"Error fetching chunk: {chunk} - {e}")
        exit(1)
    return []

def fetch_data_in_chunks(gene_symbols, chunk_size=500):
    all_results = []
    total_processing_units = len(gene_symbols)
    first_chunks = [gene_symbols[i:i + chunk_size] for i in range(0, len(gene_symbols), chunk_size)]
    chunks = []
    for _, chunk in enumerate(first_chunks):
        subdivision = len(' '.join(chunk))
        if subdivision >= 4000:
            subdivision = subdivision // 4000 + 1
            k, m = divmod(len(chunk), subdivision)
            for i in range(subdivision):
                chunks.append(chunk[i * k + min(i, m):(i + 1) * k + min(i + 1, m)])
        else:
            chunks.append(chunk)
    
    for _, chunk in enumerate(chunks):
        logger.info("Establishing the connection to the server...")
        server = BiomartServer("http://www.ensembl.org/biomart")
        logger.success("Connected to the server!")
        hsapiens = server.datasets['hsapiens_gene_ensembl']
        logger.info("Selected the database from biomart.")
        result = fetch_chunk(chunk, dataset=hsapiens)
        if result:
            all_results.extend(result)
        logger.info(f"{len(chunk)}/{total_processing_units} genes are processed for Ensembl ID matching.")
    
    return all_results


def generate_ensemble_dataframe():
    # Specify the directory
    directory = 'rna_seq_analysis/transcripts'

    # Get the list of all filenames in the directory
    filenames = os.listdir(directory)

    if len(filenames) == 0:
        logger.critical(f"{directory} is empty. Please check if the transcript is correctly generated.")
        exit(1)
    
    filename = filenames[0]  # Pick a CSV transcript data file
    full_path = os.path.join(directory, filename)
    logger.info(f"Generate the gene ID database from {full_path}...")
    dataframe = pd.read_csv(full_path, sep='\s+', header=None, names=["Gene", "Expression"])
    logger.success("Successfully generated the gene ID database.")
    
    logger.info("Querying the list of entries to drop from the database...")
    entries_to_drop = ['__no_feature', '__ambiguous', '__too_low_aQual', '__not_aligned', '__alignment_not_unique']

    logger.info("Dropping rows with specific entries...")
    dataframe = dataframe[~dataframe['Gene'].isin(entries_to_drop)]
    logger.success("Successfully normalized the database.")

    logger.info("Querying the list of gene symbols to lookup...")
    gene_symbols = dataframe['Gene'].to_list()[:1000]

    logger.info("Fetching data in chunks sequentially...")
    results = fetch_data_in_chunks(gene_symbols)
    logger.success("Successfully retrieved the data chunks from the BioMart database.")

    logger.info("Creating the final dataframe...")
    df = pd.DataFrame(results, columns=["Ensembl Gene ID", "Gene Symbol", "Biotype", "Chromosome"])

    # save the dataframe as a pickle file
    df.to_pickle("rna_seq_analysis/ensemble_df.pkl")
    return 0

def generate_count_dataframe():
    # Specify the directory
    directory = 'rna_seq_analysis/transcripts'
    # Get the list of all filenames in the directory
    filenames = os.listdir(directory)
    ensemble_df = pd.read_pickle("rna_seq_analysis/ensemble_df.pkl")
    for filename in filenames:
        full_path = os.path.join(directory, filename)
        geo_id = metadata.get_geo_accession_from_srr(filename[:filename.index(".csv")])
        transcript_dataframe = pd.read_csv(full_path, sep='\s+', header=None, names=["Gene Symbol", geo_id])
        ensemble_df = ensemble_df.merge(transcript_dataframe, how='inner', on="Gene Symbol")
    # save complete count dataframe into csv
    ensemble_df.to_csv("rna_seq_analysis/counts.csv")
    return 0


if __name__ == "__init__":
    generate_ensemble_dataframe()
    generate_count_dataframe()