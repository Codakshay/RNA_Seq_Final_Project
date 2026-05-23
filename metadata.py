# Metadata utilities: SRR accession → GEO GSM accession lookups with TSV caching.
from Bio import Entrez
import re
import os
import csv
import time
from loguru import logger

Entrez.email = "dzk5572@gmail.com"


def read_file_to_list(filename):
    with open(filename, 'r') as f:
        return [line.strip() for line in f.readlines()]


def get_geo_accession_from_srr(srr_number):
    """Return the GSM accession for a single SRR run ID.

    Raises ValueError if the SRA record cannot be found or contains no GSM ID.
    """
    handle = Entrez.esearch(db="sra", term=srr_number)
    record = Entrez.read(handle)
    handle.close()

    if not record["IdList"]:
        raise ValueError(f"No SRA records found for {srr_number}")

    sra_id = record["IdList"][0]

    handle = Entrez.efetch(db="sra", id=sra_id, rettype="docsum")
    sra_record = Entrez.read(handle)
    handle.close()

    xml_content = sra_record[0]['ExpXml']
    match = re.search(r'GSM\d+', xml_content)
    if match:
        return match.group(0)

    raise ValueError(f"No GSM ID found in SRA record for {srr_number}")


def build_srr_to_gsm_map(srr_list, cache_path):
    """Build and cache a {SRR: GSM} mapping, querying Entrez only for missing entries.

    Reads an existing TSV cache at cache_path (columns: SRR<tab>GSM), fetches
    any SRRs not yet cached, writes the updated cache back, and returns a dict
    limited to the requested srr_list.  SRRs that cannot be resolved are logged
    and omitted from the returned dict.
    """
    cache = {}

    if os.path.exists(cache_path):
        with open(cache_path, newline='') as f:
            for row in csv.reader(f, delimiter='\t'):
                if len(row) == 2:
                    cache[row[0]] = row[1]
        logger.info(f"Loaded {len(cache)} cached SRR→GSM entries from {cache_path}.")

    missing = [srr for srr in srr_list if srr not in cache]
    if missing:
        logger.info(f"Querying Entrez for {len(missing)} uncached SRR IDs.")

    for srr in missing:
        try:
            gsm = get_geo_accession_from_srr(srr)
            cache[srr] = gsm
            logger.info(f"Mapped {srr} → {gsm}")
        except Exception as exc:
            logger.warning(f"Could not resolve GSM for {srr}: {exc}. Skipping.")
        # Stay within NCBI's 3 requests/second guideline
        time.sleep(0.34)

    cache_dir = os.path.dirname(os.path.abspath(cache_path))
    os.makedirs(cache_dir, exist_ok=True)
    with open(cache_path, 'w', newline='') as f:
        writer = csv.writer(f, delimiter='\t')
        for srr, gsm in cache.items():
            writer.writerow([srr, gsm])
    logger.success(f"SRR→GSM cache written to {cache_path} ({len(cache)} total entries).")

    return {srr: cache[srr] for srr in srr_list if srr in cache}
