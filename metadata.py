# this script deals with functionalities regarding the metadata file
from Bio import Entrez
import re
from loguru import logger

# Provide your email to access NCBI
Entrez.email = "dzk5572@gmail.com"

# Function to read lines from a file and save them as a list
def read_file_to_list(filename):
    with open(filename, 'r') as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]  # Remove newline characters
    return lines

# Function to get GEO accession ID from SRR number
def get_geo_accession_from_srr(srr_number):
    # Query SRA database with SRR number
    handle = Entrez.esearch(db="sra", term=srr_number)
    record = Entrez.read(handle)
    handle.close()

    # Get SRA ID
    sra_id = record["IdList"][0]

    # Fetch details from SRA database
    handle = Entrez.efetch(db="sra", id=sra_id, rettype="docsum")
    sra_record = Entrez.read(handle)
    handle.close()

    # Extract the ExpXml content
    xml_content = sra_record[0]['ExpXml']

    # Search for GSM ID
    match = re.search(r'GSM\d+', xml_content)
    
    if match:
        result = match.group(0)
    else:
        logger.warning("No match is found. Please check if the correct metadata file is provided.")
        exit(1)
    
    return result