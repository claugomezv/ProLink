import requests
import re
import logging
from Bio import SeqIO  # To properly handle FASTA files
from Bio import ExPASy, SwissProt

logger = logging.getLogger()

def check_uniprot_single(wp_code):
    """
    Verify the existence of a single WP code in UniProt.
    
    Parameters:
    wp_code (str): WP code to verify.

    Returns:
    bool: True if the WP code exists in UniProt, False otherwise.
    """
    # Primero se intenta obtener el registro desde Swiss-Prot
    try:
        handle = ExPASy.get_sprot_raw(wp_code)
        record = SwissProt.read(handle)
        return True
    except Exception as e:
        logger.info(f"{wp_code} no encontrado en Swiss‑Prot, intentando en UniProtKB REST API.")
        rest_url = f"https://rest.uniprot.org/uniprotkb/{wp_code}.json"
        params = {
            "query": f"xref:RefSeq-{wp_code}",
            "fields": "accession",
            "format": "json",
            "size": 1  # We only need to check if it exists
        }
        try:
            response = requests.get(rest_url, params=params, timeout=10)
            response.raise_for_status()
            data = response.json()
            return bool(data.get("results"))
                
        except Exception as e2:
            logger.error(f"Error al recuperar el registro para {wp_code} en UniProtKB: {e2}")
            return False

def filter_valid_sequences(input_fasta, output_fasta,  block_size=100):
    """
    Filters sequences by removing those whose WP codes do not exist in UniProt.
    Sequences without a WP_ code are retained.

    Parameters:
    input_fasta (str): Input FASTA file with sequences.
    output_fasta (str): Output FASTA file with valid sequences.
    """
    sequences = list(SeqIO.parse(input_fasta, "fasta"))

    # Extract WP_ codes from sequence descriptions
    wp_data = {}
    for seq in sequences:
        match = re.search(r'(WP_\d{9}\.\d)', seq.description)
        if match:
            wp_data[seq.description] = match.group(1)
    
    #print(f"Códigos WP extraídos: {list(wp_data.values())}")
    logger.info(f"Número total de secuencias: {len(sequences)}")
    logger.info(f"Número de códigos WP encontrados: {len(wp_data)}")

    # Procesar los códigos WP en bloques
    all_wp_codes = list(set(wp_data.values()))
    valid_wp_codes = set()
    for i in range(0, len(all_wp_codes), block_size):
        block = all_wp_codes[i:i+block_size]
        valid_codes_block = {code for code in block if check_uniprot_single(code)}
        valid_wp_codes.update(valid_codes_block)

    # Filter valid sequences
    valid_sequences = [
        seq for seq in sequences 
        if seq.description not in wp_data or wp_data[seq.description] in valid_wp_codes
    ]
   
    # Write the valid sequences to the new FASTA file
    SeqIO.write(valid_sequences, output_fasta, "fasta")
    print(f"Secuencias válidas después del filtrado: {len(valid_sequences)}")  # Debug: Show number of valid sequences
    logger.info(f"Resultados guardados en {output_fasta}")
