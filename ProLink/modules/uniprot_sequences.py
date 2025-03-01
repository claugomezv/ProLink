import requests
import re
import logging
from Bio import SeqIO  # To properly handle FASTA files

logger = logging.getLogger()

def check_uniprot_batch(wp_codes):
    """
    Verify the existence of multiple WP codes in UniProt in a single request.

    Parameters:
    wp_codes (list): List of WP codes to verify.

    Returns:
    set: Set of WP codes found in UniProt.
    """
    url = "https://rest.uniprot.org/uniprotkb/search"
    # Prepend each WP code with "accession:" so that the search looks in the accession field
    query = " OR ".join(f"xref_refseq:{wp}" for wp in wp_codes)
    params = {
        "query": query,
        "fields": "accession, xref_refseq", #Include xref_refseq to see the matches RefSeq IDs
        "format": "json",
        "size": len(wp_codes)  # Ensure we get as many results as the number of queries in the batch
    }

    print(f"Consulta a UniProt: {query}")  # Debug: check what is being sent to UniProt

    try:
        response = requests.get(url, params=params)
        response.raise_for_status()
        data = response.json()
        valid_entries = {entry['primaryAccession'] for entry in data.get("results", [])}
        print(f"Códigos WP encontrados en UniProt: {valid_entries}")  # Debug: show WP codes returned by UniProt
        return valid_entries
    except Exception as e:
        logger.error(f"Error al conectar con UniProt: {e}")
        return set()

def filter_valid_sequences(input_fasta, output_fasta):
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
        match = re.search(r'WP_\d+\.\d+', seq.description)
        #logger.info(f"Código WP encontrado en {seq.id}: {match.group(1) if match else 'Ninguno'}")
        if match:
            wp_data[seq.id] = match.group(1)
    
    print(f"Códigos WP extraídos: {list(wp_data.values())}")  # Debug: show extracted WP codes
    
    logger.info(f"Número total de secuencias: {len(sequences)}")
    logger.info(f"Número de códigos WP encontrados: {len(wp_data)}")

    # Verify in UniProt in batches (using batch size of 100)
    wp_codes = list(set(wp_data.values()))  # Remove duplicate codes
    logger.info(f"Total códigos WP extraídos: {len(wp_codes)}")
    #logger.info(f"Códigos WP únicos a consultar: {wp_codes}")
    valid_wp_codes = set()
    batch_size = 100  

    for i in range(0, len(wp_codes), batch_size):
        batch = wp_codes[i:i+batch_size]
        valid_wp_codes.update(check_uniprot_batch(batch))

    # Filter valid sequences
    valid_sequences = [
        seq for seq in sequences 
        if (seq.id in wp_data and wp_data[seq.id] in valid_wp_codes) or seq.id not in wp_data
    ]
   
    #logger.info(f"Secuencias válidas encontradas después del filtrado: {[seq.id for seq in valid_sequences]}")

    # Write the valid sequences to the new FASTA file
    SeqIO.write(valid_sequences, output_fasta, "fasta")
    #logger.info(f"Archivo filtrado guardado en {output_fasta} con {len(valid_sequences)} secuencias.")

    print(f"Secuencias válidas después del filtrado: {len(valid_sequences)}")  # Debug: show number of valid sequences
    logger.info(f"Resultados guardados en {output_fasta}")
