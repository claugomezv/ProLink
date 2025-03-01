import requests
from Bio import SeqIO

def check_uniprot_entries(wp_codes):
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    results = {}

    for wp in wp_codes:
        query = f"xref:RefSeq-{wp}"
        params = {"query": query, "fields": "accession", "format": "json"}
        response = requests.get(base_url, params=params)

        if response.status_code == 200:
            data = response.json()
            results[wp] = bool(data.get("results"))  # True si tiene entrada, False si no
        else:
            results[wp] = None  # Error en la consulta

    return results

def filter_sequences_by_uniprot(input_fasta, output_fasta):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    wp_codes = [seq.id.split()[0] for seq in sequences]

    valid_entries = check_uniprot_entries(wp_codes)

    # Filtrar secuencias que tienen entrada en UniProt
    filtered_sequences = [seq for seq in sequences if valid_entries.get(seq.id.split()[0], False)]

    # Guardar nuevo FASTA filtrado
    SeqIO.write(filtered_sequences, output_fasta, "fasta")
    print(f"Filtrado completado: {len(filtered_sequences)} secuencias guardadas en {output_fasta}")
