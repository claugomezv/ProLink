
import logging
import subprocess
import re
from .. import ProLink_path


logger = logging.getLogger()

def extract_species_name(description: str) -> str:
    '''
    Extrae solo el nombre de la especie eliminando códigos WP, nombre de la proteína, número de cluster y "Same Domains".
    '''
    match = re.search(r"([A-Za-z0-9_]+(?:_[A-Za-z0-9_]+)*)", description)  # Search for species name'
    if match:
        return match.group(1)  # Devuelve solo el nombre de la especie
    return description  # Si no encuentra, devuelve el texto original


def modify_newick_with_species_only(newick_input: str, newick_output: str) -> None:
    '''
    Modifica el archivo Newick para dejar solo el nombre de la especie en cada nodo.
    '''
    with open(newick_input, 'r') as file:
        newick_data = file.read()

    # Expresión regular para limpiar los nombres de los taxones
    modified_newick = re.sub(
        r"(WP_\d+\.\d+|MULTISPECIES:|alkene_reductase|---C\d+---|Same_Domains|[\"'])",
        "", 
        newick_data
    )
    
    # Aplicar la función para extraer solo el nombre de la especie
    modified_newick = re.sub(
        r"([A-Za-z0-9_]+)", 
        lambda match: extract_species_name(match.group(1)), 
        modified_newick
    )

    # Guardar el archivo Newick modificado
    with open(newick_output, 'w') as file:
        file.write(modified_newick)


    logger.info(f"✅ Árbol modificado guardado en: {newick_output}")



def align(muscle_input:str, muscle_output:str) -> None:
    '''
    Run a local alignment with MUSCLE v5

    Parameters
    ----------
    muscle_input : str
        Path of the input MUSCLE file
    muscle_output : str
        Path of the output MUSCLE file
    '''
    logging.info(f"\n-- Aligning sequences with MUSCLE")
    muscle_cmd = ['muscle', '-super5', muscle_input, '-output', muscle_output]
    logging.debug(f"Running MUSCLE alignment: {' '.join(muscle_cmd)}")
    muscle_run = subprocess.run(muscle_cmd)
    if muscle_run.returncode != 0:
        logger.error(f"ERROR: MUSCLE failed")
        raise RuntimeError(f"MUSCLE failed")

def tree(tree_type:str, bootstrap_replications:int, muscle_output:str, mega_output:str) -> None:
    '''
    Run MEGA-CC to generate a phylogenetic tree

    Parameters
    ----------
    tree_type : str
        Type of tree to generate
    bootstrap_replications : int
        Number of bootstrap replications
    muscle_output : str
        Path of the input file (FASTA format from MUSCLE)
    mega_output : str
        Path of the MEGA-CC output file
    '''
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    logging.info(f"\n-- Generating phylogenetic tree with MEGA-CC")
    
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output, '-log', 'mega_log.txt']
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")
    
    mega_run = subprocess.run(mega_cmd, capture_output=True, text=True)
    
    if mega_run.returncode != 0:
        logger.error(f"ERROR: MEGA-CC failed with error: {mega_run.stderr}")
        raise RuntimeError(f"MEGA-CC failed with error: {mega_run.stderr}")
    
    # Llamar a la función que modificará el archivo Newick
    newick_output = mega_output.replace('.tree', '_modified.tree')  # Crear un nuevo nombre de archivo para el árbol modificado
    modify_newick_with_species_only(mega_output, newick_output)
    logger.info(f"Modified tree saved as: {newick_output}")
