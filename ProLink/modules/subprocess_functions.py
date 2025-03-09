
import logging
import subprocess
import re
from .. import ProLink_path


logger = logging.getLogger()

# Al principio del archivo, después de las importaciones
def extract_species_name(description: str) -> str:
    '''
    Extrae solo el nombre de la especie de una cadena de descripción.
    '''
    match = re.search(r"_(\w+_\w+)_", description)  # Buscar el primer grupo de caracteres antes de un espacio
    if match:
        return match.group(1)  # Devuelve el primer grupo, que es el nombre de la especie
    return description  # Si no se encuentra, devuelve la descripción original


def modify_newick_with_species_only(newick_input: str, newick_output: str) -> None:
    '''
    Modifica el archivo Newick para que solo contenga el nombre de la especie en cada nodo.
    '''
    with open(newick_input, 'r') as file:
        newick_data = file.read()

    # Buscar todos los taxones en el árbol
    def replace_taxon(match):
        return extract_species_name(match.group(1))  # Extrae y reemplaza con solo el nombre de la especie

    modified_newick = re.sub(r"([A-Za-z0-9_.]+)", replace_taxon, newick_data)

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
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output]
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")
    mega_run = subprocess.run(mega_cmd)
    if mega_run.returncode != 0:
        logger.error(f"ERROR: MEGA-CC failed")
        raise RuntimeError(f"MEGA-CC failed")
    
    # Llamar a la función que modificará el archivo Newick
    newick_output = mega_output.replace('.tree', '_modified.tree')  # Crear un nuevo nombre de archivo para el árbol modificado
    modify_newick_with_species_only(mega_output, newick_output)
    logger.info(f"Modified tree saved as: {newick_output}")
