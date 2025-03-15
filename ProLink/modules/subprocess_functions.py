import logging
import subprocess
import time
import os
import re
from .. import ProLink_path
from Bio import Phylo
import io

logger = logging.getLogger()

def clean_label(label, protein_name='alkene_reductase'):
    """
    Cleans a single Newick label by extracting only the species name and the cluster marker.
    It removes unwanted substrings such as WP codes, "MULTISPECIES:", the protein name,
    and the word "unclassified". Then it extracts the species and the cluster marker.
    For example, from:
      WP_072607337.1_alkene_reductase_Aquibium_oceanicum_---C51---Same_Domains
    it returns:
      Aquibium_oceanicum_---C51
    """
    # Remove any surrounding quotes
    label = label.strip("'\"")
    # Remove WP codes
    label = re.sub(r"WP[\s_]+\d{9}\.\d", "", label, flags=re.IGNORECASE)
    # Remove "MULTISPECIES:" if present
    label = re.sub(r"MULTISPECIES:\s*", "", label, flags=re.IGNORECASE)
     # Remove any protein name that ends with "reductase" (e.g., "alkene_reductase", "oxidoreductase")
    label = re.sub(r"\b\w*reductase\b", "", label, flags=re.IGNORECASE)
    # Remove the word "unclassified"
    label = re.sub(r"\bunclassified\b", "", label, flags=re.IGNORECASE)
    # Remove any variant of "Same Domains" (e.g., "Same_Domains", "Same Domains", "Same-Domains", etc.)
    label = re.sub(r"[-_\s]*Same[-_\s]*Domains", "", label, flags=re.IGNORECASE)
    # Clean extra spaces and underscores from the beginning and end
    label = label.strip(" _")
    
    # Use re.findall para capturar todas las coincidencias del patrón:
    # Se espera que la parte que nos interesa tenga el formato: NombreDeEspecie ---C[digits]
    matches = re.findall(r"([A-Z][a-z]+(?:_[a-z]+)*)(?:[\s_-]+)(---C\d+)", label)
    if matches:
        # Tomamos la última coincidencia, que suele ser la correcta
        species, cluster = matches[-1]
        return f"{species}_{cluster}"
    else:
        return label
    
    # Now extract only the species and cluster marker.
    # This regex looks for a pattern where the species name is followed by a cluster marker (---C followed by digits)
    m = re.search(r"([A-Za-z0-9]+(?:[_\s][A-Za-z0-9\.]+)*)[\s_-]+(---C\d+)", label)
    if m:
        species = m.group(1).strip().replace(" ", "_")
        cluster = m.group(2).strip()
        return f"{species}_{cluster}"
    else:
        return label

def clean_taxa_in_tree(tree):
    """
    Traverse the tree and replace each clade name with a cleaned version
    (only the species name and cluster marker).
    """
    for clade in tree.find_clades():
        if clade.name:
            clade.name = clean_label(clade.name)
    return tree

def clean_newick_string(newick_str):
    """
    Cleans all labels in a Newick tree string by parsing the tree with Bio.Phylo,
    cleaning each clade's name using clean_taxa_in_tree, and writing the cleaned tree back to a string.
    """
    try:
        tree_obj = Phylo.read(io.StringIO(newick_str), "newick")
    except Exception as e:
        raise Exception(f"Error parsing Newick string: {e}")
    tree_obj = clean_taxa_in_tree(tree_obj)
    output_io = io.StringIO()
    Phylo.write(tree_obj, output_io, "newick")
    return output_io.getvalue()

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
    Run MEGA-CC to generate a phylogenetic tree, then clean the tree using Bio.Phylo.

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
    # Build the path to the MEGA-CC configuration file
    mega_config_input = f"{ProLink_path}/mega_configs/{tree_type}_{bootstrap_replications}.mao"
    logging.info(f"\n-- Generating phylogenetic tree with MEGA-CC")
    mega_cmd = ['megacc', '-a', mega_config_input, '-d', muscle_output, '-o', mega_output]
    logging.debug(f"Running MEGA-CC: {' '.join(mega_cmd)}")

    # Capture stdout and stderr to review MEGA-CC messages
    mega_run = subprocess.run(mega_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    logger.debug(f"MEGA-CC stdout: {mega_run.stdout}")
    logger.debug(f"MEGA-CC stderr: {mega_run.stderr}")
    if mega_run.returncode != 0:
        logger.error("ERROR: MEGA-CC failed")
        raise RuntimeError("MEGA-CC failed")
    
    # Wait for a few seconds to give time for the output file to be written
    time.sleep(5)
    
    # Verify that the output file exists before attempting to clean it
    if not os.path.exists(mega_output):
        # If the expected file with .nwk is not found, try the alternative with .mega
        alternative = mega_output.rsplit('.', 1)[0] + ".mega"
        if os.path.exists(alternative):
            logging.info(f"Using alternative output file: {alternative}")
            mega_output = alternative
        else:
            logger.error(f"ERROR: MEGA-CC did not produce the output file: {mega_output}")
            raise FileNotFoundError(f"Output file {mega_output} not found")
    

    # Read the Newick tree, clean taxa names using Bio.Phylo, and write the cleaned tree back.
    try:
        with open(mega_output, 'r') as f:
            newick = f.read()
        cleaned_newick = clean_newick_string(newick, protein_name='alkene_reductase')
        with open(mega_output, 'w') as f:
            f.write(cleaned_newick)
        logging.info(f"Cleaned Newick tree saved in '{mega_output}'")
    except Exception as e:
        logger.error(f"ERROR while cleaning the Newick file: {e}")
        raise
