
import logging
import subprocess
import time
import os
import re

from .. import ProLink_path
from Bio import Phylo
import io

logger = logging.getLogger()

def clean_taxa_in_tree(tree, protein_name='alkene_reductase'):
    """
    Traverse the tree and clean each clade name using the clean_label function.
    """
    for clade in tree.find_clades():
        if clade.name:
            # Clean the clade name using your clean_label function
            clade.name = clean_label(clade.name, protein_name)
    return tree

def clean_label(label, protein_name='alkene_reductase'):
    """
    Cleans a single Newick label by removing:
      - WP codes (pattern: WP_\d{9}\.\d)
      - The term "MULTISPECIES:" if present
      - The protein name (default 'alkene_reductase')
      - The word "unclassified" if present
      - Variants of "Same Domains" (with hyphens, underscores, or spaces)
    Then, it extracts and returns only the species name and the cluster marker.
    
    For example, from:
      WP_072607337.1_alkene_reductase_Aquibium_oceanicum_---C51---Same_Domains
    it returns:
      Aquibium_oceanicum_---C51
    """
    # Remove any surrounding quotes
    label = label.strip("'\"")
    # Remove WP codes
    label = re.sub(r"WP[\s_]\d{9}\.\d", "", label)
    # Remove "MULTISPECIES:" if present
    label = re.sub(r"MULTISPECIES:\s*", "", label, flags=re.IGNORECASE)
    # Remove the protein name (allowing both underscores and spaces)
    protein_regex = re.escape(protein_name).replace(r'\_', r'[\s_]+')
    label = re.sub(protein_regex, "", label, flags=re.IGNORECASE)
    # Remove the word "unclassified"
    label = re.sub(r"\bunclassified\b", "", label, flags=re.IGNORECASE)
    # Remove variants of "Same Domains" (accepting hyphens, underscores, or spaces)
    label = re.sub(r"[-_\s]*Same[-_\s]*Domains", "", label, flags=re.IGNORECASE)
    # Clean extra spaces and underscores
    label = label.strip(" _")
    
    # Extract species and cluster marker (assuming the cluster marker starts with '---C' followed by digits)
    pattern = re.compile(
        r"(?P<species>[A-Za-z0-9]+(?:[_\s][A-Za-z0-9\.]+)*)"  
        r"[\s_-]+(?P<cluster>---C\d+)",
        flags=re.IGNORECASE
    )
    match = pattern.search(label)
    if match:
        species = match.group('species').strip().replace(" ", "_")
        cluster = match.group('cluster').strip()
        return f"{species}_{cluster}"
    else:
        return label

def clean_newick_string(newick_str, protein_name='alkene_reductase'):
    """
    Cleans all labels in a Newick tree string by applying the clean_label function.
    It looks for labels that may include an optional branch length (e.g., ":0.13368245").
    The branch length is separated, the label cleaned, and then the branch length reattached.
    """
    # Pattern to match labels (quoted or unquoted) that include a cluster marker,
    # optionally followed by a branch length (starting with a colon)
    pattern = re.compile(
        r"('([^':]+---C\d+[^':]*)'|\"([^\":]+---C\d+[^\":]*)\"|([A-Za-z0-9 _\.\-]+---C\d+))(:[0-9.eE+-]+)?",
        flags=re.IGNORECASE
    )
    def replacer(match):
        full_match = match.group(0)
        branch = match.group(5) if match.group(5) is not None else ""
        if branch:
            label_part = full_match[:-len(branch)]
        else:
            label_part = full_match
        cleaned = clean_label(label_part, protein_name)
        return f"{cleaned}{branch}"
    return pattern.sub(replacer, newick_str)

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
    

    # Read the generated Newick tree, parse it with Bio.Phylo, clean taxa names, and write it back out
    try:
        with open(mega_output, 'r') as f:
            newick = f.read()
        # Parse the tree from the Newick string
        tree_obj = Phylo.read(io.StringIO(newick), "newick")
        # Clean each clade's name
        tree_obj = clean_taxa_in_tree(tree_obj, protein_name='alkene_reductase')
        # Write the cleaned tree to a string
        output_io = io.StringIO()
        Phylo.write(tree_obj, output_io, "newick")
        cleaned_newick = output_io.getvalue()
        # Write the cleaned Newick tree back to the file
        with open(mega_output, 'w') as f:
            f.write(cleaned_newick)
        logging.info(f"Cleaned Newick tree saved in '{mega_output}'")
    except Exception as e:
        logger.error(f"ERROR while cleaning the Newick file: {e}")
        raise

