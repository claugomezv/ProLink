
import logging
import subprocess
import time
import os
import re

from .. import ProLink_path

logger = logging.getLogger()

def clean_label(label, protein_name='alkene_reductase'):
    """
    Cleans a single Newick label by removing:
      - WP codes (pattern: WP_\d{9}\.\d)
      - The term "MULTISPECIES:" if present
      - The protein name (default 'alkene_reductase')
      - Trailing "Same Domains" (or similar variants)
    Then, it extracts and returns only the species name and the cluster marker.
    
    For example, from:
      WP_058328214.1_alkene_reductase_Sinorhizobium_sp._Sb3_---C28---Same_Domains
    or
      WP 062476070.1 MULTISPECIES: alkene reductase unclassified Rhizobium ---C22---Same Domains
    it returns:
      Sinorhizobium_sp._Sb3_---C28
    """
    # Remove any surrounding quotes
    label = label.strip("'\"")
    # Remove WP codes
    label = re.sub(r"WP[\s_]\d{9}\.\d", "", label)
    # Remove "MULTISPECIES:" if present
    label = re.sub(r"MULTISPECIES:\s*", "", label, flags=re.IGNORECASE)
    # Remove the protein name
    # Allow both underscore and space in the protein name (e.g., "alkene_reductase" or "alkene reductase")
    protein_regex = re.escape(protein_name).replace(r'\_', r'[\s_]+')
    label = re.sub(protein_regex, "", label, flags=re.IGNORECASE)
    # Remove the word "unclassified"
    label = re.sub(r"\bunclassified\b", "", label, flags=re.IGNORECASE)
    # Remove variants of "Same Domains" (with hyphens, underscores, or spaces)
    label = re.sub(r"[-_]*Same[_\s]*Domains", "", label, flags=re.IGNORECASE)
    # Clean extra spaces and underscores from the beginning and end
    label = label.strip(" _")
    
    # Use a regex to extract the species name and the cluster marker.
    # We assume that the cluster marker starts with '---C' followed by digits.
    pattern = re.compile(
        r"(?P<species>[A-Za-z0-9]+(?:[_\s][A-Za-z0-9\.]+)*)"  # species name: letters/numbers with underscores or spaces
        r"[\s_-]+(?P<cluster>---C\d+)",                       # cluster marker: ---C followed by digits
        flags=re.IGNORECASE
    )
    match = pattern.search(label)
    if match:
        species = match.group('species').strip().replace(" ", "_")
        cluster = match.group('cluster').strip()
        return f"{species}_{cluster}"
    else:
        # If the pattern is not found, return the cleaned label as is
        return label

def clean_newick_string(newick_str, protein_name='alkene_reductase'):
    """
    Cleans all labels in a Newick tree string by applying the clean_label function.
    It looks for any substring that contains a cluster marker (---C followed by digits)
    and replaces it with the cleaned version.
    """
    # This regex matches labels (quoted or unquoted) that include a cluster marker,
    # followed by an optional branch length (starting with a colon).
    pattern = re.compile(
        r"('([^':]+---C\d+[^':]*)'|\"([^\":]+---C\d+[^\":]*)\"|([A-Za-z0-9 _\.\-]+---C\d+))(:[0-9.eE+-]+)?",
        flags=re.IGNORECASE
    )
    def replacer(match):
        full_label = match.group(0)
        branch = match.group(5) if match.group(5) is not None else ""
        # Remove branch length from the label portion
        if branch:
            label_part = full_label[:-len(branch)]
        else:
            label_part = full_label
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
    

    # Read the generated Newick tree and clean its labels
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
