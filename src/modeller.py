import glob
import shutil
import pandas as pd
import numpy as np
import os
import re
import requests
import logging
from Bio.Blast import NCBIWWW  

from modeller import *  
from modeller.automodel import *  

logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    handlers=[logging.StreamHandler()]
)

ASSETS_PATH = f"assets/"
FASTA_PATH = f"{ASSETS_PATH}fasta/"
PDB_PATH = f"{ASSETS_PATH}pdbs/"
BLAST_PATH = f"{ASSETS_PATH}blast/"
PIR_PATH = f"{ASSETS_PATH}pir/"
ALI_PATH = f"{ASSETS_PATH}ali/"
FASTA_FILE = f"fastas.csv"  
MODEL = 'modeller'

env = Environ()  

def generate_blast(gene, fasta):
    logging.info(f"Generating BLAST for gene: {gene}")
    logging.debug(f"FASTA content for {gene}: {fasta[:500]}...")  # Log first 500 characters of the fasta
    result_handle = NCBIWWW.qblast("blastp", "pdb", fasta)
    blast_result = result_handle.read()  
    logging.debug(f"BLAST result for {gene}: {blast_result[:500]}...")  # Log first 500 characters of the result
    result_handle.close() 
    pdb_ids = re.findall(r'<Hit_id>pdb\|(\w+)\|', blast_result)
    logging.info(f"Found PDB IDs: {pdb_ids}")
    download_successful = False
    for pdb_id in pdb_ids: 
        pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb" 
        logging.info(f"Attempting to download {pdb_id}.pdb for gene {gene} from URL: {pdb_url}")
        response = requests.get(pdb_url) 
        if response.status_code == 200:
            blast_content = response.content
            os.makedirs(BLAST_PATH, exist_ok=True)  # Ensure BLAST_PATH exists
            blast_filename = f"{BLAST_PATH}{gene}.pdb"
            with open(blast_filename, "wb") as file: 
                file.write(blast_content)  
            logging.info(f"Successfully downloaded {gene}.pdb")
            download_successful = True 
            break  
        else:
            logging.error(f"Failed to download {pdb_id}.pdb for gene {gene}")  
    if not download_successful: 
        logging.error(f"Failed to download any PDB file for gene {gene}")
    logging.info(f"Successfully generated BLAST file for gene {gene}: {blast_filename}")
    return blast_filename

def generate_pir(id, fasta, blast_filename):
    logging.info(f"Generating PIR file for {id}")
    if os.path.isfile(blast_filename):
        ali_content = f""">P1;{id}
sequence:{id}:::::::0.00: 0.00
{fasta}*
"""
        os.makedirs(PIR_PATH, exist_ok=True)  # Ensure PIR_PATH exists
        pir_filename = f"{PIR_PATH}{id}.txt"
        with open(pir_filename, 'w') as file:
            file.write(ali_content)
        logging.info(f"Successfully generated PIR file: {pir_filename}")  
    else:
        logging.error(f"PDB file for {id} not found.")
    return pir_filename

def generate_ali(id, gene, pir_filename):
    logging.info(f"Generating alignment file for {id}")
    ali_filename = pir_filename.replace('.txt', '.ali').replace('/pir/', '/ali/')  
    os.makedirs(ALI_PATH, exist_ok=True) 
    if not os.path.exists(ali_filename):
        aln = Alignment(env)  
        md1 = Model(env, file=f'{BLAST_PATH}{gene}.pdb')
        aln.append_model(md1, align_codes=gene, atom_files=f'{BLAST_PATH}{gene}.pdb')
        aln.append(file=f'{PIR_PATH}{id}.txt', align_codes=id)
        aln.align2d()
        aln.write(file=ali_filename, alignment_format='PIR')
        logging.info(f"Successfully generated alignment file: {ali_filename}")  
    else:
        logging.info(f"Alignment file for {id} already exists. Skipping.")
    return ali_filename

def generate_pdb(id, gene, ali_filename):
    logging.info(f"Generating PDB file for {id}")
    if os.path.exists(ali_filename):  
        env = Environ() 
        env.io.atom_files_directory = ['.', BLAST_PATH]  
        a = AutoModel(env,
                      alnfile=ali_filename, 
                      knowns=gene,  
                      sequence=id)  
        a.starting_model = 1 
        a.ending_model = 1                      
        a.make()  
        logging.info(f"Successfully generated PDB file for {id}")
        move_pdb_file(id)  
        delete_other_files(id)  
    else:
        logging.error(f"ALI file for {id} not found.")

def move_pdb_file(variant):
    logging.info(f"Moving PDB file for variant: {variant}")
    pdb_files = glob.glob(os.path.join(f"{variant}*.pdb"))
    if not pdb_files:
        logging.warning(f"No PDB files found for variant {variant}")  
        return
    os.makedirs(PDB_PATH, exist_ok=True)  # Ensure PDB_PATH exists
    file_to_move = pdb_files[0]  
    destination_file = f"{PDB_PATH}{variant}.pdb" 
    shutil.move(file_to_move, destination_file) 
    logging.info(f"Moved file: {file_to_move} to {destination_file}") 

def delete_other_files(variant):
    logging.info(f"Deleting other files for variant: {variant}")
    files_to_delete = glob.glob(os.path.join(f"{variant}*.*"))  
    for file in files_to_delete:
        if not file.endswith(".pdb"):  
            os.remove(file)  
            logging.info(f"Deleted file: {file}")

def get_pdb(gene, variant, fasta):
    logging.info(f"Starting PDB generation for gene: {gene}, variant: {variant}")
    blast_filename = generate_blast(gene, fasta)
    id = gene if variant == 'wild' else variant
    pir_filename = generate_pir(id, fasta, blast_filename) 
    ali_filename = generate_ali(id, gene, pir_filename)
    logging.info(f"Generated alignment file: {ali_filename}")
    generate_pdb(id, gene, ali_filename)

def get_pdbs():
    df = pd.read_csv(FASTA_FILE, sep=';')
    logging.info("Starting PDB generation for all entries in the DataFrame.")
    if not os.path.exists(PDB_PATH):
        os.makedirs(PDB_PATH) 
        logging.info(f"Created PDB directory: {PDB_PATH}")
    if not df.empty:
        for _, row in df.iterrows():
            logging.info("Start modeling...")
            get_pdb(row["gene"], row["variant"], row["fasta"])
            logging.info(f"Successfully modeled, you can find the generated PDBs in: {PDB_PATH}")
    else:
        logging.warning("No models to be created.")