import os
import re
import requests
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW
from pandarallel import pandarallel
from tqdm import tqdm
from modeller import *  
from modeller.automodel import *  
from src.helper import delete_other_files, move_pdb_file, update_column_based_on_files, update_csv_file, write_file
from src.constants import FASTA_FILE, PDB_PATH, BLAST_PATH, PIR_PATH

pandarallel.initialize(progress_bar=True)

MODEL = 'modeller'
env = Environ()

def download_pdb_file(pdb_id, gene):
    """
    Downloads a PDB file from the RCSB PDB database.
    """
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    response = requests.get(pdb_url)
    if response.status_code == 200:
        blast_filename = os.path.join(BLAST_PATH, f"{gene}.pdb")
        with open(blast_filename, "wb") as file:
            file.write(response.content)
        return blast_filename
    return None

def generate_blast(gene, fasta):
    """
    Runs a BLAST search and downloads the corresponding PDB file.
    """
    try:
        result_handle = NCBIWWW.qblast("blastp", "pdb", fasta)
        blast_result = result_handle.read()
        result_handle.close()
        pdb_ids = re.findall(r'<Hit_id>pdb\|(\w+)\|', blast_result)
        for pdb_id in pdb_ids:
            blast_filename = download_pdb_file(pdb_id, gene)
            if blast_filename:
                return blast_filename
        return None
    except Exception as e:
        return f"Error in generate_blast: {str(e)}"

def generate_pir(id, fasta, blast_filename):
    """
    Generates a PIR file for the given ID and FASTA sequence.
    """
    if not os.path.isfile(blast_filename):
        return None
    ali_content = f""">P1;{id}
sequence:{id}:::::::0.00: 0.00
{fasta}*
"""
    pir_filename = os.path.join(PIR_PATH, f"{id}.txt")
    return write_file(ali_content, pir_filename)

def generate_ali(id, gene, pir_filename):
    """
    Generates an ALI file using the PIR file and BLAST results.
    """
    try:
        ali_filename = pir_filename.replace('.txt', '.ali').replace('/pir/', '/ali/')
        if not os.path.exists(ali_filename):
            aln = Alignment(env)
            md1 = Model(env, file=os.path.join(BLAST_PATH, f"{id}.pdb"))
            aln.append_model(md1, align_codes=gene, atom_files=os.path.join(BLAST_PATH, f"{gene}.pdb"))
            aln.append(file=pir_filename, align_codes=id)
            aln.align2d()
            aln.write(file=ali_filename, alignment_format='PIR')
        return ali_filename
    except Exception as e:
        return f"Error in generate_ali: {str(e)}"

def generate_pdb(id, gene, ali_filename):
    """
    Generates a PDB file using the ALI file.
    """
    try:
        if not os.path.exists(ali_filename):
            return f"Error: ALI file for {id} not found."
        env.io.atom_files_directory = ['.', BLAST_PATH]
        a = AutoModel(env, alnfile=ali_filename, knowns=gene, sequence=id)
        a.starting_model = 1
        a.ending_model = 1
        a.make()
        move_pdb_file(id)
        delete_other_files(id)
    except Exception as e:
        return f"Error in generate_pdb: {str(e)}"

def process_row(row):
    """
    Processes a single row to generate a PDB file for the given gene and variant.
    """
    try:
        gene = row["gene"]
        variant = row["variant"]
        fasta = row["fasta"]
        blast_filename = generate_blast(gene, fasta)
        if not blast_filename:
            return f"Error: Failed to generate BLAST file for gene {gene}"
        id = gene if variant == 'wild' else variant
        pir_filename = generate_pir(id, fasta, blast_filename)
        if not pir_filename:
            return f"Error: Failed to generate PIR file for {id}"
        ali_filename = generate_ali(id, gene, pir_filename)
        if not ali_filename:
            return f"Error: Failed to generate ALI file for {id}"
        pdb_error = generate_pdb(id, gene, ali_filename)
        if pdb_error:
            return pdb_error
        return "ok"
    except Exception as e:
        return f"Error in process_row: {str(e)}"

def get_pdbs():
    """
    Processes all genes and variants in the FASTA file to generate PDB files.
    """
    try:
        df = pd.read_csv(FASTA_FILE, sep=';')
        os.makedirs(PDB_PATH, exist_ok=True)
        if MODEL not in df.columns:
            df[MODEL] = np.nan
        df = update_column_based_on_files(df, MODEL)
        todo_df = df[df[MODEL] != 'ok']
        if todo_df.empty:
            return "No models to be created." 
        todo_df[MODEL] = todo_df.parallel_apply(process_row, axis=1)
        update_csv_file(df, MODEL)
    except Exception as e:
        return f"Error in get_pdbs: {str(e)}"