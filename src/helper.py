import os
import pandas as pd
import numpy as np

from src.constants import FASTA_FILE, PDB_PATH

def update_column_based_on_files(df, column):
    for index, row in df.iterrows():
        file_name = f"{row['variant']}.pdb"
        current_dir = os.getcwd()
        file_path = f"{current_dir}/{PDB_PATH}{column}/{file_name}"
        if os.path.exists(file_path):
            df.at[index, column] = 'ok'
        else:
            df.at[index, column] = np.nan
    return df

def update_csv_file(df, column):  
    df = update_column_based_on_files(df, column)
    original_df = pd.read_csv(FASTA_FILE, sep=';')
    for index, row in df.iterrows():
        if index in original_df.index:
            original_df.loc[index, column] = row[column]
    original_df.to_csv(FASTA_FILE, sep=';', index=False)

def write_file(content, filepath):
    """
    Writes content to a file.
    """
    try:
        with open(filepath, 'w') as file:
            file.write(content)
        return filepath
    except Exception as e:
        return f"Error writing file {filepath}: {str(e)}"
    
def move_pdb_file(variant):
    """
    Moves the generated PDB file to the PDB_PATH directory.
    """
    try:
        pdb_files = glob.glob(f"{variant}*.pdb")
        if not pdb_files:
            return f"Error: No PDB files found for variant {variant}"
        shutil.move(pdb_files[0], os.path.join(PDB_PATH, f"{variant}.pdb"))
    except Exception as e:
        return f"Error in move_pdb_file: {str(e)}"

def delete_other_files(variant):
    """
    Deletes all files related to the variant except the PDB file.
    """
    try:
        files_to_delete = glob.glob(f"{variant}*.*")
        for file in files_to_delete:
            if not file.endswith(".pdb"):
                os.remove(file)
    except Exception as e:
        return f"Error in delete_other_files: {str(e)}"