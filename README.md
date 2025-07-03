# Protein Modeling Project

---

## 1. Overview

This project is designed for **modeling protein structures** by generating **3D PDB files** from **FASTA sequences** using the **Modeller** tool. The workflow involves:

1. Using **BLAST** to search for protein templates.
2. Employing **Modeller** to construct 3D models based on these templates.

---

## 2. Setup

### 2.1 Prerequisites
- **Docker**: Ensure Docker is installed on your system.

---

### 2.2 Data Preparation

To prepare the data for modeling, create a ```fasta.csv``` file in the ```data``` directory with the following structure:

| Column      | Description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| ```gene```      | The gene name or identifier.                                               |
| ```variant```   | The variant type (e.g., wild or mutation).                                 |
| ```fasta```     | The FASTA sequence of the protein.                                         |
| ```modeller```  | *(Optional)* This column will be filled with ```OK``` after successful modeling.|

#### Example:
```csv
gene;variant;fasta;modeller
atpE;wild;MDPTIAAGALIGGGLIMAGGAIGAGIGDGVAGNALISGVARQPEAQGRLFTPFFITVGLVEAAYFINLAFMALFVFATPVK;
Rv0678;wild;VSVNDGVDQMGAEPDIMEFVEQMGGYFESRSLTRLAGRLLGWLLVCDPERQSSEELATALAASSGGISTNARMLIQFGFIERLAVAGDRRTYFRLRPNAFAAGERERIRAMAELQDLADVGLRALGDAPPQRSRRLREMRDLLAYMENVVSDALGRYSQRTGEDD;
ddn;ddn_Leu49Pro;MPKSPPRFLNSPLSDFFIKWMSRINTWMYRRNDGEGLGGTFQKIPVALPTTTGRKTGQPRVNPLYFLRDGGRVIVAASKGGAEKNPMWYLNLKANPKVQVQIKKEVLDLTARDATDEERAEYWPQLVTMYPSYQDYQSWTDRTIPIVVCEP;
```

#### Notes:
1. The ```modeller``` column is optional and will be automatically updated with ```OK``` after successful modeling.
2. If you need to add new instances to an existing dataset, simply append new rows to ```fasta.csv``` and move old results to ```data/assets/pdbs``` folder. The script will check and update the ```modeller``` column for existing files.

---

### 2.3 Installation

Build the Docker environment by running the following command:

```sh
docker build --no-cache --build-arg KEY_MODELLER=<Modeller key> -t modeller:1.0.0 .
```

Replace <```Modeller key```> with your Modeller license key.

---

## 3. Running the Project

Run the modeling script using Docker:

```sh
docker run --name modeller modeller:1.0.0
```

---

## 4. Downloading Results

After running the modeling script, the generated **PDB files** will be available in the ```data/results``` directory. To download them, use the following command:

```sh
docker cp modeller:/app/data/assets/pdbs ./data/results
```

---

## 5. Key Concepts

1. **FASTA**: A text format representing nucleotide or peptide sequences.
2. **BLAST**: Basic Local Alignment Search Tool used for finding regions of similarity between sequences.
3. **PDB**: Protein Data Bank files containing 3D coordinates of protein structures.
4. **Modeller**: A tool for comparative protein structure modeling.

---

## 6. Contact

For any questions or suggestions, feel free to reach out:  
Email: **[rvlampert@gmail.com](mailto:rvlampert@gmail.com)**

---

### Happy Modeling!