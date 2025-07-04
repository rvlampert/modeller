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

To prepare the data for modeling, create a ```fasta.csv``` file in the ```root``` directory with the following structure:

| Column      | Description                                                                 |
|-------------|-----------------------------------------------------------------------------|
| ```gene```      | The gene name or identifier.                                               |
| ```variant```   | The variant type (e.g., wild or mutation).                                 |
| ```fasta```     | The FASTA sequence of the protein.                                         |

#### Example:
```csv
gene;variant;fasta
atpE;wild;MDPTIAAGALIGGGLIMAGGAIGAGIGDGVAGNALISGVARQPEAQGRLFTPFFITVGLVEAAYFINLAFMALFVFATPVK;
Rv0678;wild;VSVNDGVDQMGAEPDIMEFVEQMGGYFESRSLTRLAGRLLGWLLVCDPERQSSEELATALAASSGGISTNARMLIQFGFIERLAVAGDRRTYFRLRPNAFAAGERERIRAMAELQDLADVGLRALGDAPPQRSRRLREMRDLLAYMENVVSDALGRYSQRTGEDD;
ddn;ddn_Leu49Pro;MPKSPPRFLNSPLSDFFIKWMSRINTWMYRRNDGEGLGGTFQKIPVALPTTTGRKTGQPRVNPLYFLRDGGRVIVAASKGGAEKNPMWYLNLKANPKVQVQIKKEVLDLTARDATDEERAEYWPQLVTMYPSYQDYQSWTDRTIPIVVCEP;
```
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
docker run --name modeller -v modeller:/app modeller:1.0.0
```

---

## 4. Downloading Results

After running the modeling script, the generated **PDB files** will be available in the ```data/results``` directory. To download them, use the following command:

```sh
docker cp modeller:/app/assets/pdbs results
```

## 5. Cleaning Up
To clean up the Docker container and remove it after use, run:

```sh
docker rm -f modeller
docker volume rm modeller
```

## 6 Complete script
```sh
docker build --no-cache --build-arg KEY_MODELLER=<Modeller key> -t modeller:1.0.0 .
docker run --name modeller -v modeller:/app modeller:1.0.0
docker cp modeller:/app/assets/pdbs results
docker rm -f modeller
docker volume rm modeller
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