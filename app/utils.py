import subprocess
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

def run_pepsmi(sequence):
    try:
        result = subprocess.run(["pepsmi", sequence], capture_output=True, text=True, check=True)
        smiles = result.stdout.strip()
        return smiles
    except subprocess.CalledProcessError as e:
        raise ValueError(f"PepSMI conversion failed: {e.stderr}")

def run_alphafold_convert(input_cif, output_pdb):
    try:
        subprocess.run(["alphafold_convert", input_cif, output_pdb], check=True)
        return output_pdb
    except subprocess.CalledProcessError as e:
        raise ValueError(f"AlphaFold conversion failed: {e.stderr}")

def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        return True, None
    else:
        return False, "Invalid SMILES string. Please check the syntax and try again."

def validate_pdb_file(pdb_file):
    try:
        content = pdb_file.read().decode('utf-8')
        for line in content.splitlines():
            if not line.startswith(("ATOM", "HETATM", "END", "HEADER")):
                return False, "Invalid PDB format."
        return True, None
    except Exception as e:
        return False, f"Error reading PDB file: {str(e)}"

def convert_pdb_to_pdbqt(input_pdb, output_pdbqt):
    try:
        subprocess.run(["prepare_receptor4.py", "-r", input_pdb, "-o", output_pdbqt], check=True)
        return output_pdbqt
    except subprocess.CalledProcessError as e:
        raise ValueError(f"PDB to PDBQT conversion failed: {e.stderr}")

def convert_smiles_to_pdbqt(smiles, output_pdbqt):
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string.")
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.UFFOptimizeMolecule(mol)
    Chem.MolToPDBFile(mol, "temp.pdb")
    try:
        subprocess.run(["prepare_ligand4.py", "-l", "temp.pdb", "-o", output_pdbqt], check=True)
        return output_pdbqt
    except subprocess.CalledProcessError as e:
        raise ValueError(f"SMILES to PDBQT conversion failed: {e.stderr}")

def parse_vina_output(output_file):
    scores = []
    with open(output_file, 'r') as f:
        for line in f:
            if "REMARK VINA RESULT" in line:
                parts = line.strip().split()
                score = parts[-1]
                scores.append(score)
    return scores

def calculate_tanimoto_similarity(input_smiles, known_drugs):
    input_mol = Chem.MolFromSmiles(input_smiles)
    if not input_mol:
        raise ValueError("Invalid input SMILES string.")
    input_fp = AllChem.GetMorganFingerprintAsBitVect(input_mol, 2)
    results = []
    for drug in known_drugs:
        drug_mol = Chem.MolFromSmiles(drug.smiles)
        if drug_mol:
            drug_fp = AllChem.GetMorganFingerprintAsBitVect(drug_mol, 2)
            similarity = DataStructs.TanimotoSimilarity(input_fp, drug_fp)
            results.append({"name": drug.name, "smiles": drug.smiles, "score": similarity, "description": drug.description})
    return sorted(results, key=lambda x: x["score"], reverse=True)
