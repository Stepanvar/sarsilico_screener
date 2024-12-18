import subprocess
from rdkit import Chem

def validate_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None

def validate_pdb_file(pdb_file):
    content = pdb_file.read().decode('utf-8')
    for line in content.splitlines():
        if not line.startswith(("ATOM", "HETATM", "END")):
            return False, "Invalid PDB format."
    return True, None

def run_pepsmi(sequence):
    result = subprocess.run(["pepsmi", sequence], capture_output=True, text=True)
    if result.returncode == 0:
        return result.stdout.strip()
    else:
        raise ValueError("PepSMI conversion failed.")

def run_alphafold_convert(input_cif, output_pdb):
    subprocess.run(["alphafold_convert", input_cif, output_pdb], check=True)

def convert_to_pdbqt(input_pdb, output_pdbqt):
    subprocess.run(["prepare_receptor4.py", "-r", input_pdb, "-o", output_pdbqt], check=True)

def run_vina(receptor_pdbqt, ligand_pdbqt, output, config):
    command = [
        "vina",
        "--receptor", receptor_pdbqt,
        "--ligand", ligand_pdbqt,
        "--out", output,
        "--config", config,
    ]
    subprocess.run(command, check=True)

def parse_vina_output(output_file):
    scores = []
    with open(output_file, 'r') as f:
        for line in f:
            if "REMARK VINA RESULT" in line:
                parts = line.strip().split()
                # Example: REMARK VINA RESULT:     -7.5
                # Usually the score is after "RESULT:"
                # Adjust parsing as per actual output format
                # Assuming the score is parts[-1]
                score = parts[-1]
                scores.append(score)
    return scores

def calculate_tanimoto_similarity(input_smiles, known_drugs):
    from rdkit.Chem import AllChem, DataStructs
    input_mol = Chem.MolFromSmiles(input_smiles)
    input_fp = AllChem.GetMorganFingerprintAsBitVect(input_mol, 2)
    results = []
    for drug in known_drugs:
        drug_mol = Chem.MolFromSmiles(drug.smiles)
        if drug_mol:
            drug_fp = AllChem.GetMorganFingerprintAsBitVect(drug_mol, 2)
            similarity = DataStructs.TanimotoSimilarity(input_fp, drug_fp)
            results.append({"name": drug.name, "score": similarity})
    return sorted(results, key=lambda x: x["score"], reverse=True)
