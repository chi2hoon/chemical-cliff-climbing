# File: activity_cliff_app/structure_diff.py
from rdkit import Chem
from rdkit.Chem import Descriptors

def detect_diff_type(mol1, mol2):
    chiral1 = Chem.FindMolChiralCenters(mol1, includeUnassigned=True)
    chiral2 = Chem.FindMolChiralCenters(mol2, includeUnassigned=True)
    if chiral1 != chiral2:
        return "stereochemistry"

    if "n" in Chem.MolToSmiles(mol2) and "c1ccccc1" in Chem.MolToSmiles(mol1):
        return "electronic effect"

    if abs(Descriptors.MolWt(mol1) - Descriptors.MolWt(mol2)) > 50:
        return "lipophilicity"

    return "unknown"
