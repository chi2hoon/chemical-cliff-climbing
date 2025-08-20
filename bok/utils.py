# File: activity_cliff_app/utils.py
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdFMCS import FindMCS
import numpy as np
import difflib
from typing import Tuple

def compute_pIC50(ic50_nM):
    return -np.log10(ic50_nM * 1e-9)

def smiles_to_image(smiles, size=(300, 300)):
    mol = Chem.MolFromSmiles(smiles)
    return Draw.MolToImage(mol, size=size)

def smiles_diff_to_images(smiles1, smiles2, size=(300, 300)):
    mol1 = Chem.MolFromSmiles(smiles1)
    mol2 = Chem.MolFromSmiles(smiles2)

    if mol1 is None or mol2 is None:
        return None, None # Handle invalid SMILES

    # Find Maximum Common Substructure (MCS)
    mcs_result = FindMCS([mol1, mol2])
    mcs_smarts = mcs_result.smartsString
    
    highlight_atoms1 = []
    highlight_atoms2 = []

    if not mcs_smarts:
        # If no MCS found, highlight all atoms
        highlight_atoms1 = list(range(mol1.GetNumAtoms()))
        highlight_atoms2 = list(range(mol2.GetNumAtoms()))
    else:
        mcs_mol = Chem.MolFromSmarts(mcs_smarts)
        
        # Get atom indices in mol1 that are part of MCS
        mol1_mcs_matches = mol1.GetSubstructMatches(mcs_mol)
        mol1_mcs_atoms_indices = set()
        for match in mol1_mcs_matches:
            mol1_mcs_atoms_indices.update(match)

        # Get atom indices in mol2 that are part of MCS
        mol2_mcs_matches = mol2.GetSubstructMatches(mcs_mol)
        mol2_mcs_atoms_indices = set()
        for match in mol2_mcs_matches:
            mol2_mcs_atoms_indices.update(match)

        # Identify differing atoms (not part of MCS)
        highlight_atoms1 = [atom.GetIdx() for atom in mol1.GetAtoms() if atom.GetIdx() not in mol1_mcs_atoms_indices]
        highlight_atoms2 = [atom.GetIdx() for atom in mol2.GetAtoms() if atom.GetIdx() not in mol2_mcs_atoms_indices]

    # Generate images with highlighting
    img1 = Draw.MolToImage(mol1, size=size, highlightAtoms=highlight_atoms1, highlightColor=(1, 0, 0))
    img2 = Draw.MolToImage(mol2, size=size, highlightAtoms=highlight_atoms2, highlightColor=(1, 0, 0))

    return img1, img2


def _html_escape(text: str) -> str:
    return (
        text.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
    )


def highlight_canonical_smiles_diff(a: str, b: str) -> Tuple[str, str]:
    """Produce HTML with colored differences for two canonical SMILES strings.

    Returns two HTML strings (for a and b) with differing segments highlighted.
    """
    if a is None:
        a = ""
    if b is None:
        b = ""

    matcher = difflib.SequenceMatcher(a=a, b=b)

    a_out_parts = []
    b_out_parts = []

    for tag, i1, i2, j1, j2 in matcher.get_opcodes():
        a_seg = _html_escape(a[i1:i2])
        b_seg = _html_escape(b[j1:j2])
        if tag == 'equal':
            a_out_parts.append(a_seg)
            b_out_parts.append(b_seg)
        elif tag == 'replace':
            a_out_parts.append(f"<span style=\"background:#ffdddd;color:#a40000\">{a_seg}</span>")
            b_out_parts.append(f"<span style=\"background:#ddeeff;color:#003380\">{b_seg}</span>")
        elif tag == 'delete':
            a_out_parts.append(f"<span style=\"background:#ffdddd;color:#a40000\">{a_seg}</span>")
        elif tag == 'insert':
            b_out_parts.append(f"<span style=\"background:#ddeeff;color:#003380\">{b_seg}</span>")

    a_html = "".join(a_out_parts)
    b_html = "".join(b_out_parts)

    # Monospace formatting and preserved spacing
    a_html = f"<div style=\"font-family:monospace; white-space:pre-wrap;\">{a_html}</div>"
    b_html = f"<div style=\"font-family:monospace; white-space:pre-wrap;\">{b_html}</div>"

    return a_html, b_html