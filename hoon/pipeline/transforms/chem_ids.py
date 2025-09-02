def derive_chem_ids(smiles):
    """Args: smiles(str|None) -> dict

    SMILES에서 캐노니컬 SMILES, InChI, InChIKey를 생성한다.
    구조를 만들 수 없으면 has_structure=False로 반환.
    """
    if smiles is None:
        return {"has_structure": False}
    s = str(smiles).strip()
    if s == "":
        return {"has_structure": False}
    try:
        from rdkit import Chem
        from rdkit.Chem import inchi
    except Exception:
        return {"has_structure": False}
    try:
        mol = Chem.MolFromSmiles(s)
        if mol is None:
            return {"has_structure": False}
        smiles_canonical = Chem.MolToSmiles(mol, canonical=True)
        ik = ""
        try:
            ik = inchi.MolToInchiKey(mol)
        except Exception:
            ik = ""
        out = {
            "has_structure": True,
            "smiles_canonical": smiles_canonical,
            "compound_key": ik,
            "inchikey": ik,
            "inchikey14": ik[:14] if ik else "",
        }
        return out
    except Exception:
        return {"has_structure": False}
