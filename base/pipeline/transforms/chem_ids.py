def derive_chem_ids(smiles, iupac_name=None):
    """Args: smiles(str|None), iupac_name(str|None) -> dict

    SMILES에서 캐노니컬 SMILES, InChIKey를 생성한다.
    - RDKit가 nitro 표기(N(O)=O, N(=O)O 등)로 실패하는 경우를 보정해 재시도한다.
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

    def _fix_common_smiles_issues(text):
        """Args: text(str) -> str

        자주 발생하는 nitro 그룹 표기를 정전하 표기로 교체한다.
        예: N(O)=O, N(=O)O, N(=O)(O) → [N+](=O)[O-]
        """
        import re as _re
        pat = _re.compile(r"N\(=O\)=O|N\(O\)=O|N\(=O\)O|N\(=O\)\(O\)")
        return pat.sub("[N+](=O)[O-]", text)

    def _mk(mol):
        smiles_canonical = Chem.MolToSmiles(mol, canonical=True)
        ik = ""
        try:
            ik = inchi.MolToInchiKey(mol)
        except Exception:
            ik = ""
        if not ik:
            try:
                import hashlib as _h
                ik = "SMI:" + _h.sha1(smiles_canonical.encode("utf-8")).hexdigest()
            except Exception:
                ik = ""
        return {
            "has_structure": True,
            "smiles_canonical": smiles_canonical,
            "compound_key": ik,
            "inchikey": ik,
            "inchikey14": (ik[:14] if ik and not ik.startswith("SMI:") else ""),
        }

    try:
        mol = Chem.MolFromSmiles(s)
        if mol is not None:
            return _mk(mol)
    except Exception:
        mol = None
    s2 = _fix_common_smiles_issues(s)
    if s2 != s:
        try:
            mol2 = Chem.MolFromSmiles(s2)
            if mol2 is not None:
                return _mk(mol2)
        except Exception:
            pass
    if iupac_name:
        smi_from_name = None
        try:
            from opsin import opsin as _opsin
            try:
                smi_from_name = _opsin.iupac_to_smiles(iupac_name)
            except Exception:
                smi_from_name = None
        except Exception:
            try:
                import opsin as _ops
                smi_from_name = _ops.parseIUPACName(iupac_name)
            except Exception:
                smi_from_name = None
        if smi_from_name:
            smi_from_name = _fix_common_smiles_issues(str(smi_from_name).strip())
            try:
                mol3 = Chem.MolFromSmiles(smi_from_name)
                if mol3 is not None:
                    return _mk(mol3)
            except Exception:
                pass
    return {"has_structure": False}

