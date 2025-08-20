# File: activity_cliff_app/cliff_detector.py
# RDKit을 기반으로 분자 유사도(Tanimoto)와 생물학적 활성(IC50 또는 pIC50 등)의 차이를 이용하여 Activity Cliff 후보 쌍을 자동으로 탐지
import streamlit as st # Import streamlit
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from itertools import combinations
import pandas as pd

# Correct import for GetMorganGenerator
from rdkit.Chem import rdFingerprintGenerator
from functools import lru_cache
import math

# Initialize the MorganGenerator once
morgan_gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def fingerprint(mol):
    # Use the MorganGenerator to get the fingerprint
    return morgan_gen.GetFingerprint(mol)


def to_canonical_smiles(mol: Chem.Mol) -> str:
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None

def similarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

@st.cache_data # Add the cache decorator
def detect_activity_cliffs(df, sim_thres=0.85, act_thres=1.0):
    df['mol'] = df['SMILES'].apply(Chem.MolFromSmiles)
    df['fp'] = df['mol'].apply(fingerprint)
    df['canonical_smiles'] = df['mol'].apply(to_canonical_smiles)

    results = []
    for i, j in combinations(df.index, 2):
        # skip identical canonical molecules
        if df.loc[i, 'canonical_smiles'] is not None and df.loc[i, 'canonical_smiles'] == df.loc[j, 'canonical_smiles']:
            continue
        sim = similarity(df.loc[i, 'fp'], df.loc[j, 'fp'])
        if sim >= sim_thres:
            diff = abs(df.loc[i, 'pIC50'] - df.loc[j, 'pIC50'])
            if diff >= act_thres:
                results.append({
                    'mol1_idx': i,
                    'mol2_idx': j,
                    'mol1_smiles': df.loc[i, 'SMILES'],
                    'mol2_smiles': df.loc[j, 'SMILES'],
                    'mol1_canonical': df.loc[i, 'canonical_smiles'],
                    'mol2_canonical': df.loc[j, 'canonical_smiles'],
                    'mol1_activity': df.loc[i, 'pIC50'],
                    'mol2_activity': df.loc[j, 'pIC50'],
                    'sim': sim,
                    'activity_diff': diff,
                })
    return pd.DataFrame(results)


# --- Alternative detection strategies ---

@lru_cache(maxsize=8)
def get_morgan_generator(radius: int = 2, fp_size: int = 2048):
    return rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)


@st.cache_data
def detect_activity_cliffs_sali(
    df: pd.DataFrame,
    top_n: int = 200,
    radius: int = 2,
    fp_bits: int = 2048,
    sim_floor: float = 0.5,
    scaffold_constrained: bool = False,
):
    """Rank pairs by SALI = delta_pIC50 / (1 - similarity).

    Returns a DataFrame with the same schema as the default method, plus 'sali'.
    """
    # Prepare molecules and fingerprints
    mols = df['SMILES'].apply(Chem.MolFromSmiles)
    gen = get_morgan_generator(radius=radius, fp_size=fp_bits)
    fps = mols.apply(gen.GetFingerprint)
    canon = mols.apply(to_canonical_smiles)

    # Optional scaffold strings for constraint
    scaffold_smiles = None
    if scaffold_constrained:
        scaffold_smiles = mols.apply(lambda m: MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m is not None else None)

    rows = []
    for i in df.index:
        fpi = fps[i]
        for j in df.index:
            if j <= i:
                continue
            # Skip identical canonical molecules
            if canon[i] is not None and canon[i] == canon[j]:
                continue
            # Scaffold filter
            if scaffold_constrained and scaffold_smiles[i] != scaffold_smiles[j]:
                continue

            sim = DataStructs.TanimotoSimilarity(fpi, fps[j])
            if sim < sim_floor:
                continue
            dp = abs(df.at[i, 'pIC50'] - df.at[j, 'pIC50'])
            denom = max(1e-6, 1.0 - sim)
            sali = dp / denom
            rows.append((i, j, sim, dp, sali))

    if not rows:
        return pd.DataFrame(columns=[
            'mol1_idx','mol2_idx','mol1_smiles','mol2_smiles','mol1_activity','mol2_activity','sim','activity_diff','sali'
        ])

    res = pd.DataFrame(rows, columns=['i', 'j', 'sim', 'dp', 'sali'])
    res = res.sort_values('sali', ascending=False).head(top_n)
    out = pd.DataFrame({
        'mol1_idx': res['i'],
        'mol2_idx': res['j'],
        'mol1_smiles': res['i'].map(df['SMILES']),
        'mol2_smiles': res['j'].map(df['SMILES']),
        'mol1_canonical': res['i'].map(canon),
        'mol2_canonical': res['j'].map(canon),
        'mol1_activity': res['i'].map(df['pIC50']),
        'mol2_activity': res['j'].map(df['pIC50']),
        'sim': res['sim'],
        'activity_diff': res['dp'],
        'sali': res['sali'],
    })
    return out.reset_index(drop=True)


@st.cache_data
def detect_activity_cliffs_knn(
    df: pd.DataFrame,
    k: int = 8,
    act_thres: float = 1.0,
    radius: int = 2,
    fp_bits: int = 2048,
    min_sim: float = 0.0,
    scaffold_constrained: bool = False,
):
    """Nearest-neighbor based detection: for each molecule, consider its top-k most
    similar neighbors and keep those with delta pIC50 >= threshold.
    """
    # Prepare molecules and fingerprints
    mols = df['SMILES'].apply(Chem.MolFromSmiles)
    gen = get_morgan_generator(radius=radius, fp_size=fp_bits)
    fps = mols.apply(gen.GetFingerprint)
    canon = mols.apply(to_canonical_smiles)

    # Optional scaffold strings for constraint
    scaffold_smiles = None
    if scaffold_constrained:
        scaffold_smiles = mols.apply(lambda m: MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m is not None else None)

    # Precompute similarity matrix (upper triangle sufficient)
    index_list = list(df.index)
    sim_rows = {i: [] for i in index_list}
    for idx_a, i in enumerate(index_list):
        fpi = fps[i]
        for j in index_list:
            if j == i:
                continue
            sim = DataStructs.TanimotoSimilarity(fpi, fps[j])
            if sim >= min_sim:
                sim_rows[i].append((j, sim))

    # For each i, take top-k neighbors
    pairs = {}
    for i in index_list:
        neighbors = sorted(sim_rows[i], key=lambda t: t[1], reverse=True)[:max(0, k)]
        for j, sim in neighbors:
            if scaffold_constrained and scaffold_smiles is not None and scaffold_smiles[i] != scaffold_smiles[j]:
                continue
            # Skip identical canonical molecules
            if canon[i] is not None and canon[i] == canon[j]:
                continue
            dp = abs(df.at[i, 'pIC50'] - df.at[j, 'pIC50'])
            if dp >= act_thres:
                a, b = (i, j) if i < j else (j, i)
                # keep the best sim/delta for the pair if seen twice
                prev = pairs.get((a, b))
                if (prev is None) or (sim > prev[0]):
                    pairs[(a, b)] = (sim, dp)

    if not pairs:
        return pd.DataFrame(columns=[
            'mol1_idx','mol2_idx','mol1_smiles','mol2_smiles','mol1_activity','mol2_activity','sim','activity_diff'
        ])

    rows = []
    for (i, j), (sim, dp) in pairs.items():
        rows.append({
            'mol1_idx': i,
            'mol2_idx': j,
            'mol1_smiles': df.at[i, 'SMILES'],
            'mol2_smiles': df.at[j, 'SMILES'],
            'mol1_activity': df.at[i, 'pIC50'],
            'mol1_canonical': canon[i],
            'mol2_canonical': canon[j],
            'mol2_activity': df.at[j, 'pIC50'],
            'sim': sim,
            'activity_diff': dp,
        })

    out = pd.DataFrame(rows)
    # Sort by activity_diff descending, then similarity descending
    if not out.empty:
        out = out.sort_values(['activity_diff', 'sim'], ascending=[False, False]).reset_index(drop=True)
    return out