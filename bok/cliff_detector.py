# File: activity_cliff_app/cliff_detector.py
# RDKit을 기반으로 분자 유사도(Tanimoto)와 생물학적 활성(IC50 또는 pIC50 등)의 차이를 이용하여 Activity Cliff 후보 쌍을 자동으로 탐지
import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs, rdFMCS
from rdkit.Chem.Scaffolds import MurckoScaffold
from itertools import combinations
import pandas as pd
import time

# Correct import for GetMorganGenerator
from rdkit.Chem import rdFingerprintGenerator
from functools import lru_cache
import math


# --- Similarity Calculation ---

def mcs_similarity(mol1, mol2):
    """Calculates Dice similarity based on the Maximum Common Substructure."""
    if mol1 is None or mol2 is None:
        return 0.0
    mcs_result = rdFMCS.FindMCS([mol1, mol2], timeout=5) # 5-second timeout
    if mcs_result.numAtoms == 0:
        return 0.0
    
    # Dice Similarity formula
    return (2 * mcs_result.numAtoms) / (mol1.GetNumAtoms() + mol2.GetNumAtoms())

def fingerprint_similarity(fp1, fp2):
    """Calculates Tanimoto similarity for fingerprints."""
    return DataStructs.TanimotoSimilarity(fp1, fp2)



def to_canonical_smiles(mol: Chem.Mol) -> str:
    if mol is None:
        return None
    try:
        return Chem.MolToSmiles(mol, isomericSmiles=True, canonical=True)
    except Exception:
        return None

def similarity(fp1, fp2):
    return DataStructs.TanimotoSimilarity(fp1, fp2)

def detect_activity_cliffs(df, sim_thres=0.85, act_thres=1.0, fingerprint_type='ECFP', similarity_metric='fingerprint', progress_callback=None):
    df['mol'] = df['SMILES'].apply(Chem.MolFromSmiles)
    df['canonical_smiles'] = df['mol'].apply(to_canonical_smiles)

    # Pre-calculate fingerprints if needed
    if similarity_metric == 'fingerprint':
        fp_gen = get_fingerprint_generator(radius=2, fp_size=2048, fp_type=fingerprint_type)
        df['fp'] = df['mol'].apply(fp_gen.GetFingerprint)

    results = []
    start_time = time.time()
    pair_combinations = list(combinations(df.index, 2))
    total_pairs = len(pair_combinations)
    for idx, (i, j) in enumerate(pair_combinations):
        if df.loc[i, 'canonical_smiles'] is not None and df.loc[i, 'canonical_smiles'] == df.loc[j, 'canonical_smiles']:
            continue

        sim = 0.0
        if similarity_metric == 'fingerprint':
            sim = fingerprint_similarity(df.loc[i, 'fp'], df.loc[j, 'fp'])
        else: # mcs
            sim = mcs_similarity(df.loc[i, 'mol'], df.loc[j, 'mol'])

        if sim >= sim_thres:
            diff = abs(df.loc[i, 'pIC50'] - df.loc[j, 'pIC50'])
            if diff >= act_thres:
                results.append({
                    'mol1_idx': i, 'mol2_idx': j,
                    'mol1_smiles': df.loc[i, 'SMILES'], 'mol2_smiles': df.loc[j, 'SMILES'],
                    'mol1_canonical': df.loc[i, 'canonical_smiles'], 'mol2_canonical': df.loc[j, 'canonical_smiles'],
                    'mol1_activity': df.loc[i, 'pIC50'], 'mol2_activity': df.loc[j, 'pIC50'],
                    'sim': sim, 'activity_diff': diff,
                })
        if progress_callback:
            progress_callback(idx + 1, total_pairs)
            
    end_time = time.time()
    runtime = end_time - start_time
    return pd.DataFrame(results), runtime


# --- Fingerprint Generation ---

@lru_cache(maxsize=8)
def get_fingerprint_generator(radius: int = 2, fp_size: int = 2048, fp_type: str = 'ECFP'):
    """
    Initializes and returns a fingerprint generator.
    In RDKit, ECFP is equivalent to Morgan with useFeatures=False.
    The 'Morgan' type here can be considered a more generic circular fingerprint.
    ECFP4 (radius=2) is a common standard.
    """
    if fp_type.upper() == 'MORGAN':
        # Morgan can optionally use feature invariants (FCFP)
        return rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)
    elif fp_type.upper() == 'ECFP':
        # ECFP is the default and standard
        return rdFingerprintGenerator.GetMorganGenerator(radius=radius, fpSize=fp_size)
    else:
        raise ValueError(f"Unsupported fingerprint type: {fp_type}")


def detect_activity_cliffs_sali(
    df: pd.DataFrame,
    top_n: int = 200,
    radius: int = 2,
    fp_bits: int = 2048,
    sim_floor: float = 0.5,
    scaffold_constrained: bool = False,
    fingerprint_type: str = 'ECFP',
    similarity_metric: str = 'fingerprint',
    progress_callback=None,
):
    """Rank pairs by SALI = delta_pIC50 / (1 - similarity)."""
    mols = df['SMILES'].apply(Chem.MolFromSmiles)
    canon = mols.apply(to_canonical_smiles)

    fps = None
    if similarity_metric == 'fingerprint':
        gen = get_fingerprint_generator(radius=radius, fp_size=fp_bits, fp_type=fingerprint_type)
        fps = mols.apply(gen.GetFingerprint)

    scaffold_smiles = None
    if scaffold_constrained:
        scaffold_smiles = mols.apply(lambda m: MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m is not None else None)

    rows = []
    start_time = time.time()
    total_pairs = len(df.index) * (len(df.index) - 1) // 2
    count = 0
    for i in df.index:
        for j in df.index:
            if j <= i: continue
            count += 1
            if canon[i] is not None and canon[i] == canon[j]: continue
            if scaffold_constrained and scaffold_smiles[i] != scaffold_smiles[j]: continue

            sim = 0.0
            if similarity_metric == 'fingerprint':
                sim = fingerprint_similarity(fps[i], fps[j])
            else: # mcs
                sim = mcs_similarity(mols[i], mols[j])

            if sim < sim_floor: continue
            
            dp = abs(df.at[i, 'pIC50'] - df.at[j, 'pIC50'])
            denom = max(1e-6, 1.0 - sim)
            sali = dp / denom
            rows.append((i, j, sim, dp, sali))
            if progress_callback:
                progress_callback(count, total_pairs)
    
    end_time = time.time()
    runtime = end_time - start_time

    if not rows:
        return pd.DataFrame(columns=[
            'mol1_idx','mol2_idx','mol1_smiles','mol2_smiles','mol1_activity','mol2_activity','sim','activity_diff','sali'
        ]), runtime

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
    return out.reset_index(drop=True), runtime


def detect_activity_cliffs_knn(
    df: pd.DataFrame,
    k: int = 8,
    act_thres: float = 1.0,
    radius: int = 2,
    fp_bits: int = 2048,
    min_sim: float = 0.0,
    scaffold_constrained: bool = False,
    fingerprint_type: str = 'ECFP',
    similarity_metric: str = 'fingerprint',
    progress_callback=None,
):
    """Nearest-neighbor based detection."""
    mols = df['SMILES'].apply(Chem.MolFromSmiles)
    canon = mols.apply(to_canonical_smiles)

    scaffold_smiles = None
    if scaffold_constrained:
        scaffold_smiles = mols.apply(lambda m: MurckoScaffold.MurckoScaffoldSmiles(mol=m) if m is not None else None)

    index_list = list(df.index)
    sim_rows = {i: [] for i in index_list}

    # Similarity calculation loop
    fps = None
    if similarity_metric == 'fingerprint':
        gen = get_fingerprint_generator(radius=radius, fp_size=fp_bits, fp_type=fingerprint_type)
        fps = mols.apply(gen.GetFingerprint)
    
    start_time = time.time()
    total_calcs = len(index_list) * len(index_list)
    count = 0
    for i in index_list:
        for j in index_list:
            if j == i: continue
            count += 1
            sim = 0.0
            if similarity_metric == 'fingerprint':
                sim = fingerprint_similarity(fps[i], fps[j])
            else: # mcs
                sim = mcs_similarity(mols[i], mols[j])

            if sim >= min_sim:
                sim_rows[i].append((j, sim))
            if progress_callback:
                progress_callback(count, total_calcs)
    
    end_time = time.time()
    runtime = end_time - start_time

    pairs = {}
    for i in index_list:
        neighbors = sorted(sim_rows[i], key=lambda t: t[1], reverse=True)[:max(0, k)]
        for j, sim in neighbors:
            if scaffold_constrained and scaffold_smiles is not None and scaffold_smiles[i] != scaffold_smiles[j]: continue
            if canon[i] is not None and canon[i] == canon[j]: continue
            
            dp = abs(df.at[i, 'pIC50'] - df.at[j, 'pIC50'])
            if dp >= act_thres:
                a, b = (i, j) if i < j else (j, i)
                prev = pairs.get((a, b))
                if (prev is None) or (sim > prev[0]):
                    pairs[(a, b)] = (sim, dp)

    if not pairs:
        return pd.DataFrame(columns=[
            'mol1_idx','mol2_idx','mol1_smiles','mol2_smiles','mol1_activity','mol2_activity','sim','activity_diff'
        ]), runtime

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
    return out, runtime
