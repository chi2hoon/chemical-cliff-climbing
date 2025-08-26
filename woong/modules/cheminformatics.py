import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import rdFingerprintGenerator
from itertools import combinations

def find_activity_cliffs(
    df: pd.DataFrame,
    smiles_col: str,
    activity_col: str,
    similarity_threshold: float,
    activity_diff_threshold: float,
    fingerprint: str = "morgan2",
):
    """
    데이터프레임에서 Activity Cliff 쌍을 찾습니다.

    Args:
        df (pd.DataFrame): SMILES와 활성도 데이터가 포함된 데이터프레임
        smiles_col (str): SMILES 데이터가 있는 컬럼명
        activity_col (str): 활성도 데이터가 있는 컬럼명
        similarity_threshold (float): 구조 유사도 임계값 (0.0 ~ 1.0)
        activity_diff_threshold (float): 활성도 차이 임계값

    Returns:
        pd.DataFrame: 찾은 Activity Cliff 쌍들의 정보를 담은 데이터프레임
    """
    # 1) 입력 정리: SMILES를 문자열로 캐스팅하고 공백/결측 제거, 활성도 숫자화
    work = df.copy()
    work[smiles_col] = work[smiles_col].astype(str).str.strip()
    work[activity_col] = pd.to_numeric(work[activity_col], errors='coerce')
    invalid_smiles = work[smiles_col].str.lower().isin(["", "nan", "none", "null"])
    work = work[(~invalid_smiles) & (work[activity_col].notna())].reset_index(drop=True)

    # 2) 유효한 분자만 선택하면서 인덱스 매핑 유지
    mols = [Chem.MolFromSmiles(s) for s in work[smiles_col]]
    valid_idx = [k for k, m in enumerate(mols) if m is not None]
    if len(valid_idx) < 2:
        return pd.DataFrame(columns=['SMILES_1','Activity_1','SMILES_2','Activity_2','Similarity','Activity_Diff'])

    # Prefer the new MorganGenerator to avoid deprecation warnings. Fallback to legacy API if needed.
    morgan_generator = None
    if fingerprint == "morgan2":
        try:
            morgan_generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
        except Exception:
            morgan_generator = None

    def _make_fp(mol):
        if fingerprint == "morgan2":
            if morgan_generator is not None:
                return morgan_generator.GetFingerprint(mol)
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        return Chem.RDKFingerprint(mol)

    fps = [_make_fp(mols[k]) for k in valid_idx]

    # 3) 쌍 비교 시 원본 행 위치를 매핑하여 정합성 유지
    cliff_pairs = []
    for a, b in combinations(range(len(fps)), 2):
        i = valid_idx[a]
        j = valid_idx[b]
        similarity = DataStructs.TanimotoSimilarity(fps[a], fps[b])
        if similarity < similarity_threshold:
            continue

        activity_i = float(work.loc[i, activity_col])
        activity_j = float(work.loc[j, activity_col])
        activity_diff = abs(activity_i - activity_j)
        if activity_diff < activity_diff_threshold:
            continue

        cliff_pairs.append({
            'SMILES_1': work.loc[i, smiles_col],
            'Activity_1': activity_i,
            'SMILES_2': work.loc[j, smiles_col],
            'Activity_2': activity_j,
            'Similarity': similarity,
            'Activity_Diff': activity_diff
        })

    return pd.DataFrame(cliff_pairs)
