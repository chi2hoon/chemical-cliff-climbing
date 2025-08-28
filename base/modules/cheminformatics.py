import pandas as pd
from rdkit import Chem
from rdkit.DataStructs import BulkTanimotoSimilarity
from itertools import combinations

def find_activity_cliffs(df: pd.DataFrame, smiles_col: str, activity_col: str, similarity_threshold: float, activity_diff_threshold: float, higher_is_better: bool = True):
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
    mols = [Chem.MolFromSmiles(smi) for smi in df[smiles_col]]
    fps = [Chem.RDKFingerprint(mol) for mol in mols if mol is not None]

    cliff_pairs = []

    # 모든 분자 쌍에 대해 반복
    for i, j in combinations(range(len(fps)), 2):
        # Tanimoto 유사도 계산
        similarity = BulkTanimotoSimilarity(fps[i], [fps[j]])[0]

        if similarity >= similarity_threshold:
            # 활성도 차이 계산
            activity_i = df.loc[i, activity_col]
            activity_j = df.loc[j, activity_col]
            activity_diff = abs(activity_i - activity_j)

            if activity_diff >= activity_diff_threshold:
                cliff_pairs.append({
                    'SMILES_1': df.loc[i, smiles_col],
                    'Activity_1': activity_i,
                    'SMILES_2': df.loc[j, smiles_col],
                    'Activity_2': activity_j,
                    'Similarity': similarity,
                    'Activity_Diff': activity_diff
                })

    return pd.DataFrame(cliff_pairs)
