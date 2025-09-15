import pandas as pd
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
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
    # 1) SMILES → Mol 변환과 유효 인덱스 매핑(인덱스 정합성 보장)
    smiles_series = df[smiles_col].astype(str).tolist()
    mols = [Chem.MolFromSmiles(smi) for smi in smiles_series]
    valid_idx = [i for i, m in enumerate(mols) if m is not None]
    if not valid_idx:
        return pd.DataFrame(columns=['SMILES_1','Activity_1','SMILES_2','Activity_2','Similarity','Activity_Diff'])
    fps = [Chem.RDKFingerprint(mols[i]) for i in valid_idx]

    # 2) 대규모 쌍 계산 방지용 간단한 가드(과도한 계산 시 일부만 샘플링)
    # 유사도 계산은 O(n^2)이므로 n이 큰 경우 스트림릿 UI가 멈출 수 있음
    max_pairs = 400_000  # 약 40만 쌍 상한
    n = len(fps)
    est_pairs = n * (n - 1) // 2
    use_sampling = est_pairs > max_pairs

    cliff_pairs = []

    # 3) 모든 분자 쌍에 대해 반복(인덱스 매핑 적용)
    import random
    index_range = list(range(n))
    if use_sampling:
        random.seed(42)
    for ii in range(n):
        # 샘플링 모드에서는 ii 이후 j 후보를 무작위 일부만 검사
        if use_sampling:
            js = random.sample(index_range[ii + 1 :], k=min(512, max(0, n - ii - 1)))
        else:
            js = index_range[ii + 1 :]
        for jj in js:
            similarity = TanimotoSimilarity(fps[ii], fps[jj])
            if similarity < similarity_threshold:
                continue
            # 원본 df의 실제 행 인덱스 매핑
            i_row = valid_idx[ii]
            j_row = valid_idx[jj]
            activity_i = df.iloc[i_row][activity_col]
            activity_j = df.iloc[j_row][activity_col]
            try:
                activity_diff = abs(float(activity_i) - float(activity_j))
            except Exception:
                # 숫자 변환 실패 시 스킵
                continue
            if activity_diff >= activity_diff_threshold:
                cliff_pairs.append({
                    'SMILES_1': df.iloc[i_row][smiles_col],
                    'Activity_1': float(activity_i) if activity_i is not None else None,
                    'SMILES_2': df.iloc[j_row][smiles_col],
                    'Activity_2': float(activity_j) if activity_j is not None else None,
                    'Similarity': similarity,
                    'Activity_Diff': activity_diff
                })

    return pd.DataFrame(cliff_pairs)
