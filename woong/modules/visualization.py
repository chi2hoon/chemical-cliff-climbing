from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.rdFMCS import FindMCS
from PIL import Image

def visualize_structure_difference(smiles1: str, smiles2: str, legend1: str = "Molecule 1", legend2: str = "Molecule 2") -> Image:
    """
    두 SMILES 문자열을 비교하여 구조적 차이점을 하이라이트한 이미지를 생성합니다.

    Args:
        smiles1 (str): 첫 번째 분자의 SMILES 문자열
        smiles2 (str): 두 번째 분자의 SMILES 문자열
        legend1 (str): 첫 번째 분자의 레전드
        legend2 (str): 두 번째 분자의 레전드

    Returns:
        PIL.Image: 두 분자가 나란히 그려지고 차이점이 하이라이트된 이미지 객체
    """
    mols = []
    mol1 = None
    mol2 = None

    try:
        mol1 = Chem.MolFromSmiles(smiles1)
    except Exception as e:
        print(f"ERROR: smiles1 ({smiles1})로부터 mol1 생성 중 오류 발생: {e}")
        raise ValueError(f"유효하지 않은 SMILES 문자열 (Molecule 1): {smiles1}. 오류: {e}")

    try:
        mol2 = Chem.MolFromSmiles(smiles2)
    except Exception as e:
        print(f"ERROR: smiles2 ({smiles2})로부터 mol2 생성 중 오류 발생: {e}")
        raise ValueError(f"유효하지 않은 SMILES 문자열 (Molecule 2): {smiles2}. 오류: {e}")

    # 분자 객체의 유효성 검사 및 원자 존재 여부 확인
    if mol1 is None or mol1.GetNumAtoms() == 0:
        raise ValueError(f"유효하지 않거나 빈 SMILES 문자열 (Molecule 1): {smiles1}")
    if mol2 is None or mol2.GetNumAtoms() == 0:
        raise ValueError(f"유효하지 않거나 빈 SMILES 문자열 (Molecule 2): {smiles2}")

    mols.append(mol1)
    mols.append(mol2)

    # 하이라이트할 원자 목록 초기화 (항상 두 개의 빈 리스트로 시작)
    highlight_atom_lists = [[], []]

    # 두 분자 간의 최대 공통 부분구조(MCS)를 찾습니다.
    mcs_result = None
    try:
        mcs_result = FindMCS(mols, timeout=5)
    except Exception as e:
        print(f"ERROR: FindMCS 실행 중 오류 발생: {e}")
        # MCS 찾기 실패 시 하이라이트 없이 진행하거나, 오류를 다시 발생시킬 수 있습니다.
        # 여기서는 하이라이트 없이 진행하도록 합니다.
        pass # MCS 찾기 실패 시 하이라이트 없이 진행

    if mcs_result and mcs_result.numAtoms >= 1: # MCS가 존재할 경우에만 하이라이트 로직 실행
        mcs_pattern = None
        try:
            mcs_pattern = Chem.MolFromSmarts(mcs_result.smartsString)
        except Exception as e:
            print(f"ERROR: smartsString로부터 mcs_pattern 생성 중 오류 발생: {e}")
            mcs_pattern = None # 오류 발생 시 패턴을 None으로 설정하여 하이라이트 로직 건너뛰기

        if mcs_pattern is not None: # MCS 패턴이 유효한 경우에만 매칭 시도
            try:
                match1 = mols[0].GetSubstructMatch(mcs_pattern)
                match2 = mols[1].GetSubstructMatch(mcs_pattern)

                # MCS에 포함되지 않는 원자들을 하이라이트 목록에 추가
                highlight_atom_lists[0] = [i for i in range(mols[0].GetNumAtoms()) if i not in match1]
                highlight_atom_lists[1] = [i for i in range(mols[1].GetNumAtoms()) if i not in match2]
            except Exception as e:
                print(f"ERROR: SubstructMatch 또는 하이라이트 목록 생성 중 오류 발생: {e}")

    # 두 분자를 나란히 그리고, 차이점을 하이라이트 처리합니다.
    img = None
    try:
        img = Draw.MolsToGridImage(
            mols,
            molsPerRow=2,
            subImgSize=(350, 350),
            legends=[legend1, legend2],
            highlightAtomLists=highlight_atom_lists
        )
    except Exception as e:
        print(f"ERROR: MolsToGridImage 생성 중 오류 발생: {e}")
        raise RuntimeError(f"이미지 생성 중 오류 발생: {e}")

    return img