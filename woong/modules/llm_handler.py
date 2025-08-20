
from openai import OpenAI

def generate_hypothesis(
    api_key: str,
    smiles1: str,
    activity1: float,
    smiles2: str,
    activity2: float,
    structural_difference_description: str,
    model: str = "gpt-4-turbo",
    max_tokens: int = 4096,
) -> str:
    """
    두 분자의 구조-활성 관계(SAR)에 대한 화학적 가설을 생성합니다.

    LLM에 다음을 요청합니다:
    - 제공된 두 화합물 데이터를 기반으로, 반증 가능하고 실험적으로 검증 가능한 가설을 생성합니다.
    - 모든 주장은 구조적 사실(원자, 결합, 작용기 등의 변화)에 근거해야 합니다.
    - 출력은 엄격한 JSON 형식이어야 합니다.

    Args:
        api_key (str): OpenAI API 키.
        smiles1 (str): 첫 번째 화합물의 SMILES 문자열.
        activity1 (float): 첫 번째 화합물의 활성도.
        smiles2 (str): 두 번째 화합물의 SMILES 문자열.
        activity2 (float): 두 번째 화합물의 활성도.
        structural_difference_description (str): 두 화합물 간의 구조적 차이 요약.
        model (str): 사용할 OpenAI 모델 이름.
        max_tokens (int): 생성할 최대 토큰 수.

    Returns:
        str: LLM이 생성한 가설 (JSON 형식의 문자열).
    """
    client = OpenAI(api_key=api_key)

    # LLM에 전달할 프롬프트. 명확성과 가독성을 위해 여러 줄 문자열로 구성합니다.
    prompt = f'''
[STRICT SAR HYPOTHESIS REQUEST]

**역할**: 당신은 숙련된 의약화학자입니다. 아래 데이터만을 근거로 구조-활성 상관(SAR)에 대한
'반증 가능'하고 '실험으로 검증 가능한' 가설을 생성하세요. 일반론적 추측은 금지하며,
모든 주장은 반드시 구조적 사실(변경된 원자/결합/치환기)과 연결되어야 합니다.

**데이터**
- **화합물 1**: SMILES={smiles1}; activity={activity1}
- **화합물 2**: SMILES={smiles2}; activity={activity2}
- **구조 차이 요약(MMP 관점)**: {structural_difference_description}
- **지표/단위**: 미지정(사용자 제공 가정 필요). 값이 낮을수록 좋다는(IC50/Ki) 가정을 기본으로 하되, 다르면 '가정:'에 명시
- **Tanimoto(유사도/거리)**: 미제공 시 N/A
- **Assay 맥락**: 미제공(필요시 가정 명시)

**분석 지침** (각 항목은 '주장 → 근거(구조적 사실) → 예상 효과 방향' 형태로 1~2문장)
1)  **활성 비교**: 어느 화합물이 더 유리한 활성인지 먼저 선언하고, 가능하면 Δ(pAct)를 계산·표기.
    - pAct = -log10(IC50[M]) 또는 동등한 로그척도. 단위 불명확 시 '가정:'으로 단위 명시 후 계산
2)  **기전 분석**: 변경된 부분(코어/링/링커/치환기) 중심으로 다음 기전 관점을 검토.
    - 전자효과(유도/공명), 입체효과(부피/혼잡/접근성), 수소결합(HBD/HBA 변화),
      소수성/극성(logP/TPSA 방향 추정), 이온화상태(pKa 방향 추정),
      컨포메이션(회전 자유도/토션 페널티), π-상호작용/아로마틱성
3)  **결합 모드**: 결합 모드 추정은 일반적 포켓 유형 수준으로 제한(특정 잔기 단정 금지).
4)  **반증 가능성**: 최소 2개의 대안 가설과 각 가설의 반례 조건 제시.
5)  **설계 제안(검증용 MMP 스캔)**: 최소 3개 치환 제안.
    - design(치환/변경) · 예상 효과(↑/↓, 약/중/강) · 근거 · 검증지표(예: pIC50)
6)  **ADMET 위험 신호**: 용해도/투과성/P-gp/CYP/반응성 등 1~3개 위험 신호와 구조적 근거.
7)  **가정과 한계**: 사용한 가정(단위, pKa 추정 등)과 데이터 한계를 명시.
8)  **추론 범위**: 제공 데이터와 일반 화학 상식 범위 내에서만 추정(환각 금지).

**출력 형식** (JSON만; 키 유지, 한국어 값):
{{
  "more_active": "<compound_1 또는 compound_2>",
  "delta_pAct_explained": "<가능하면 ΔpAct 계산 또는 불가 사유>",
  "primary_hypothesis": "<가장 설득력 있는 단일 가설(1~3문장)>",
  "mechanistic_rationale": {{
    "electronic": "<주장. 근거: ...>",
    "steric": "<주장. 근거: ...>",
    "hydrogen_bonding": "<주장. 근거: ...>",
    "hydrophobic_polar": "<주장. 근거: ...>",
    "ionization_pKa": "<주장. 근거: ...>",
    "conformation": "<주장. 근거: ...>",
    "aromatic_pi": "<선택>"
  }},
  "observations_from_pair": [
    "<구조 차이에서 직접 관찰된 사실 1>",
    "<사실 2>"
  ],
  "counter_hypotheses": [
    "<대안 가설 1 + 반례 조건>",
    "<대안 가설 2 + 반례 조건>"
  ],
  "design_suggestions": [
    {{"design": "<치환 제안 1>", "expected_effect": "<↑/↓, 약/중/강>", "rationale": "<근거>", "validation_metric": "<예: pIC50>"}},
    {{"design": "<치환 제안 2>", "expected_effect": "<↑/↓, 약/중/강>", "rationale": "<근거>", "validation_metric": "<예: pIC50>"}},
    {{"design": "<치환 제안 3>", "expected_effect": "<↑/↓, 약/중/강>", "rationale": "<근거>", "validation_metric": "<예: pIC50>"}}
  ],
  "admet_flags": ["<우려 1>", "<우려 2>"],
  "assumptions_and_limits": ["<가정 1>", "<한계 1>"],
  "confidence": 0.0
}}

지침을 엄격히 준수하여 **JSON만** 출력하세요(설명·코드블록 금지).
'''

    system_message = "You are a helpful assistant that is an expert in medicinal chemistry and drug discovery."

    try:
        # OpenAI 최신 API는 `messages` 파라미터를 사용합니다.
        chat_result = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user", "content": prompt},
            ],
            max_tokens=max_tokens,
            temperature=0.7,
            response_format={"type": "json_object"},  # JSON 모드 활성화
        )
        return chat_result.choices[0].message.content.strip()
    except Exception as e:
        return f"LLM 가설 생성 중 오류 발생: {e}"
