
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
    prompt = f'''### 핵심 지침 (Core Instructions)
* **엄격성이 최우선:** 당신의 1차 목표는 **완전하고 엄격하게 정당화된** SAR 가설을 제시하는 것입니다. 논리 비약 없이 각 단계의 주장에 **구조적 사실**과 **물리·유기화학적 원리**를 명시적으로 연결하세요. 단순한 직감이나 일반론적 상투구는 실패로 간주됩니다.
* **완결성에 대한 정직성:** 완전한 가설(충분히 검증 가능한 일관된 설명)을 제시할 수 없으면 **추측으로 채우지 마십시오**. 그 대신, **유의미한 부분 결과**만 제시하세요.

### 문제 데이터 (Problem Data)
- 화합물 1 SMILES: {smiles1}
- 화합물 2 SMILES: {smiles2}
- 활성도 지표 및 값: = {activity1} (화합물 1), {activity2}(화합물 2)
- 유사도: Tanimoto = {tanimoto_or_NA}
- 구조 차이 요약(MMP 관점): {structural_difference_description}
- 규칙: 값이 낮을수록 좋은 지표(IC50/Ki 등)라면 이를 명기하고, 그렇지 않으면 해당 가정을 **Summary**의 “가정(Assumptions)”에 명확히 적으세요.

### 출력 형식 (Output Format)
**지침: 다음의 키를 사용하는 JSON 객체만을 반환하세요. 값은 한국어로 작성합니다. `Generation.md`의 논리 구조를 따르되, 이 JSON 형식에 맞춰 내용을 채워주세요.**

```json
{{
  "more_active": "<activity_2의 값이 더 높으면 '화합물 2', 그렇지 않으면 '화합물 1'>",
  "delta_pAct_explained": "<가능하면 pAct를 계산하고 그 차이(ΔpAct)를 설명. 불가능하면 그 사유를 명시>",
  "primary_hypothesis": "<가장 설득력 있는 단일 가설(1~3문장). Generation.md의 'Summary > a. Verdict'에 해당>",
  "mechanistic_rationale": {{
    "electronic": "<전자 효과(유도/공명)에 대한 주장과 근거>",
    "steric": "<입체 효과(부피/혼잡)에 대한 주장과 근거>",
    "hydrogen_bonding": "<수소 결합(HBD/HBA) 변화에 대한 주장과 근거>",
    "hydrophobic_polar": "<소수성/극성(logP/TPSA) 변화에 대한 주장과 근거>",
    "ionization_pKa": "<이온화 상태(pKa) 변화에 대한 주장과 근거>",
    "conformation": "<컨포메이션(회전 자유도/토션) 변화에 대한 주장과 근거>",
    "aromatic_pi": "<π-상호작용/평면성 변화에 대한 주장과 근거>"
  }},
  "counter_hypotheses": [
    "<대안 가설 1 + 해당 가설이 성립하지 않는 반례 조건>",
    "<대안 가설 2 + 해당 가설이 성립하지 않는 반례 조건>"
  ],
  "design_suggestions": [
    {{"design": "<MMP 스캔 기반 설계 제안 1(SMILES)>", "expected_effect": "<예상 효과 ↑/↓, 약/중/강>", "rationale": "<설계의 근거>", "validation_metric": "<검증에 사용할 지표, 예: pIC50>"}},
    {{"design": "<설계 제안 2(SMILES)>", "expected_effect": "<예상 효과 ↑/↓, 약/중/강>", "rationale": "<설계의 근거>", "validation_metric": "<검증에 사용할 지표, 예: pIC50>"}},
    {{"design": "<설계 제안 3(SMILES)>", "expected_effect": "<예상 효과 ↑/↓, 약/중/강>", "rationale": "<설계의 근거>", "validation_metric": "<검증에 사용할 지표, 예: pIC50>"}}
  ],
  "admet_flags": [
      "<예상되는 ADMET 위험 1과 그 구조적 근거>",
      "<예상되는 ADMET 위험 2와 그 구조적 근거>"
  ],
  "assumptions_and_limits": [
    "<분석에 사용된 가정 1 (예: 활성도 단위 가정)>",
    "<분석의 한계 1>"
  ],
  "confidence": 0.85
}}
```

### 자기 교정 지침 (Self-Correction Instruction)
제출 전, 응답이 위의 JSON 스키마를 정확히 따르는지, 그리고 `Generation.md`의 엄격성/완결성 기준을 만족하는지 최종 점검하세요.
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
