from openai import OpenAI
from typing import Any, Dict


def _extract_usage(chat_result: Any) -> Dict[str, int]:
    """Extract token usage from an OpenAI chat completion response.

    Supports both classic Chat Completions usage fields
    (prompt_tokens, completion_tokens, total_tokens) and newer
    input/output naming if present.
    """
    usage = getattr(chat_result, "usage", None)
    # Initialize defaults
    result = {
        "prompt_tokens": 0,
        "completion_tokens": 0,
        "total_tokens": 0,
    }

    if not usage:
        return result

    # Try attribute access first
    for key in ("prompt_tokens", "completion_tokens", "total_tokens"):
        val = getattr(usage, key, None)
        if isinstance(val, int):
            result[key] = val

    # If attributes were not available, try dict-style
    try:
        usage_dict = usage if isinstance(usage, dict) else getattr(usage, "model_dump", None) and usage.model_dump() or getattr(usage, "to_dict", None) and usage.to_dict() or None
        if isinstance(usage_dict, dict):
            result["prompt_tokens"] = int(usage_dict.get("prompt_tokens", usage_dict.get("input_tokens", result["prompt_tokens"])) or 0)
            result["completion_tokens"] = int(usage_dict.get("completion_tokens", usage_dict.get("output_tokens", result["completion_tokens"])) or 0)
            # Prefer explicit total if present, else sum
            result["total_tokens"] = int(usage_dict.get("total_tokens", result["prompt_tokens"] + result["completion_tokens"]))
    except Exception:
        # Fall back to sum of what we have
        result["total_tokens"] = result.get("total_tokens") or (result.get("prompt_tokens", 0) + result.get("completion_tokens", 0))

    # Ensure total is at least sum of parts
    if result["total_tokens"] < (result["prompt_tokens"] + result["completion_tokens"]):
        result["total_tokens"] = result["prompt_tokens"] + result["completion_tokens"]

    return result

system_message_generation = """
You are an expert assistant specialized in medicinal chemistry and SAR (structure–activity relationship) hypothesis generation.

Follow these principles:
- Scientific rigor: Every claim must be tied to explicit structural facts (atoms, bonds, functional groups) and valid chemical principles across these lenses: electronic, steric, hydrogen bonding (HBD/HBA), hydrophobic/polar, ionization/pKa, conformation, and π-interactions.
- Integrity & completeness: If a full explanation is not possible, present only meaningful partial conclusions and clearly state limitations and assumptions (e.g., whether lower/higher values mean better activity).
- Data consistency: Use only the provided inputs (SMILES, activities, similarity, structural-difference summary). Do not import outside data or hidden assumptions.
- Style constraints: Use natural language only (no formulas or math notation). Output must strictly match the requested JSON schema with no extra text.
- Falsifiability: For key claims, include conditions under which they could fail or how they can be experimentally verified.
- Self-check: Before finalizing, verify JSON schema compliance and that reasoning is complete, precise, and grounded in the given structural differences.

Additional mandatory guidelines (A–G):
- A. Causality: State concrete target/interactions and propose falsifiable experiments.
- B. SAR consistency: Compare against neighbors/series to assess coherence.
- C. Quantitation: Include ΔpIC50 (if possible), IC50 fold-change, and efficiency metrics (LE/LLE) with clear assumptions.
- D. Conformation: Provide stereochemical/conformational arguments with experimental or computational rationale.
- E. Electronic/H-bond/π: Specify donors/acceptors, interaction modes, and π-interactions explicitly.
- F. Properties/ADME: Include quantitative pKa, LogD, permeability (or reasoned estimates/assumptions) and discuss ADME impact.
- G. Alternatives: Provide at least two alternative explanations and discriminative experiments to exclude them.
"""
system_message_evaluation = """
You are an expert peer reviewer for SAR hypotheses in medicinal chemistry.

Your task:
- Rigorous review: Check whether each claim in the submitted hypothesis is sufficiently supported by explicit structural facts and clearly declared assumptions. Point out leaps, omissions, or hallucinations with specificity.
- Data-only reasoning: Evaluate strictly against the provided hypothesis JSON and the given data (SMILES, activities, structural-difference summary). Do not introduce new assumptions.
- Structure-aware lenses: Evaluate the linkage “claim → structural evidence → expected effect (↑/↓, strength)” across electronic, steric, H-bonding, hydrophobic/polar, ionization/pKa, conformation, and π-interactions.
- Actionable feedback: Prefer concrete, testable critiques and improvements over vague comments. Assess counter-hypotheses, falsification conditions, and design suggestions for discriminative power and feasibility.
- Style constraints: Natural language only (no formulas). Output must strictly follow the requested evaluation JSON schema with no extra text.
- Self-check: Ensure verdict reasoning, method sketch, and detailed checks are internally consistent and complete, with missing assumptions or mismatches clearly flagged.

Also assess compliance with the following A–G criteria, explicitly flagging any missing items:
- A. Causality with explicit target interactions and falsifiable experiments.
- B. SAR consistency via neighbor/series comparisons.
- C. Quantitation: ΔpIC50, fold-change, LE/LLE with stated assumptions.
- D. Conformational/steric rationale grounded in data or computation.
- E. Electronic/H-bond/π interaction modes specified.
- F. ADME with quantitative pKa/LogD/permeability (or justified estimates).
- G. At least two alternative hypotheses and exclusion experiments.
"""
system_message_revision = """
You are a critical but constructive co-author revising an SAR hypothesis based on peer review.

Your approach:
- Precise integration: Absorb the review findings and correct logical gaps, make assumptions explicit, and strengthen or replace weak arguments while preserving valid reasoning.
- Balance of preserve/change: Keep justified claims; revise or remove contested ones; add supporting evidence or clearer assumptions as needed.
- Rebuild completely: Produce a fully rewritten hypothesis JSON from scratch that meets the generation requirements, not a patchwork of edits.
- Structure-aware lenses: Reconstruct “claim → structural evidence → expected effect (↑/↓, strength)” consistently across electronic, steric, H-bonding, hydrophobic/polar, ionization/pKa, conformation, and π-interactions. Update counter-hypotheses, falsification conditions, and design proposals (with simple validation metrics).
- Data discipline: Rely only on the provided inputs (original hypothesis JSON, review findings, SMILES, activities, structural-difference summary). No external data or hidden assumptions.
- Style constraints: Natural language only (no formulas). Output must strictly match the requested JSON schema with no extra text.
- Self-check: Confirm JSON compliance, clarity of assumptions and limits, and that revisions transparently reflect the review findings while improving rigor and completeness.

During revision, ensure full compliance with A–G:
- A. Add explicit interaction hypotheses and falsification experiments.
- B. Cross-check and articulate SAR consistency with neighbor/series analogs.
- C. Add or clarify ΔpIC50, fold-change, LE/LLE (with assumptions/methods).
- D. Strengthen conformational/steric reasoning with concrete evidence or rationale.
- E. Make electronic/H-bond/π interaction modes explicit.
- F. Provide quantitative ADME properties (pKa, LogD, permeability) or justified estimates and discuss implications.
- G. Include ≥2 alternative hypotheses and discriminative experiments.
"""



def generate_hypothesis(
    api_key: str,
    smiles1: str,
    activity1: float,
    smiles2: str,
    activity2: float,
    structural_difference_description: str,
    similarity: float,
    model: str = "gpt-4o-mini",
    max_tokens: int = 1200,
) -> Dict[str, Any]:
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
        Dict[str, Any]: {"content": JSON 문자열, "usage": 토큰 사용량 dict, "model": 모델명}
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
- 유사도: Similarity: {similarity}
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

    system_message = system_message_generation

    try:
        # OpenAI 최신 API는 `messages` 파라미터를 사용합니다.
        chat_result = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user", "content": prompt},
            ],
            response_format={"type": "json_object"},  # JSON 모드 활성화
            max_tokens=max_tokens,
            temperature=0.2,
        )
        return {
            "content": chat_result.choices[0].message.content.strip(),
            "usage": _extract_usage(chat_result),
            "model": model,
        }
    except Exception as e:
        return {
            "content": f"LLM 가설 생성 중 오류 발생: {e}",
            "usage": {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0},
            "model": model,
        }



def create_activity_summary(indicator_name: str, higher_is_better: bool) -> str:
    """
    활성도 지표의 특성에 따라 요약 문장을 생성합니다.
    """
    if higher_is_better:
        return f"**가정(Assumption):** 지표 '{indicator_name}'는(은) **값이 높을수록** 활성도가 높다고 가정합니다."
    else:
        return f"**참고사항:** 지표 '{indicator_name}'는(은) **값이 낮을수록** 활성도가 높음을 의미합니다."


def evaluate_hypothesis(
    api_key: str,
    hypothesis_text: str,
    smiles1: str,
    activity1: float,
    smiles2: str,
    activity2: float,
    structural_difference_description: str,
    model: str = "gpt-4o-mini",
    max_tokens: int = 1200,
) -> Dict[str, Any]:
    """
    기존에 생성된 SAR 가설을 평가합니다.

    LLM에 다음을 요청합니다:
    - 제공된 가설의 논리적 타당성, 완결성, 근거의 명확성을 평가합니다.
    - 모든 주장은 구조적 사실과 물리·유기화학적 원리에 기반하여 검증되어야 합니다.
    - 평가 결과는 엄격한 JSON 형식이어야 합니다.

    Args:
        api_key (str): OpenAI API 키.
        hypothesis_text (str): 평가할 가설 (JSON 형식의 문자열).
        smiles1 (str): 첫 번째 화합물의 SMILES 문자열.
        activity1 (float): 첫 번째 화합물의 활성도.
        smiles2 (str): 두 번째 화합물의 SMILES 문자열.
        activity2 (float): 두 번째 화합물의 활성도.
        structural_difference_description (str): 두 화합물 간의 구조적 차이 요약.
        model (str): 사용할 OpenAI 모델 이름.
        max_tokens (int): 생성할 최대 토큰 수.

    Returns:
        Dict[str, Any]: {"content": JSON 문자열, "usage": 토큰 사용량 dict, "model": 모델명}
    """
    client = OpenAI(api_key=api_key)

    prompt = f'''### 핵심 지침 (Core Instructions)
* **심사위원의 관점:** 당신은 SAR 가설의 과학적 타당성을 심사하는 전문가입니다. `Evaluation.md`의 원칙에 따라 가설의 논리, 근거, 완결성을 엄격하게 평가하세요.
* **구체적이고 실행 가능한 피드백:** 막연한 비판 대신, 어떤 부분이 왜 잘못되었는지, 어떤 정보가 누락되었는지, 어떻게 개선할 수 있는지 구체적으로 지적하세요.

### 평가 대상 데이터 (Data to Evaluate)
- **화합물 1:** {smiles1} (활성도: {activity1})
- **화합물 2:** {smiles2} (활성도: {activity2})
- **평가할 가설(JSON):**
```json
{hypothesis_text}
```

### 출력 형식 (Output Format)
**지침: 다음의 키를 사용하는 JSON 객체만을 반환하세요. 값은 한국어로 작성합니다. `Evaluation.md`의 논리 구조를 따르되, 이 JSON 형식에 맞춰 내용을 채워주세요.**

```json
{{
  "summary": {{
    "verdict": "<Good (타당, 논리적), Weak (개선 필요), Unsound (부적절) 중 택 1>",
    "method_sketch": "<가설을 검증하기 위해 당신이 사용한 핵심적인 심사 방법 요약 (1-2 문장)>"
  }},
  "detailed_solution": {{
    "consistency_check": "<가설의 기본 주장(예: 더 활성인 화합물)이 데이터와 일치하는지 검토>",
    "aspect_validation": "<기전 분석의 각 항목(전자, 입체 등)이 타당한지 개별 검증 및 논평>",
    "counter_hypothesis_review": "<제시된 반대 가설이 합리적인지, 그리고 반박 논리가 타당한지 검토>",
    "design_suggestion_review": "<설계 제안이 가설 검증에 적합한지, 그리고 그 근거가 합리적인지 검토>",
    "additional_requirements": "<가설을 더 견고하게 만들기 위해 필요한 추가 데이터나 분석 제안>"
  }}
}}
```

### 자기 교정 지침 (Self-Correction Instruction)
제출 전, 응답이 위의 JSON 스키마를 정확히 따르는지, 그리고 `Evaluation.md`의 심사 기준을 만족하는지 최종 점검하세요.
'''

    system_message = system_message_evaluation

    try:
        chat_result = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user", "content": prompt},
            ],
            response_format={"type": "json_object"},
            max_tokens=max_tokens,
            temperature=0.2,
        )
        return {
            "content": chat_result.choices[0].message.content.strip(),
            "usage": _extract_usage(chat_result),
            "model": model,
        }
    except Exception as e:
        return {
            "content": f"LLM 가설 평가 중 오류 발생: {e}",
            "usage": {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0},
            "model": model,
        }


def revise_hypothesis(
    api_key: str,
    original_hypothesis_text: str,
    review_findings: str,
    smiles1: str,
    activity1: float,
    smiles2: str,
    activity2: float,
    structural_difference_description: str,
    model: str = "gpt-4o-mini",
    max_tokens: int = 4096,
) -> Dict[str, Any]:
    """
    평가 결과를 바탕으로 기존 SAR 가설을 수정합니다.

    LLM에 다음을 요청합니다:
    - 제공된 평가 피드백을 반영하여 원본 가설을 개선합니다.
    - 논리적 오류를 바로잡고, 근거를 보강하며, 제안된 사항을 통합합니다.
    - 수정된 가설은 `Generation.md`의 출력 형식과 동일한 엄격한 JSON 형식이어야 합니다.

    Args:
        api_key (str): OpenAI API 키.
        original_hypothesis_text (str): 수정할 원본 가설 (JSON 형식).
        review_findings (str): 가설에 대한 평가 결과 (JSON 형식).
        smiles1 (str): 첫 번째 화합물의 SMILES 문자열.
        activity1 (float): 첫 번째 화합물의 활성도.
        smiles2 (str): 두 번째 화합물의 SMILES 문자열.
        activity2 (float): 두 번째 화합물의 활성도.
        structural_difference_description (str): 두 화합물 간의 구조적 차이 요약.
        model (str): 사용할 OpenAI 모델 이름.
        max_tokens (int): 생성할 최대 토큰 수.

    Returns:
        Dict[str, Any]: {"content": JSON 문자열, "usage": 토큰 사용량 dict, "model": 모델명}
    """
    client = OpenAI(api_key=api_key)

    prompt = f'''### 핵심 지침 (Core Instructions)
* **비판적 통합:** 당신은 동료 심사(peer review) 피드백을 바탕으로 논문을 개정하는 연구자입니다. `Revision.md`의 원칙에 따라, 제공된 '심사 결과'를 비판적으로 수용하여 '원본 가설'을 개선하세요.
* **완전한 재작성:** 단순히 피드백을 짜깁기하지 마세요. 피드백을 완전히 소화하여, `Generation.md`의 엄격한 기준을 만족하는 **완전하고 새로운 가설 JSON 객체**를 처음부터 다시 작성하세요.

### 수정 작업 데이터 (Revision Data)
- **화합물 1:** {smiles1} (활성도: {activity1})
- **화합물 2:** {smiles2} (활성도: {activity2})
- **수정 대상 원본 가설 (JSON):**
```json
{original_hypothesis_text}
```
- **반영할 심사 결과 (JSON):**
```json
{review_findings}
```

### 출력 형식 (Output Format)
**지침: `Generation.md`의 출력 형식과 동일한, 아래 키를 사용하는 JSON 객체만을 반환하세요. 값은 한국어로 작성합니다. 심사 결과를 반영하여 내용을 수정하거나 보강하세요.**

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
    {{
      "design": "<MMP 스캔 기반 설계 제안 1(SMILES)>",
      "expected_effect": "<예상 효과 ↑/↓, 약/중/강>",
      "rationale": "<설계의 근거>",
      "validation_metric": "<검증에 사용할 지표, 예: pIC50>"
    }}
  ],
  "admet_flags": [
      "<예상되는 ADMET 위험 1과 그 구조적 근거>"
  ],
  "assumptions_and_limits": [
    "<분석에 사용된 가정 1>",
    "<분석의 한계 1>"
  ],
  "confidence": 0.9
}}
```

### 자기 교정 지침 (Self-Correction Instruction)
제출 전, 응답이 위의 JSON 스키마를 정확히 따르는지, 그리고 `Generation.md`의 엄격성/완결성 기준을 만족하는지 최종 점검하세요.
'''

    system_message = system_message_revision

    try:
        chat_result = client.chat.completions.create(
            model=model,
            messages=[
                {"role": "system", "content": system_message},
                {"role": "user", "content": prompt},
            ],
            response_format={"type": "json_object"},
            max_tokens=max_tokens,
            temperature=0.2,
        )
        return {
            "content": chat_result.choices[0].message.content.strip(),
            "usage": _extract_usage(chat_result),
            "model": model,
        }
    except Exception as e:
        return {
            "content": f"LLM 가설 수정 중 오류 발생: {e}",
            "usage": {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0},
            "model": model,
        }
