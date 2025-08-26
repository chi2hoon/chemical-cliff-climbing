
from openai import OpenAI
from typing import Optional
import re

def generate_hypothesis(
    api_key: str,
    smiles1: str,
    activity1: float,
    smiles2: str,
    activity2: float,
    structural_difference_description: str,
    similarity: Optional[float] = None,
    model: str = "gpt-5-nano",
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

    tanimoto_text = f"{similarity:.3f}" if similarity is not None else "N/A"

    # LLM에 전달할 프롬프트. 명확성과 가독성을 위해 여러 줄 문자열로 구성합니다.
    prompt = f'''### 핵심 지침 (Core Instructions)
* **엄격성이 최우선:** 당신의 1차 목표는 **완전하고 엄격하게 정당화된** SAR 가설을 제시하는 것입니다. 논리 비약 없이 각 단계의 주장에 **구조적 사실**과 **물리·유기화학적 원리**를 명시적으로 연결하세요. 단순한 직감이나 일반론적 상투구는 실패로 간주됩니다.
* **완결성에 대한 정직성:** 완전한 가설(충분히 검증 가능한 일관된 설명)을 제시할 수 없으면 **추측으로 채우지 마십시오**. 그 대신, **유의미한 부분 결과**만 제시하세요.

### 문제 데이터 (Problem Data)
- 화합물 1 SMILES: {smiles1}
- 화합물 2 SMILES: {smiles2}
- 활성도 지표 및 값: = {activity1} (화합물 1), {activity2}(화합물 2)
- 유사도: Tanimoto = {tanimoto_text}
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
        # gpt-5-nano는 JSON 강제 옵션과 매개변수 제약이 있으므로, 별도의 컴팩트 프롬프트/토큰 제한을 적용합니다.
        is_nano = ("gpt-5-nano" in model)
        compact_prompt = None
        if is_nano:
            compact_prompt = (
                "아래 키만 포함한 단일 JSON 객체로만 답변하세요. 백틱/설명/주석 금지. "
                "모호하면 'N/A' 또는 빈 배열을 사용하세요. 반드시 '{'로 시작하고 '}'로 끝내세요. "
                "내부 추론은 최소화하고, 체인-오브-소트는 절대 노출하지 마세요. 출력은 오직 JSON 하나입니다.\n"
                "문서:\n"
                "{"
                "\"more_active\": \"화합물 1 또는 화합물 2\","
                "\"delta_pAct_explained\": \"N/A 또는 간단 설명\","
                "\"primary_hypothesis\": \"1~2문장\","
                "\"mechanistic_rationale\": {"
                "\"electronic\": \"N/A\", \"steric\": \"N/A\", \"hydrogen_bonding\": \"N/A\","
                "\"hydrophobic_polar\": \"N/A\", \"ionization_pKa\": \"N/A\", \"conformation\": \"N/A\", \"aromatic_pi\": \"N/A\"},"
                "\"counter_hypotheses\": [], \"design_suggestions\": [], \"admet_flags\": [],"
                "\"assumptions_and_limits\": [], \"confidence\": 0.7" "}"
                "\n데이터: "
                f"smiles1={smiles1}, activity1={activity1}, smiles2={smiles2}, activity2={activity2}, tanimoto={tanimoto_text}, diff={structural_difference_description}"
            )

        params = {
            "model": model,
            "messages": [
                {"role": "system", "content": system_message},
                {"role": "user", "content": (compact_prompt if compact_prompt else prompt)},
            ],
        }
        # GPT‑5 계열(특히 nano)은 내부 reasoning 토큰이 "max_completion_tokens" 예산을 함께 소모합니다.
        # nano에서는 reasoning.effort를 최소로 두고, 토큰 예산을 넉넉히 잡아 JSON이 끊기지 않게 합니다.
        if is_nano:
            params["max_completion_tokens"] = min(max_tokens, 4096)
        else:
            # 비‑reasoning 모델은 기존 "max_tokens"를 사용
            params["max_tokens"] = max_tokens
        # 일부 경량 모델(gpt-5-nano)은 response_format을 지원하지 않거나 비표준 응답을 낼 수 있음
        if not is_nano:
            params["response_format"] = {"type": "json_object"}
        # gpt-5-nano는 temperature 커스터마이즈 미지원 → 파라미터 제거
        if not is_nano:
            params["temperature"] = 0.7

        chat_result = client.chat.completions.create(**params)
        raw = chat_result.choices[0].message.content if chat_result and chat_result.choices else ""
        # 우선 JSON 추출 시도, 실패 시 원문(raw) 그대로 반환하여 UI에서 표시되게 함
        extracted = _extract_json_string(raw)
        if _is_json_object_string(extracted):
            return extracted.strip()
        # synthesize minimal JSON to avoid blank UI
        return _coerce_hypothesis_json(
            raw_text=(raw or "").strip(),
            smiles1=smiles1,
            activity1=activity1,
            smiles2=smiles2,
            activity2=activity2,
            similarity=tanimoto_text,
        )
    except Exception as e:
        return f"LLM 가설 생성 중 오류 발생: {e}"


def evaluate_hypothesis(
    api_key: str,
    hypothesis_text: str,
    smiles1: str,
    activity1: float,
    smiles2: str,
    activity2: float,
    structural_difference_description: str = "",
    model: str = "gpt-4-turbo",
    max_tokens: int = 4096,
) -> str:
    """
    기존 가설을 엄격성/일관성 기준으로 평가하여 JSON 스키마로 반환합니다.

    반환 스키마:
    {
      "summary": {"verdict": str, "method_sketch": str},
      "detailed_solution": {
        "consistency_check": str,
        "aspect_validation": str,
        "counter_hypothesis_review": str,
        "design_suggestion_review": str,
        "additional_requirements": str
      }
    }
    """
    client = OpenAI(api_key=api_key)

    prompt = f'''### 평가 대상 데이터
- 화합물 1: {smiles1} | activity = {activity1}
- 화합물 2: {smiles2} | activity = {activity2}
- 구조 차이 요약: {structural_difference_description}

### 평가 대상 가설 본문
{hypothesis_text}

### 과제
위 가설을 다음 기준으로 엄정 평가하고, 아래 JSON 스키마에 맞춰 결과만 반환하세요.
- 기본 일치성: 데이터/가정/논리의 일관성, pAct/단위 가정 명확성
- 관점별 검증: 전자/입체/HBD-HBA/소수성-극성/pKa/컨포메이션/π 상호작용의 적절성
- 반대 가설 검토: 적절한 대안 및 반례 조건 제시 여부
- 설계 제안 검토: 제안의 실험 가능성과 검증 지표의 타당성
- 추가 필요 요소: 데이터/분석 보강 필요 사항

JSON 스키마:
{{
  "summary": {{"verdict": "승인/보류/거부", "method_sketch": "평가 접근 요약"}},
  "detailed_solution": {{
    "consistency_check": "문장",
    "aspect_validation": "문장",
    "counter_hypothesis_review": "문장",
    "design_suggestion_review": "문장",
    "additional_requirements": "문장"
  }}
}}
'''

    try:
        is_nano = ("gpt-5-nano" in model)
        params = {
            "model": model,
            "messages": [
                {"role": "system", "content": "You are a rigorous medicinal chemistry reviewer."},
                {"role": "user", "content": prompt},
            ],
        }
        if is_nano:
            params["max_completion_tokens"] = min(max_tokens, 4096)
        else:
            params["max_tokens"] = max_tokens
            params["response_format"] = {"type": "json_object"}
            params["temperature"] = 0.3

        chat_result = client.chat.completions.create(**params)
        raw = chat_result.choices[0].message.content if chat_result and chat_result.choices else ""
        extracted = _extract_json_string(raw)
        if _is_json_object_string(extracted):
            return extracted.strip()
        return _coerce_evaluation_json((raw or "").strip())
    except Exception as e:
        return f"LLM 가설 평가 중 오류 발생: {e}"


def revise_hypothesis(
    api_key: str,
    original_hypothesis_text: str,
    review_findings: str,
    smiles1: str,
    activity1: float,
    smiles2: str,
    activity2: float,
    structural_difference_description: str = "",
    model: str = "gpt-4-turbo",
    max_tokens: int = 4096,
) -> str:
    """
    평가 결과를 반영해 가설을 수정하여, 생성 가설과 동일한 JSON 스키마로 반환합니다.
    """
    client = OpenAI(api_key=api_key)

    prompt = f'''### 문제 데이터
- 화합물 1: {smiles1} | activity = {activity1}
- 화합물 2: {smiles2} | activity = {activity2}
- 구조 차이 요약: {structural_difference_description}

### 원본 가설
{original_hypothesis_text}

### 리뷰 결과 요약(JSON 또는 텍스트)
{review_findings}

### 과제
위 리뷰의 지적 사항을 모두 반영하여 가설을 수정하세요. 출력은 "생성 가설"과 동일한 JSON 스키마여야 합니다.
필수 키: more_active, delta_pAct_explained, primary_hypothesis, mechanistic_rationale(7개 키), counter_hypotheses(2+), design_suggestions(3개), admet_flags(2+), assumptions_and_limits(2+), confidence.
'''

    try:
        is_nano = ("gpt-5-nano" in model)
        params = {
            "model": model,
            "messages": [
                {"role": "system", "content": "You are a medicinal chemistry assistant that writes concise, well-justified SAR hypotheses."},
                {"role": "user", "content": prompt},
            ],
        }
        if is_nano:
            params["max_completion_tokens"] = min(max_tokens, 4096)
        else:
            params["max_tokens"] = max_tokens
            params["response_format"] = {"type": "json_object"}
            params["temperature"] = 0.7

        chat_result = client.chat.completions.create(**params)
        raw = chat_result.choices[0].message.content if chat_result and chat_result.choices else ""
        extracted = _extract_json_string(raw)
        if _is_json_object_string(extracted):
            return extracted.strip()
        return _coerce_hypothesis_json(
            raw_text=(raw or "").strip(),
            smiles1=smiles1,
            activity1=activity1,
            smiles2=smiles2,
            activity2=activity2,
            similarity="N/A",
        )
    except Exception as e:
        return f"LLM 가설 수정 중 오류 발생: {e}"


# --- Helpers ---
def _extract_json_string(text: str) -> str:
    """Best-effort extraction of a JSON object string from arbitrary text.

    Strategy:
    1) If fenced code (```json ... ``` or ``` ... ```) exists, prefer it.
    2) Otherwise, locate the first top-level JSON object via brace matching.
    3) Fallback to original text (strip), which may still be valid JSON.
    """
    if not text:
        return ""

    s = text.strip()

    # 1) Fenced code block with optional json language tag
    match = re.search(r"```(?:json)?\s*([\s\S]*?)```", s, re.IGNORECASE)
    if match:
        candidate = match.group(1).strip()
        if candidate.startswith("{") and candidate.endswith("}"):
            return candidate
        # Try to find JSON object inside candidate
        inner = _find_first_json_object(candidate)
        if inner:
            return inner

    # 2) Find first top-level JSON object in the string
    inner = _find_first_json_object(s)
    if inner:
        return inner

    # 3) Fallback to the original (trimmed) text
    return s


def _find_first_json_object(text: str) -> str:
    in_string = False
    escape = False
    depth = 0
    start_index = None
    for i, ch in enumerate(text):
        if ch == "\\" and not escape:
            escape = True
            continue
        if ch == '"' and not escape:
            in_string = not in_string
        if not in_string:
            if ch == '{':
                if depth == 0:
                    start_index = i
                depth += 1
            elif ch == '}':
                if depth > 0:
                    depth -= 1
                    if depth == 0 and start_index is not None:
                        return text[start_index:i+1].strip()
        if escape:
            escape = False
    return ""


def _is_json_object_string(s: str) -> bool:
    if not s:
        return False
    t = s.strip()
    return t.startswith("{") and t.endswith("}")


def _coerce_hypothesis_json(raw_text: str, smiles1: str, activity1: float, smiles2: str, activity2: float, similarity: str) -> str:
    # Minimal, schema-compatible JSON so UI renders something useful
    safe_text = (raw_text or "")[:4000]
    return (
        '{"more_active": "N/A", "delta_pAct_explained": "N/A", '
        f'"primary_hypothesis": {json_dumps_safe(safe_text)}, '
        '"mechanistic_rationale": {"electronic": "N/A", "steric": "N/A", "hydrogen_bonding": "N/A", '
        '"hydrophobic_polar": "N/A", "ionization_pKa": "N/A", "conformation": "N/A", "aromatic_pi": "N/A"}, '
        '"counter_hypotheses": [], "design_suggestions": [], "admet_flags": [], '
        '"assumptions_and_limits": [], "confidence": 0.5}'
    )


def _coerce_evaluation_json(raw_text: str) -> str:
    safe_text = (raw_text or "")[:4000]
    return (
        '{"summary": {"verdict": "보류", "method_sketch": "자동 폴백"}, '
        f'"detailed_solution": {{"consistency_check": {json_dumps_safe(safe_text)}, '
        '"aspect_validation": "N/A", "counter_hypothesis_review": "N/A", '
        '"design_suggestion_review": "N/A", "additional_requirements": "N/A"}}}'
    )


def json_dumps_safe(text: str) -> str:
    # very small, dependency-free string escaper for embedding into JSON
    if text is None:
        return '""'
    escaped = (
        text.replace("\\", "\\\\")
        .replace("\"", "\\\"")
        .replace("\n", "\\n")
        .replace("\r", "\\r")
        .replace("\t", "\\t")
    )
    return '"' + escaped + '"'
