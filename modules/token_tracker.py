from __future__ import annotations

from typing import Any, Dict, List, Optional
from dataclasses import dataclass, field
import datetime as _dt


def _extract_usage_from_response(resp: Any) -> Dict[str, int]:
    """Best-effort extraction of token usage from an OpenAI response object.

    Supports classic Chat Completions usage fields (prompt_tokens, completion_tokens,
    total_tokens) and newer input/output naming. Missing values default to 0.
    """
    usage = getattr(resp, "usage", None)
    result = {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0}

    if not usage:
        return result

    # Try attribute access
    for k in ("prompt_tokens", "completion_tokens", "total_tokens"):
        v = getattr(usage, k, None)
        if isinstance(v, int):
            result[k] = v

    # If attribute access didn't work fully, try dict-like access/model_dump
    usage_dict = None
    if isinstance(usage, dict):
        usage_dict = usage
    else:
        to_dict = getattr(usage, "to_dict", None)
        model_dump = getattr(usage, "model_dump", None)
        if callable(model_dump):
            try:
                usage_dict = model_dump()
            except Exception:
                usage_dict = None
        if usage_dict is None and callable(to_dict):
            try:
                usage_dict = to_dict()
            except Exception:
                usage_dict = None

    if isinstance(usage_dict, dict):
        result["prompt_tokens"] = int(
            usage_dict.get("prompt_tokens", usage_dict.get("input_tokens", result["prompt_tokens"])) or 0
        )
        result["completion_tokens"] = int(
            usage_dict.get("completion_tokens", usage_dict.get("output_tokens", result["completion_tokens"])) or 0
        )
        total = usage_dict.get("total_tokens")
        result["total_tokens"] = int(total) if isinstance(total, int) else (result["prompt_tokens"] + result["completion_tokens"]) 

    if result["total_tokens"] < (result["prompt_tokens"] + result["completion_tokens"]):
        result["total_tokens"] = result["prompt_tokens"] + result["completion_tokens"]

    return result


@dataclass
class TokenUsageRecord:
    phase: Optional[str]
    model: Optional[str]
    prompt_tokens: int
    completion_tokens: int
    total_tokens: int
    at: _dt.datetime = field(default_factory=lambda: _dt.datetime.utcnow())


class TokenTracker:
    """Simple in-memory tracker to accumulate token usage across calls."""

    def __init__(self) -> None:
        self._calls: List[TokenUsageRecord] = []
        self._totals: Dict[str, int] = {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0}
        self._per_model: Dict[str, Dict[str, int]] = {}

    def reset(self) -> None:
        self._calls.clear()
        self._totals = {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0}
        self._per_model.clear()

    def add(self, prompt_tokens: int, completion_tokens: int, total_tokens: Optional[int] = None,
            *, model: Optional[str] = None, phase: Optional[str] = None) -> None:
        pt = int(prompt_tokens or 0)
        ct = int(completion_tokens or 0)
        tt = int(total_tokens if total_tokens is not None else (pt + ct))

        self._calls.append(TokenUsageRecord(phase=phase, model=model, prompt_tokens=pt, completion_tokens=ct, total_tokens=tt))
        self._totals["prompt_tokens"] += pt
        self._totals["completion_tokens"] += ct
        self._totals["total_tokens"] += tt

        if model:
            m = self._per_model.setdefault(model, {"prompt_tokens": 0, "completion_tokens": 0, "total_tokens": 0})
            m["prompt_tokens"] += pt
            m["completion_tokens"] += ct
            m["total_tokens"] += tt

    def add_from_response(self, resp: Any, *, model: Optional[str] = None, phase: Optional[str] = None) -> Dict[str, int]:
        usage = _extract_usage_from_response(resp)
        self.add(usage.get("prompt_tokens", 0), usage.get("completion_tokens", 0), usage.get("total_tokens"), model=model, phase=phase)
        return usage

    def totals(self) -> Dict[str, int]:
        return dict(self._totals)

    def per_model(self) -> Dict[str, Dict[str, int]]:
        return {k: dict(v) for k, v in self._per_model.items()}

    def calls(self) -> List[TokenUsageRecord]:
        return list(self._calls)


# Convenience singleton
default_tracker = TokenTracker()

