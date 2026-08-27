from __future__ import annotations

import os
import random
import time
from typing import Any

from scagent.config import load_config
from scagent.cache import llm_key, load_json, save_json
from scagent.logutil import get_logger, timed


def _prompts():
    from scagent.config import REPO_ROOT

    return REPO_ROOT / "prompts"


def read_prompt(name: str) -> str:
    return (_prompts() / f"{name}.md").read_text(encoding="utf-8")


log = get_logger("llm")
_last_call = 0.0
_tokens = {"input": 0, "output": 0, "total": 0}


def token_usage() -> dict[str, int]:
    return dict(_tokens)


def get_llm(cfg: dict | None = None):
    cfg = cfg or load_config()
    model_cfg = cfg["model"]
    key_env = model_cfg.get("api_key_env") or "OPENAI_API_KEY"
    api_key = os.getenv(key_env) or os.getenv("OPENAI_API_KEY")
    if not api_key:
        return None
    from langchain_openai import ChatOpenAI

    kwargs: dict[str, Any] = {
        "model": os.getenv("OPENAI_MODEL") or model_cfg["name"],
        "api_key": api_key,
        "temperature": float(model_cfg.get("temperature") or 0),
        "max_retries": 0,
    }
    base_url = os.getenv("OPENAI_BASE_URL") or model_cfg.get("base_url")
    if base_url:
        kwargs["base_url"] = base_url
    max_tokens = model_cfg.get("max_tokens")
    if max_tokens:
        kwargs["max_tokens"] = int(max_tokens)
    return ChatOpenAI(**kwargs)


def _retryable(exc: BaseException) -> bool:
    name = type(exc).__name__.lower() + str(exc).lower()
    keys = ("ratelimit", "rate limit", "timeout", "connection", "429", "503", "502", "overloaded", "unavailable")
    return any(k in name for k in keys)


def _rate_limit(cfg: dict) -> None:
    global _last_call
    rpm = float((cfg.get("model") or {}).get("rate_limit_rpm") or 0)
    if rpm <= 0:
        return
    wait = 60.0 / rpm
    gap = time.monotonic() - _last_call
    if gap < wait:
        time.sleep(wait - gap)
    _last_call = time.monotonic()


def _log_tokens(ai: Any) -> None:
    meta = getattr(ai, "usage_metadata", None) or {}
    inp = int(meta.get("input_tokens") or meta.get("prompt_tokens") or 0)
    out = int(meta.get("output_tokens") or meta.get("completion_tokens") or 0)
    tot = int(meta.get("total_tokens") or (inp + out))
    _tokens["input"] += inp
    _tokens["output"] += out
    _tokens["total"] += tot
    if tot:
        log.info("tokens this call in=%s out=%s total=%s | session total=%s", inp, out, tot, _tokens["total"])


def invoke_llm(model, messages, cfg: dict | None = None):
    cfg = cfg or load_config()
    mcfg = cfg.get("model") or {}
    use_cache = True
    if cfg.get("performance") is not None:
        use_cache = bool((cfg.get("performance") or {}).get("cache", True))
    else:
        from scagent.cache import cache_enabled

        use_cache = cache_enabled()
    ck = llm_key(messages)
    if use_cache:
        cached = load_json(ck)
        if isinstance(cached, dict) and "content" in cached:

            class _Hit:
                content = cached["content"]
                usage_metadata = {}
                tool_calls = []

            log.info("LLM cache hit %s", ck)
            return _Hit()
    retries = int(mcfg.get("max_retries") or 4)
    base = float(mcfg.get("retry_backoff_seconds") or 1.0)
    cap = float(mcfg.get("retry_backoff_max") or 30.0)
    last: BaseException | None = None
    for attempt in range(retries + 1):
        try:
            _rate_limit(cfg)
            with timed("llm.invoke", log):
                ai = model.invoke(messages)
            _log_tokens(ai)
            calls = getattr(ai, "tool_calls", None) or []
            content = ai.content
            if use_cache and not calls and isinstance(content, str):
                save_json(ck, {"content": content})
            return ai
        except Exception as exc:
            last = exc
            if attempt >= retries or not _retryable(exc):
                log.error("LLM invoke failed: %s", exc)
                raise
            sleep = min(cap, base * (2**attempt)) + random.uniform(0, 0.3)
            log.warning("LLM retry %s/%s in %.1fs (%s)", attempt + 1, retries, sleep, type(exc).__name__)
            time.sleep(sleep)
    raise last  # pragma: no cover


def invoke_json(system_prompt: str, user_text: str, cfg: dict | None = None) -> dict | None:
    """LLM JSON object via response_format. None if no key or parse failure."""
    import json
    import re

    llm = get_llm(cfg)
    if llm is None:
        return None
    from langchain_core.messages import HumanMessage, SystemMessage

    bound = llm
    try:
        bound = llm.bind(response_format={"type": "json_object"})
    except Exception:
        bound = llm
    ai = invoke_llm(bound, [SystemMessage(content=system_prompt), HumanMessage(content=user_text)], cfg)
    content = ai.content if isinstance(ai.content, str) else str(ai.content)
    try:
        data = json.loads(content)
        return data if isinstance(data, dict) else None
    except json.JSONDecodeError:
        m = re.search(r"\{.*\}", content, re.S)
        if not m:
            return None
        try:
            data = json.loads(m.group())
            return data if isinstance(data, dict) else None
        except json.JSONDecodeError:
            return None


def run_specialist(system_prompt: str, user_text: str, cfg: dict | None = None) -> str | None:
    """Bounded tool loop. Returns None if no API key (caller uses fallback)."""
    llm = get_llm(cfg)
    if llm is None:
        return None
    from langchain_core.messages import HumanMessage, SystemMessage, ToolMessage

    from agents.tools import TOOLS, TOOLS_BY_NAME

    cfg = cfg or load_config()
    max_rounds = int(cfg["model"].get("max_tool_rounds") or 4)
    model = llm.bind_tools(TOOLS)
    messages: list = [SystemMessage(content=system_prompt), HumanMessage(content=user_text)]
    for _ in range(max_rounds):
        ai = invoke_llm(model, messages, cfg)
        messages.append(ai)
        calls = getattr(ai, "tool_calls", None) or []
        if not calls:
            content = ai.content
            return content if isinstance(content, str) else str(content)
        for call in calls:
            if isinstance(call, dict):
                name, args, tid = call["name"], call.get("args") or {}, call["id"]
            else:
                name, args, tid = call.name, call.args or {}, call.id
            tool = TOOLS_BY_NAME.get(name)
            observation = tool.invoke(args) if tool else f"unknown tool: {name}"
            messages.append(ToolMessage(content=str(observation), tool_call_id=tid))
    last = messages[-1]
    content = getattr(last, "content", last)
    return content if isinstance(content, str) else str(content)
