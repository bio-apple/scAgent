from __future__ import annotations

import os
from typing import Any

from scagent.config import REPO_ROOT, load_config

PROMPTS = REPO_ROOT / "prompts"


def read_prompt(name: str) -> str:
    path = PROMPTS / f"{name}.md"
    return path.read_text(encoding="utf-8")


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
    }
    base_url = os.getenv("OPENAI_BASE_URL") or model_cfg.get("base_url")
    if base_url:
        kwargs["base_url"] = base_url
    return ChatOpenAI(**kwargs)


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
        ai = model.invoke(messages)
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
