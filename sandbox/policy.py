"""Static gate for LLM-generated scripts. First line of defense, not a complete sandbox."""

from __future__ import annotations

import re

# Obvious host-escape / destructive APIs. Generated Scanpy templates do not need these.
_BLOCK = (
    r"\bos\.system\b",
    r"\bos\.popen\b",
    r"\bos\.execv",
    r"\bos\.execl",
    r"\bos\.spawn",
    r"\bos\.fork\b",
    r"\bsubprocess\b",
    r"\bshutil\.rmtree\b",
    r"\bctypes\b",
    r"\bpty\b",
    r"\bsocket\.socket\b",
    r"\beval\s*\(",
    r"\bexec\s*\(",
    r"\bcompile\s*\(",
    r"__import__\s*\(",
    r"\bpickle\.loads\b",
    r"\bmarshal\.loads\b",
    r"\bwebbrowser\b",
    r"\bpty\.spawn\b",
    r"\bimportlib\.import_module\b",
)

_COMPILED = [re.compile(p) for p in _BLOCK]


def policy_violations(code: str) -> list[str]:
    hits: list[str] = []
    for rx in _COMPILED:
        m = rx.search(code or "")
        if m:
            hits.append(m.group(0))
    return hits
