from __future__ import annotations

from langchain_core.tools import tool

from rag.retriever import format_hits, retrieve
from scagent.skills_loader import load_skill_text, skill_catalog_text


@tool
def retrieve_papers(query: str, collection: str = "papers") -> str:
    """Search the local RAG corpus. Default collection is knowledge/papers."""
    hits = retrieve(query, collection=collection)
    if not hits and collection == "papers":
        hits = retrieve(query, collections=["methods", "best_practices"])
    return format_hits(hits)


@tool
def list_analysis_skills() -> str:
    """List project skills under skills/ (existing SciAgent-style SKILL.md files)."""
    return skill_catalog_text()


@tool
def load_skill(name: str) -> str:
    """Load a skill SOP by name, e.g. scanpy-scrna-seq, celltypist-cell-annotation."""
    return load_skill_text(name, include_references=False)


TOOLS = [retrieve_papers, list_analysis_skills, load_skill]
TOOLS_BY_NAME = {t.name: t for t in TOOLS}
