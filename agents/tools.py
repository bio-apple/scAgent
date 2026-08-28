from __future__ import annotations

from langchain_core.tools import tool

from rag.retriever import format_hits, retrieve, retrieve_fused
from scagent.best_practices_loader import load_practice_text, practices_catalog_text
from scagent.kb import format_records, lookup_structured
from scagent.skills_loader import load_skill_text, skill_catalog_text


@tool
def retrieve_papers(query: str, collection: str = "") -> str:
    """Search the fused RAG corpus (best_practices + papers + methods + sops + upstream)."""
    if collection:
        hits = retrieve(query, collection=collection)
        if not hits:
            hits = retrieve_fused(query)
    else:
        hits = retrieve_fused(query)
    return format_hits(hits)


@tool
def list_analysis_skills() -> str:
    """List project skills under skills/ (existing SciAgent-style SKILL.md files)."""
    return skill_catalog_text()


@tool
def load_skill(name: str) -> str:
    """Load a skill SOP by name, e.g. scanpy-scrna-seq, celltypist-cell-annotation."""
    return load_skill_text(name, include_references=False)


@tool
def list_best_practices() -> str:
    """List step SOPs under knowledge/best_practices (Heumos 2023 / sc-best-practices)."""
    return practices_catalog_text()


@tool
def load_best_practice(name: str) -> str:
    """Load a best-practice SOP by stem, e.g. qc, clustering, pseudobulk-de, integration."""
    return load_practice_text(name)


@tool
def lookup_knowledge(query: str, collection: str = "") -> str:
    """Lookup structured KB records (cell ontology, markers, pathways, disease signatures, tissues)."""
    cols = [collection] if collection.strip() else None
    return format_records(lookup_structured(query, collections=cols))


TOOLS = [retrieve_papers, lookup_knowledge, list_analysis_skills, load_skill, list_best_practices, load_best_practice]
TOOLS_BY_NAME = {t.name: t for t in TOOLS}
