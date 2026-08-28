from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from scagent.config import REPO_ROOT

REFERENCE_DIR = REPO_ROOT / "knowledge" / "best_practices"

# Analysis step → reference filenames (stem, no .md)
_PHASE_DOCS: dict[str, tuple[str, ...]] = {
    "qc": ("qc", "doublet-detection", "normalization", "feature-selection"),
    "clustering": ("clustering", "dimensionality-reduction"),
    "annotation": ("cell-annotation", "marker-genes"),
    "integration": ("integration",),
    "deg": ("pseudobulk-de", "marker-genes"),
    "enrichment": ("pathway-enrichment",),
    "trajectory": ("trajectory",),
    "downstream": (
        "clustering",
        "dimensionality-reduction",
        "cell-annotation",
        "marker-genes",
        "integration",
        "pseudobulk-de",
        "trajectory",
    ),
}

_ROUTE_DOCS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("qc", ("qc", "doublet-detection", "normalization", "feature-selection")),
    ("normalize", ("normalization", "feature-selection")),
    ("hvg", ("feature-selection",)),
    ("pca", ("dimensionality-reduction",)),
    ("neighbors", ("dimensionality-reduction", "clustering")),
    ("leiden", ("clustering",)),
    ("umap", ("dimensionality-reduction",)),
    ("harmony", ("integration",)),
    ("scvi", ("integration",)),
    ("integrate", ("integration",)),
    ("annotate", ("cell-annotation", "marker-genes")),
    ("annotation", ("cell-annotation", "marker-genes")),
    ("deg", ("pseudobulk-de", "marker-genes")),
    ("gsea", ("pathway-enrichment",)),
    ("enrichment", ("pathway-enrichment",)),
    ("trajectory", ("trajectory",)),
    ("paga", ("trajectory",)),
    ("dpt", ("trajectory",)),
)


@dataclass(frozen=True)
class Practice:
    name: str
    title: str
    path: Path
    body: str


def list_practices() -> list[Practice]:
    if not REFERENCE_DIR.exists():
        return []
    out: list[Practice] = []
    for path in sorted(REFERENCE_DIR.glob("*.md")):
        if path.name.lower() == "readme.md":
            continue
        text = path.read_text(encoding="utf-8")
        title = path.stem
        for line in text.splitlines():
            if line.startswith("# "):
                title = line[2:].strip()
                break
        out.append(Practice(name=path.stem, title=title, path=path, body=text))
    return out


def get_practice(name: str) -> Practice | None:
    key = name.strip().lower().removesuffix(".md")
    for p in list_practices():
        if p.name.lower() == key:
            return p
    return None


def load_practice_text(name: str, *, max_chars: int = 4000) -> str:
    p = get_practice(name)
    if p is None:
        available = ", ".join(x.name for x in list_practices()) or "none"
        return f"Best-practice '{name}' not found. Available: {available}"
    body = p.body.strip()
    if len(body) > max_chars:
        body = body[: max_chars - 1].rstrip() + "…"
    return f"# {p.name}\n\n{body}"


def practices_for_route(
    route: list[str] | None = None,
    intents: list[str] | None = None,
    query: str | None = None,
) -> list[str]:
    selected: list[str] = []

    def add(*names: str) -> None:
        for n in names:
            if n not in selected and get_practice(n):
                selected.append(n)

    steps = list(route or []) + list(intents or [])
    q = (query or "").lower()
    add("qc", "doublet-detection")
    for step in steps:
        key = str(step).lower()
        for needle, docs in _ROUTE_DOCS:
            if needle in key:
                add(*docs)
    if any(k in q for k in ("轨迹", "拟时序", "pseudotime", "trajectory", "velocity", "paga")):
        add("trajectory")
    if any(k in q for k in ("差异", "deg", "deseq", "edger", "pseudobulk", "条件")):
        add("pseudobulk-de")
    if any(k in q for k in ("通路", "gsea", "enrich", "hallmark")):
        add("pathway-enrichment")
    if any(k in q for k in ("整合", "batch", "harmony", "scvi", "批次")):
        add("integration")
    return selected


def practices_for_phase(phase: str, *, route: list[str] | None = None, query: str | None = None) -> list[str]:
    if phase == "qc":
        names = list(_PHASE_DOCS["qc"])
    elif phase in {"interpret", "enrichment"}:
        names = list(_PHASE_DOCS["enrichment"]) + list(_PHASE_DOCS["deg"])
    else:
        names = practices_for_route(route, query=query)
        for must in ("clustering", "cell-annotation", "marker-genes"):
            if must not in names:
                names.append(must)
    seen: list[str] = []
    for n in names:
        if n not in seen and get_practice(n):
            seen.append(n)
    return seen


def practices_context(
    phase: str,
    *,
    route: list[str] | None = None,
    query: str | None = None,
    max_chars: int = 8000,
) -> str:
    names = practices_for_phase(phase, route=route, query=query)
    if not names:
        names = [p.name for p in list_practices()]
    parts: list[str] = []
    used = 0
    per = max(1200, max_chars // max(len(names), 1))
    for name in names:
        text = load_practice_text(name, max_chars=per)
        if used + len(text) > max_chars:
            remain = max_chars - used
            if remain > 400:
                parts.append(text[:remain].rstrip() + "…")
            break
        parts.append(text)
        used += len(text)
    return "\n\n".join(parts)


def practices_catalog_text() -> str:
    rows = list_practices()
    if not rows:
        return "(no knowledge/best_practices docs)"
    lines = [f"knowledge/best_practices: {len(rows)} step SOPs (Heumos 2023 / sc-best-practices)"]
    for p in rows:
        lines.append(f"- {p.name}: {p.title}")
    return "\n".join(lines)
