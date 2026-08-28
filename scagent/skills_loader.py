"""Discover and recommend bundled analysis skills under skills/*/SKILL.md."""

from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from scagent.config import REPO_ROOT, load_config, resolve_path

_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---\s*\n(.*)$", re.DOTALL)
_MANIFEST_PATH = REPO_ROOT / "skills" / "awesome_single_cell_manifest.json"
_TOKEN = re.compile(r"[a-z][a-z0-9_-]{2,}")
_STOP = {
    "the",
    "and",
    "for",
    "with",
    "from",
    "this",
    "that",
    "use",
    "when",
    "data",
    "using",
    "analysis",
    "skill",
    "workflow",
}

# Lean core for a standard scRNA run; topic skills (CellChat, census, …) added on demand.
_CORE_SKILLS = (
    "anndata-data-structure",
    "scanpy-scrna-seq",
    "harmony-batch-correction",
    "scvi-tools-single-cell",
    "single-cell-annotation-guide",
    "single-cell-annotation",
    "celltypist-cell-annotation",
)

# Archived packs live under skills/_archive/ and are not discovered.
_ARCHIVE_DIRNAME = "_archive"

# Planner prompt: recommended first, then fill up to this many lines.
_DEFAULT_CATALOG_LIMIT = 40

_CATEGORIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("trajectory", ("trajectory", "pseudotime", "velocity", "lineage", "monocle", "palantir", "拟时序", "轨迹")),
    ("communication", ("communication", "cellchat", "nichenet", "liana", "ligand", "cellphonedb", "通讯", "配体")),
    ("annotation", ("annotat", "celltypist", "azimuth", "singler", "popv", "cell-type", "注释")),
    ("integration", ("batch", "harmony", "scvi", "integrat", "整合", "批次")),
    ("grn", ("grn", "scenic", "regulon", "arboreto", "gene-regulatory", "调控")),
    ("qc", ("qc", "preprocess", "doublet", "normaliz", "质控", "预处理", "sparse", "io")),
    ("clustering", ("cluster", "leiden", "pca", "umap", "聚类")),
    ("atlas", ("atlas", "census", "mapping", "lamindb")),
)

_TOPIC_KEYWORDS: tuple[tuple[tuple[str, ...], tuple[str, ...]], ...] = tuple(
    (needles, needles) for _, needles in _CATEGORIES
)


@dataclass(frozen=True)
class Skill:
    name: str
    description: str
    path: Path
    body: str
    references: tuple[Path, ...]


def _strip_leading_noise(text: str) -> str:
    """Drop HTML comments / BOM so YAML frontmatter can be parsed."""
    text = (text or "").lstrip("\ufeff")
    while True:
        s = text.lstrip()
        if s.startswith("<!--"):
            end = s.find("-->")
            if end < 0:
                break
            text = s[end + 3 :]
            continue
        break
    return text.lstrip()


def _parse_frontmatter(text: str) -> tuple[dict[str, str], str]:
    text = _strip_leading_noise(text)
    m = _FRONTMATTER.match(text)
    if not m:
        return {}, text
    meta: dict[str, str] = {}
    for line in m.group(1).splitlines():
        if ":" not in line:
            continue
        k, v = line.split(":", 1)
        meta[k.strip()] = v.strip().strip('"').strip("'")
    return meta, m.group(2).strip()


def skills_dir(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    return resolve_path(cfg, "skills")


def list_skills(root: Path | None = None) -> list[Skill]:
    base = root or skills_dir()
    skills: list[Skill] = []
    if not base.exists():
        return skills
    for skill_md in sorted(base.glob("*/SKILL.md")):
        if skill_md.parent.name.startswith("_"):
            continue
        text = skill_md.read_text(encoding="utf-8")
        meta, body = _parse_frontmatter(text)
        name = meta.get("name") or skill_md.parent.name
        refs = tuple(sorted(skill_md.parent.joinpath("references").glob("*.md")))
        skills.append(
            Skill(
                name=name,
                description=meta.get("description", ""),
                path=skill_md,
                body=body,
                references=refs,
            )
        )
    return skills


def get_skill(name: str, root: Path | None = None) -> Skill | None:
    key = name.strip().lower()
    for skill in list_skills(root):
        if skill.name.lower() == key or skill.path.parent.name.lower() == key:
            return skill
    return None


def skill_category(skill: Skill) -> str:
    hay = f"{skill.name} {skill.description}".lower()
    for cat, needles in _CATEGORIES:
        if any(n.lower() in hay for n in needles):
            return cat
    return "other"


def _truncate(text: str, limit: int = 140) -> str:
    text = " ".join(text.split())
    if len(text) <= limit:
        return text
    return text[: limit - 1].rstrip() + "…"


def _topic_query(metadata: dict | None, query: str | None) -> str:
    meta = metadata or {}
    parts = [
        str(meta.get("task") or ""),
        str(meta.get("query") or ""),
        str(meta.get("user_query") or ""),
        str(query or ""),
        str(meta.get("tissue") or ""),
        str(meta.get("platform") or ""),
        str(meta.get("integrator") or ""),
        " ".join(str(x) for x in (meta.get("intents") or [])),
        " ".join(str(x) for x in (meta.get("route") or [])),
    ]
    return " ".join(parts).lower()


def _tokens(text: str) -> set[str]:
    return {t for t in _TOKEN.findall(text.lower()) if t not in _STOP}


def _skill_score(skill: Skill, topic: str) -> int:
    hay = f"{skill.name} {skill.description}".lower()
    score = 0
    if skill.name in _CORE_SKILLS:
        score += 8
    for triggers, needles in _TOPIC_KEYWORDS:
        if any(t in topic for t in triggers) and any(n in hay for n in needles):
            score += 6
            break
    overlap = _tokens(topic) & _tokens(hay)
    score += min(5, len(overlap))
    for word in _tokens(topic):
        if len(word) >= 5 and word in hay:
            score += 1
    return score


def recommend_skills(metadata: dict, language: str = "python") -> list[str]:
    """Score bundled skills against the task; always keep lean SciAgent core."""
    skills = list_skills()
    if not skills:
        return []
    names = {s.name for s in skills}
    selected: list[str] = []

    def add(*cands: str) -> None:
        for c in cands:
            if c in names and c not in selected:
                selected.append(c)

    add(*_CORE_SKILLS)
    tissue = str(metadata.get("tissue") or "").lower()
    n_samples = int(metadata.get("n_samples") or 1)
    need_batch = bool(metadata.get("need_batch_correction")) or n_samples > 1
    integrator = metadata.get("integrator")
    if need_batch or integrator:
        add("harmony-batch-correction", "scvi-tools-single-cell", "bio-single-cell-batch-integration")
    add("bio-single-cell-markers-annotation", "bio-single-cell-cell-annotation", "bio-single-cell-preprocessing")
    if tissue in {"pbmc", "blood", "immune"} or metadata.get("use_popv"):
        add("popv-cell-annotation")
    if metadata.get("use_census") or "atlas" in str(metadata.get("task") or "").lower():
        add("cellxgene-census", "bio-machine-learning-atlas-mapping")
    if language != "python":
        add("scanpy-scrna-seq", "Single-Cell RNA-seq Core Analysis (Seurat)")

    topic = _topic_query(metadata, metadata.get("user_query") or metadata.get("task"))
    ranked = sorted((_skill_score(s, topic), s.name) for s in skills)
    for score, name in reversed(ranked):
        if score >= 6:
            add(name)
        if len(selected) >= 24:
            break
    return selected


def skill_catalog_text(
    root: Path | None = None,
    metadata: dict | None = None,
    query: str | None = None,
    *,
    limit: int | None = _DEFAULT_CATALOG_LIMIT,
) -> str:
    """Catalog for planner prompts: recommended first, then fill up to limit."""
    skills = list_skills(root) if root else list_skills()
    if not skills:
        return "(no skills found)"
    md = dict(metadata or {})
    if query:
        md.setdefault("user_query", query)
        md.setdefault("task", query)
    rec_list = recommend_skills(md) if (metadata or query) else []
    recommended = set(rec_list)
    by_name = {s.name: s for s in skills}
    topic = _topic_query(md, query)
    lines = [f"bundled skills: {len(skills)} (showing up to {limit or len(skills)})"]
    if rec_list:
        lines.append("recommended for this task: " + ", ".join(rec_list))

    ordered: list[Skill] = []
    seen: set[str] = set()
    for name in rec_list:
        if name in by_name and name not in seen:
            ordered.append(by_name[name])
            seen.add(name)
    for skill in sorted(skills, key=lambda s: -_skill_score(s, topic)):
        if skill.name not in seen:
            ordered.append(skill)
            seen.add(skill.name)

    shown = 0
    last_cat = None
    for skill in ordered:
        if limit is not None and shown >= limit:
            remaining = len(skills) - shown
            lines.append(f"... +{remaining} more (`python -m scagent skills`)")
            break
        cat = skill_category(skill)
        if cat != last_cat:
            lines.append(f"## {cat}")
            last_cat = cat
        mark = "*" if skill.name in recommended else "-"
        lines.append(f"{mark} {skill.name}: {_truncate(skill.description or skill.name, 100)}")
        shown += 1
    return "\n".join(lines)


def load_skill_text(name: str, include_references: bool = False, root: Path | None = None) -> str:
    skill = get_skill(name, root)
    if skill is None:
        available = ", ".join(s.name for s in list_skills(root)) or "none"
        return f"Skill '{name}' not found. Available: {available}"
    parts = [f"# {skill.name}\n\n{skill.description}\n\n{skill.body}"]
    if include_references:
        for ref in skill.references:
            parts.append(f"\n\n## Reference: {ref.name}\n\n{ref.read_text(encoding='utf-8')}")
    return "\n".join(parts)


def awesome_manifest() -> dict | None:
    if not _MANIFEST_PATH.exists():
        return None
    try:
        return json.loads(_MANIFEST_PATH.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


PHASE_SKILLS = {
    "qc": ["anndata-data-structure", "scanpy-scrna-seq", "bio-single-cell-preprocessing", "bio-single-cell-doublet-detection"],
    "downstream": [
        "scanpy-scrna-seq",
        "harmony-batch-correction",
        "scvi-tools-single-cell",
        "single-cell-annotation-guide",
        "single-cell-annotation",
        "celltypist-cell-annotation",
        "bio-single-cell-clustering",
        "bio-single-cell-markers-annotation",
    ],
}

# Word-boundary style hints (avoid bare "io" matching "annotation" / "communication").
_PHASE_HINTS = {
    "qc": ("qc", "preprocess", "anndata", "scanpy", "doublet", "data-io", "sparse", "normaliz", "filter"),
    "downstream": (
        "cluster",
        "annotat",
        "harmony",
        "scvi",
        "marker",
        "celltypist",
        "popv",
        "trajectory",
        "cellchat",
        "communication",
        "grn",
        "scenic",
        "velocity",
        "lineage",
        "azimuth",
        "atlas",
        "census",
    ),
}


def _hint_match(name: str, hints: tuple[str, ...]) -> bool:
    low = name.lower()
    for h in hints:
        if len(h) <= 3:
            if re.search(rf"(^|[-_]){re.escape(h)}($|[-_])", low):
                return True
        elif h in low:
            return True
    return False


def skills_for_phase(phase: str, plan_skills: list[str] | None = None, *, max_extra: int = 6) -> list[str]:
    """Core phase skills plus task-selected skills from the catalog."""
    available = {s.name for s in list_skills()}
    wanted = [n for n in (PHASE_SKILLS.get(phase) or []) if n in available]
    hints = _PHASE_HINTS.get(phase) or ()
    extra: list[str] = []
    for name in plan_skills or []:
        if name in wanted or name in extra:
            continue
        if phase in {"qc", "downstream"} and not _hint_match(name, hints):
            continue
        extra.append(name)
        if len(extra) >= max_extra:
            break
    return wanted + extra


SKILLS_ROOT = REPO_ROOT / "skills"
