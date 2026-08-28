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

# Always activated for a standard scRNA run.
_CORE_SKILLS = (
    "anndata-data-structure",
    "scanpy-scrna-seq",
    "harmony-batch-correction",
    "scvi-tools-single-cell",
    "single-cell-annotation-guide",
    "single-cell-annotation",
    "celltypist-cell-annotation",
    "popv-cell-annotation",
    "cellxgene-census",
    "cellchat-cell-communication",
)

_CATEGORIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("spatial", ("spatial", "visium", "giotto", "squidpy", "xenium", "merfish", "deconv", "空转")),
    ("trajectory", ("trajectory", "pseudotime", "velocity", "lineage", "monocle", "palantir", "拟时序", "轨迹")),
    ("communication", ("communication", "cellchat", "nichenet", "liana", "ligand", "cellphonedb", "通讯", "配体")),
    ("annotation", ("annotat", "celltypist", "azimuth", "singler", "popv", "cell-type", "注释")),
    ("integration", ("batch", "harmony", "scvi", "integrat", "整合", "批次")),
    ("atac", ("atac", "multiome", "signac", "archr", "chromatin", "cicero")),
    ("grn", ("grn", "scenic", "regulon", "arboreto", "gene-regulatory", "调控")),
    ("perturb", ("perturb", "crispr", "crop-seq", "screen")),
    ("qc", ("qc", "preprocess", "doublet", "normaliz", "质控", "预处理")),
    ("clustering", ("cluster", "leiden", "pca", "umap", "聚类")),
    ("multimodal", ("multimodal", "cite-seq", "muon", "protein", "adt")),
    ("repertoire", ("tcr", "bcr", "repertoire", "scirpy")),
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


def _parse_frontmatter(text: str) -> tuple[dict[str, str], str]:
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
    """Score every bundled skill against the task; always keep SciAgent core."""
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

    topic = _topic_query(metadata, metadata.get("user_query"))
    ranked = sorted((_skill_score(s, topic), s.name) for s in skills)
    for score, name in reversed(ranked):
        if score >= 6:
            add(name)
        if len(selected) >= 24:
            break
    return selected


def _ranked_skills(metadata: dict | None = None, query: str | None = None) -> list[Skill]:
    skills = list_skills()
    if not skills:
        return []
    topic = _topic_query(metadata, query)
    recommended = set(recommend_skills(metadata or {}, language=str((metadata or {}).get("language") or "python")))
    ranked: list[Skill] = []
    seen: set[str] = set()

    def push(skill: Skill) -> None:
        if skill.name in seen:
            return
        seen.add(skill.name)
        ranked.append(skill)

    for skill in skills:
        if skill.name in recommended:
            push(skill)
    for skill in sorted(skills, key=lambda s: -_skill_score(s, topic)):
        push(skill)
    return ranked


def skill_catalog_text(
    root: Path | None = None,
    metadata: dict | None = None,
    query: str | None = None,
    *,
    limit: int | None = None,
) -> str:
    """Full catalog grouped by topic so the planner can see every bundled skill."""
    skills = list_skills(root) if root else list_skills()
    if not skills:
        return "(no skills found)"
    md = dict(metadata or {})
    if query:
        md.setdefault("user_query", query)
        md.setdefault("task", query)
    rec_list = recommend_skills(md) if (metadata or query) else []
    recommended = set(rec_list)
    grouped: dict[str, list[Skill]] = defaultdict(list)
    for skill in skills:
        grouped[skill_category(skill)].append(skill)
    order = [c for c, _ in _CATEGORIES] + ["other"]
    lines = [f"bundled skills: {len(skills)}"]
    if rec_list:
        lines.append("recommended for this task: " + ", ".join(rec_list))
    shown = 0
    for cat in order:
        bucket = grouped.get(cat) or []
        if not bucket:
            continue
        lines.append(f"## {cat}")
        for skill in bucket:
            if limit is not None and shown >= limit:
                remaining = len(skills) - shown
                lines.append(f"... 另有 {remaining} 个 skill（`python -m scagent skills`）")
                return "\n".join(lines)
            mark = "*" if skill.name in recommended else "-"
            lines.append(f"{mark} {skill.name}: {_truncate(skill.description, 100)}")
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

_PHASE_HINTS = {
    "qc": ("qc", "preprocess", "anndata", "scanpy", "doublet", "data-io", "sparse", "normaliz", "io"),
    "downstream": (
        "cluster",
        "annotat",
        "harmony",
        "scvi",
        "marker",
        "celltypist",
        "popv",
        "trajectory",
        "spatial",
        "cellchat",
        "communication",
        "perturb",
        "atac",
        "grn",
        "scenic",
        "velocity",
        "lineage",
        "azimuth",
        "multimodal",
        "cite",
        "tcr",
        "repertoire",
        "deconv",
    ),
}


def skills_for_phase(phase: str, plan_skills: list[str] | None = None, *, max_extra: int = 6) -> list[str]:
    """Core phase skills plus task-selected skills from the full catalog."""
    available = {s.name for s in list_skills()}
    wanted = [n for n in (PHASE_SKILLS.get(phase) or []) if n in available]
    hints = _PHASE_HINTS.get(phase) or ()
    extra: list[str] = []
    for name in plan_skills or []:
        if name in wanted or name in extra:
            continue
        low = name.lower()
        if phase == "qc" and not any(h in low for h in hints):
            continue
        if phase == "downstream" and not any(h in low for h in hints):
            continue
        extra.append(name)
        if len(extra) >= max_extra:
            break
    return wanted + extra


# Keep a stable import for tests that want the repo root.
SKILLS_ROOT = REPO_ROOT / "skills"
