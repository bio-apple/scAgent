"""Discover and recommend scientific-task skills (one skill ≈ one analysis stage)."""

from __future__ import annotations

import json
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from scagent.config import REPO_ROOT, load_config, resolve_path

_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---\s*\n(.*)$", re.DOTALL)
_MANIFEST_PATH = REPO_ROOT / "skills" / "awesome_single_cell_manifest.json"

# Canonical scientific-task skills (order = default analysis DAG).
TASK_SKILLS: tuple[str, ...] = (
    "dataset_loader",
    "qc_preprocessing",
    "integration_batch",
    "clustering_embedding",
    "cell_annotation",
    "deg_pathway",
    "trajectory",
    "cell_communication",
    "visualization",
    "report_generation",
)

# Always include for a standard scRNA run.
_CORE_SKILLS = (
    "dataset_loader",
    "qc_preprocessing",
    "clustering_embedding",
    "cell_annotation",
    "visualization",
    "report_generation",
)

_DEFAULT_CATALOG_LIMIT = 20

_CATEGORIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("io", ("dataset_loader", "load", "10x", "h5ad")),
    ("qc", ("qc_preprocessing", "preprocess", "doublet", "质控")),
    ("integration", ("integration_batch", "harmony", "batch", "整合")),
    ("clustering", ("clustering_embedding", "leiden", "umap", "pca", "聚类")),
    ("annotation", ("cell_annotation", "celltypist", "marker", "注释")),
    ("deg", ("deg_pathway", "differential", "gsea", "pathway", "差异")),
    ("trajectory", ("trajectory", "pseudotime", "velocity", "轨迹")),
    ("communication", ("cell_communication", "cellchat", "ligand", "通讯")),
    ("viz", ("visualization", "plot", "figure")),
    ("report", ("report_generation", "report", "methods")),
)


@dataclass(frozen=True)
class Skill:
    name: str
    description: str
    path: Path
    body: str
    references: tuple[Path, ...]


def _strip_leading_noise(text: str) -> str:
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
    # Stable DAG order for known task skills.
    rank = {n: i for i, n in enumerate(TASK_SKILLS)}
    skills.sort(key=lambda s: (rank.get(s.name, 1000), s.name))
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


def recommend_skills(metadata: dict, language: str = "python") -> list[str]:
    """Select scientific-task skills from route / design (not algorithm primitives)."""
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
    route = [str(x).lower() for x in (metadata.get("route") or [])]
    intents = [str(x).lower() for x in (metadata.get("intents") or [])]
    topic = " ".join(
        [
            str(metadata.get("user_query") or ""),
            str(metadata.get("task") or ""),
            " ".join(route),
            " ".join(intents),
        ]
    ).lower()

    n_samples = int(metadata.get("n_samples") or 1)
    need_batch = bool(metadata.get("need_batch_correction")) or n_samples > 1
    # Sample≡condition collinearity: skip integration skill (over-correction risk).
    confounded = bool(metadata.get("batch_condition_confounded"))
    want_integration = (
        (need_batch or metadata.get("integrator") or "harmony" in topic or "integrat" in topic)
        and not confounded
    )
    if want_integration:
        add("integration_batch")

    if any(k in topic for k in ("deg", "differen", "marker", "gsea", "pathway", "pseudobulk", "差异", "通路")) or any(
        x in route for x in ("deg", "gsea", "enrichment")
    ):
        add("deg_pathway")

    if any(k in topic for k in ("traject", "pseudotime", "velocity", "monocle", "slingshot", "命运", "轨迹")) or any(
        x in route for x in ("trajectory", "paga", "dpt")
    ):
        add("trajectory")

    if any(k in topic for k in ("cellchat", "communication", "ligand", "nichenet", "通讯", "配体")):
        add("cell_communication")

    # Always finish with viz + report when core ran.
    add("visualization", "report_generation")
    return selected


def skill_catalog_text(
    root: Path | None = None,
    metadata: dict | None = None,
    query: str | None = None,
    *,
    limit: int | None = _DEFAULT_CATALOG_LIMIT,
) -> str:
    skills = list_skills(root) if root else list_skills()
    if not skills:
        return "(no skills found)"
    md = dict(metadata or {})
    if query:
        md.setdefault("user_query", query)
        md.setdefault("task", query)
    rec_list = recommend_skills(md) if (metadata or query) else list(TASK_SKILLS)
    recommended = set(rec_list)
    lines = [
        f"bundled scientific-task skills: {len(skills)}",
        "principle: one skill = one research task (not PCA/UMAP/Leiden primitives)",
    ]
    if rec_list:
        lines.append("recommended for this task: " + ", ".join(rec_list))
    lines.append("## tasks")
    shown = 0
    for skill in skills:
        if limit is not None and shown >= limit:
            lines.append(f"... +{len(skills) - shown} more (`python -m scagent skills`)")
            break
        mark = "*" if skill.name in recommended else "-"
        lines.append(f"{mark} {skill.name}: {_truncate(skill.description or skill.name, 110)}")
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
    "qc": ["dataset_loader", "qc_preprocessing"],
    "downstream": [
        "integration_batch",
        "clustering_embedding",
        "cell_annotation",
        "deg_pathway",
        "visualization",
    ],
}

_PHASE_HINTS = {
    "qc": ("dataset_loader", "qc_preprocessing", "qc", "preprocess", "load"),
    "downstream": (
        "integration_batch",
        "clustering_embedding",
        "cell_annotation",
        "deg_pathway",
        "trajectory",
        "cell_communication",
        "visualization",
        "report_generation",
        "cluster",
        "annotat",
        "integrat",
        "deg",
        "pathway",
        "traject",
        "communicat",
        "visual",
        "report",
    ),
}


def _hint_match(name: str, hints: tuple[str, ...]) -> bool:
    low = name.lower()
    return any(h in low for h in hints)


def skills_for_phase(phase: str, plan_skills: list[str] | None = None, *, max_extra: int = 6) -> list[str]:
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
