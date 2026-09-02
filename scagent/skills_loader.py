"""Discover and recommend capability-based skills (6 analysis capabilities)."""

from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

from scagent.config import REPO_ROOT, load_config, resolve_path

_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---\s*\n(.*)$", re.DOTALL)
_MANIFEST_PATH = REPO_ROOT / "skills" / "awesome_single_cell_manifest.json"

# Canonical capability skills (diagram order).
CAPABILITY_SKILLS: tuple[str, ...] = (
    "seurat-workflow",
    "cell-annotation",
    "differential-expression",
    "cell-communication",
    "trajectory",
    "spatial-analysis",
)

# Backward compatibility: legacy task skill names → capability.
LEGACY_SKILL_ALIASES: dict[str, str] = {
    "dataset_loader": "seurat-workflow",
    "qc_preprocessing": "seurat-workflow",
    "integration_batch": "seurat-workflow",
    "clustering_embedding": "seurat-workflow",
    "visualization": "seurat-workflow",
    "report_generation": "seurat-workflow",
    "cell_annotation": "cell-annotation",
    "deg_pathway": "differential-expression",
    "cell_communication": "cell-communication",
}

# Deprecated names kept for tests/docs.
TASK_SKILLS: tuple[str, ...] = tuple(LEGACY_SKILL_ALIASES.keys())

_CORE_CAPABILITIES = ("seurat-workflow", "cell-annotation")

_DEFAULT_CATALOG_LIMIT = 20

_CATEGORIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("workflow", ("seurat-workflow", "qc", "cluster", "integrat", "load", "umap", "leiden")),
    ("annotation", ("cell-annotation", "celltypist", "marker", "注释", "azimuth")),
    ("deg", ("differential-expression", "pseudobulk", "deg", "差异", "pathway")),
    ("trajectory", ("trajectory", "pseudotime", "velocity", "monocle", "轨迹")),
    ("communication", ("cell-communication", "cellchat", "ligand", "通讯")),
    ("spatial", ("spatial-analysis", "spatial", "visium", "squidpy", "空间")),
)


@dataclass(frozen=True)
class Skill:
    name: str
    description: str
    path: Path
    body: str
    references: tuple[Path, ...]
    deprecated: bool = False
    replaced_by: str | None = None


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


def _truthy(val: str | None) -> bool:
    return str(val or "").lower() in {"true", "1", "yes"}


def skills_dir(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    return resolve_path(cfg, "skills")


def _resolve_alias(name: str) -> str:
    key = name.strip().lower().replace("_", "-")
    for legacy, cap in LEGACY_SKILL_ALIASES.items():
        if key == legacy.lower().replace("_", "-") or key == legacy.lower():
            return cap
    return name.strip()


def list_skills(root: Path | None = None, *, include_legacy: bool = False) -> list[Skill]:
    base = root or skills_dir()
    skills: list[Skill] = []
    if not base.exists():
        return skills
    for skill_md in sorted(base.glob("*/SKILL.md")):
        if skill_md.parent.name.startswith("_"):
            continue
        text = skill_md.read_text(encoding="utf-8")
        meta, body = _parse_frontmatter(text)
        deprecated = _truthy(meta.get("deprecated"))
        if deprecated and not include_legacy:
            continue
        name = meta.get("name") or skill_md.parent.name.replace("_", "-")
        refs = tuple(sorted(skill_md.parent.joinpath("references").glob("*.md")))
        skills.append(
            Skill(
                name=name,
                description=meta.get("description", ""),
                path=skill_md,
                body=body,
                references=refs,
                deprecated=deprecated,
                replaced_by=meta.get("replaced_by"),
            )
        )
    rank = {n: i for i, n in enumerate(CAPABILITY_SKILLS)}
    skills.sort(key=lambda s: (rank.get(s.name, 1000), s.name))
    return skills


def get_skill(name: str, root: Path | None = None) -> Skill | None:
    resolved = _resolve_alias(name)
    key = resolved.lower()
    for skill in list_skills(root, include_legacy=True):
        if skill.name.lower() == key or skill.path.parent.name.lower().replace("_", "-") == key.replace("_", "-"):
            if skill.deprecated and skill.replaced_by:
                return get_skill(skill.replaced_by, root)
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
    """Select capability skills from route / design."""
    skills = list_skills()
    if not skills:
        return []
    names = {s.name for s in skills}
    selected: list[str] = []

    def add(*cands: str) -> None:
        for c in cands:
            c = _resolve_alias(c)
            if c in names and c not in selected:
                selected.append(c)

    add(*_CORE_CAPABILITIES)
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
    confounded = bool(metadata.get("batch_condition_confounded"))
    if (need_batch or metadata.get("integrator") or "harmony" in topic or "integrat" in topic) and not confounded:
        add("seurat-workflow")  # integration is part of core workflow

    if any(k in topic for k in ("deg", "differen", "pseudobulk", "差异", "pathway", "gsea")) or any(
        x in route for x in ("deg", "gsea", "enrichment", "pseudobulk")
    ):
        add("differential-expression")

    if any(k in topic for k in ("traject", "pseudotime", "velocity", "monocle", "slingshot", "命运", "轨迹")) or any(
        x in route for x in ("trajectory", "paga", "dpt")
    ):
        add("trajectory")

    if any(k in topic for k in ("cellchat", "communication", "ligand", "nichenet", "通讯", "配体")):
        add("cell-communication")

    if any(k in topic for k in ("spatial", "visium", "xenium", "squidpy", "空间", "st ", " slide-seq")):
        add("spatial-analysis")

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
    rec_list = recommend_skills(md) if (metadata or query) else list(CAPABILITY_SKILLS)
    recommended = set(rec_list)
    lines = [
        f"capability skills: {len(skills)}",
        "principle: one skill = one analysis capability (not PCA/UMAP/Leiden primitives)",
    ]
    if rec_list:
        lines.append("recommended for this task: " + ", ".join(rec_list))
    lines.append("## capabilities")
    shown = 0
    for skill in skills:
        if limit is not None and shown >= limit:
            lines.append(f"... +{len(skills) - shown} more (`python -m scagent skills`)")
            break
        mark = "*" if skill.name in recommended else "-"
        lines.append(f"{mark} {skill.name}: {_truncate(skill.description or skill.name, 110)}")
        shown += 1
    legacy = list_skills(include_legacy=True)
    dep = [s for s in legacy if s.deprecated]
    if dep:
        lines.append("## deprecated (alias → capability)")
        for s in dep[:8]:
            lines.append(f"- {s.path.parent.name} → {s.replaced_by or '?'}")
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
    "qc": ["seurat-workflow"],
    "cluster": ["seurat-workflow"],
    "downstream": ["seurat-workflow", "cell-annotation", "differential-expression"],
    "annotate": ["cell-annotation", "differential-expression"],
}

_PHASE_HINTS = {
    "qc": ("seurat-workflow", "seurat", "qc", "preprocess", "load", "workflow"),
    "cluster": ("seurat-workflow", "cluster", "leiden", "umap", "integrat"),
    "downstream": (
        "seurat-workflow",
        "cell-annotation",
        "differential-expression",
        "trajectory",
        "cell-communication",
        "spatial-analysis",
        "cluster",
        "annotat",
        "integrat",
        "deg",
        "traject",
        "communicat",
    ),
    "annotate": ("cell-annotation", "differential-expression", "annotat", "marker", "deg"),
}


def _hint_match(name: str, hints: tuple[str, ...]) -> bool:
    low = name.lower()
    return any(h in low for h in hints)


def skills_for_phase(phase: str, plan_skills: list[str] | None = None, *, max_extra: int = 6) -> list[str]:
    available = {s.name for s in list_skills()}
    wanted = [_resolve_alias(n) for n in (PHASE_SKILLS.get(phase) or [])]
    wanted = [n for n in wanted if n in available]
    hints = _PHASE_HINTS.get(phase) or ()
    extra: list[str] = []
    for name in plan_skills or []:
        resolved = _resolve_alias(name)
        if resolved in wanted or resolved in extra:
            continue
        if phase in {"qc", "cluster", "downstream", "annotate"} and not _hint_match(resolved, hints):
            continue
        if resolved not in available:
            continue
        extra.append(resolved)
        if len(extra) >= max_extra:
            break
    out: list[str] = []
    for n in wanted + extra:
        if n not in out:
            out.append(n)
    return out


SKILLS_ROOT = REPO_ROOT / "skills"
