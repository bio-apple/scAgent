from __future__ import annotations

import json
import re
from dataclasses import dataclass
from pathlib import Path

from scagent.config import REPO_ROOT, load_config, resolve_path

_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---\s*\n(.*)$", re.DOTALL)
_MANIFEST_PATH = REPO_ROOT / "skills" / "awesome_single_cell_manifest.json"

# Core SciAgent skills always considered for standard scRNA workflows.
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

_TOPIC_KEYWORDS: tuple[tuple[tuple[str, ...], tuple[str, ...]], ...] = (
    (("spatial", "visium", "squidpy", "giotto", "空间", "空转"), ("spatial", "visium", "squidpy", "giotto")),
    (("trajectory", "pseudotime", "monocle", "palantir", "velocity", "轨迹", "拟时序"), ("trajectory", "pseudotime", "velocity", "lineage")),
    (("doublet", "scrublet", "doubletfinder", "双细胞"), ("doublet",)),
    (("scatac", "atac", "multiome", "signac", "archr"), ("atac", "multiome", "signac")),
    (("perturb", "crispr screen", "crop-seq"), ("perturb", "crispr")),
    (("scenic", "regulon", "grn", "gene regulatory"), ("scenic", "regulon", "grn")),
    (("deconv", "去卷积"), ("deconv",)),
    (("cell communication", "cellchat", "nichenet", "liana", "ligand", "细胞通讯", "配体"), ("communication", "cellchat", "nichenet", "liana")),
    (("marker", "annotation", "cell type", "注释"), ("marker", "annotation", "cell-annotation")),
    (("batch", "integration", "harmony", "批次", "整合"), ("batch", "integration")),
    (("preprocess", "qc", "normalize", "质控", "预处理"), ("preprocess", "qc", "scrna", "scanpy")),
    (("differential", "deg", "pseudobulk", "差异"), ("differential", "deg", "pseudobulk", "expression")),
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


def _truncate(text: str, limit: int = 140) -> str:
    text = " ".join(text.split())
    if len(text) <= limit:
        return text
    return text[: limit - 1].rstrip() + "…"


def _topic_query(metadata: dict | None, query: str | None) -> str:
    parts = [
        str((metadata or {}).get("task") or ""),
        str((metadata or {}).get("query") or ""),
        str(query or ""),
        str((metadata or {}).get("tissue") or ""),
        str((metadata or {}).get("platform") or ""),
    ]
    return " ".join(parts).lower()


def _skill_matches_topic(skill: Skill, topic: str) -> bool:
    hay = f"{skill.name} {skill.description}".lower()
    for triggers, needles in _TOPIC_KEYWORDS:
        if any(t in topic for t in triggers):
            if any(n in hay for n in needles):
                return True
    return False


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
    if topic:
        for skill in skills:
            if _skill_matches_topic(skill, topic):
                push(skill)
    for skill in skills:
        if skill.name in _CORE_SKILLS:
            push(skill)
    for skill in skills:
        push(skill)
    return ranked


def skill_catalog_text(
    root: Path | None = None,
    metadata: dict | None = None,
    query: str | None = None,
    *,
    limit: int = 48,
) -> str:
    skills = _ranked_skills(metadata, query)
    if not skills:
        return "(no skills found)"
    lines: list[str] = []
    for skill in skills[:limit]:
        lines.append(f"- {skill.name}: {_truncate(skill.description)}")
    if len(skills) > limit:
        lines.append(f"... 另有 {len(skills) - limit} 个 skill（`python -m scagent skills` 查看完整列表）")
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


def recommend_skills(metadata: dict, language: str = "python") -> list[str]:
    """Map analysis context onto bundled skills (SciAgent core + awesome-bio single-cell)."""
    names = [s.name for s in list_skills()]
    if not names:
        return []
    selected: list[str] = []

    def add(*cands: str) -> None:
        for c in cands:
            if c in names and c not in selected:
                selected.append(c)

    add("anndata-data-structure", "scanpy-scrna-seq")
    tissue = str(metadata.get("tissue") or "").lower()
    n_samples = int(metadata.get("n_samples") or 1)
    need_batch = bool(metadata.get("need_batch_correction")) or n_samples > 1
    integrator = metadata.get("integrator")
    if need_batch or integrator:
        if integrator == "scvi":
            add("scvi-tools-single-cell", "harmony-batch-correction", "bio-single-cell-batch-integration")
        else:
            add("harmony-batch-correction", "scvi-tools-single-cell", "bio-single-cell-batch-integration")
    add("single-cell-annotation-guide", "single-cell-annotation", "celltypist-cell-annotation")
    add("bio-single-cell-markers-annotation", "bio-single-cell-cell-annotation")
    if tissue in {"pbmc", "blood", "immune"} or metadata.get("use_popv"):
        add("popv-cell-annotation")
    if metadata.get("use_census") or "atlas" in str(metadata.get("task") or "").lower():
        add("cellxgene-census", "bio-machine-learning-atlas-mapping")
    task = _topic_query(metadata, metadata.get("user_query"))
    if any(
        k in task
        for k in (
            "cellchat",
            "cell chat",
            "cell-cell",
            "cell communication",
            "ligand-receptor",
            "ligand receptor",
            "细胞通讯",
            "配体受体",
            "配体-受体",
        )
    ):
        add("cellchat-cell-communication", "bio-single-cell-cell-communication", "cell-communication")
    if any(k in task for k in ("doublet", "scrublet", "doubletfinder", "双细胞")):
        add("bio-single-cell-doublet-detection", "scrna-preprocessing-clustering")
    if any(k in task for k in ("spatial", "visium", "squidpy", "giotto", "空间", "空转")):
        add("spatial-transcriptomics", "spatial-data-io", "spatial-preprocessing")
    if any(k in task for k in ("trajectory", "pseudotime", "monocle", "palantir", "velocity", "轨迹", "拟时序")):
        add(
            "trajectory-lineage",
            "bio-single-cell-trajectory-inference",
            "Single-Cell Trajectory Inference",
            "rna-velocity-agent",
        )
    if any(k in task for k in ("scatac", "atac", "multiome", "signac")):
        add("multiome-scatac", "bio-atac-seq-single-cell-atac")
    if any(k in task for k in ("perturb", "crop-seq", "crispr screen")):
        add("bio-single-cell-perturb-seq", "bio-crispr-screens-perturb-seq-analysis")
    if any(k in task for k in ("scenic", "regulon", "grn")):
        add("bio-gene-regulatory-networks-scenic-regulons")
    if language != "python":
        add("scanpy-scrna-seq")
    return selected


# Keep a stable import for tests that want the repo root.
SKILLS_ROOT = REPO_ROOT / "skills"
