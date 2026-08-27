from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from scagent.config import REPO_ROOT, load_config, resolve_path

_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---\s*\n(.*)$", re.DOTALL)


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


def skill_catalog_text(root: Path | None = None) -> str:
    lines = []
    for s in list_skills(root):
        lines.append(f"- {s.name}: {s.description}")
    return "\n".join(lines) if lines else "(no skills found)"


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


def recommend_skills(metadata: dict, language: str = "python") -> list[str]:
    """Map analysis context onto the existing SciAgent-style skills (do not invent new ones)."""
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
    if need_batch:
        add("harmony-batch-correction", "scvi-tools-single-cell")
    add("single-cell-annotation-guide", "celltypist-cell-annotation")
    if tissue in {"pbmc", "blood", "immune"} or metadata.get("use_popv"):
        add("popv-cell-annotation")
    if metadata.get("use_census") or "atlas" in str(metadata.get("task") or "").lower():
        add("cellxgene-census")
    if language != "python":
        # Existing skills are Python-first; keep them listed as the executable SOP.
        add("scanpy-scrna-seq")
    return selected


# Keep a stable import for tests that want the repo root.
SKILLS_ROOT = REPO_ROOT / "skills"
