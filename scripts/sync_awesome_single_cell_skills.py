#!/usr/bin/env python3
"""Sync single-cell skills from awesome-bio-agent-skills into scAgent/skills/."""

from __future__ import annotations

import argparse
import csv
import json
import re
import shutil
from pathlib import Path

from scagent.config import REPO_ROOT

_FRONTMATTER = re.compile(r"^---\s*\n(.*?)\n---\s*\n", re.DOTALL)
_SKIP_DIRS = {"repo", "tests", "demo", "examples", ".git", "__pycache__", "node_modules"}
_SKIP_SUFFIXES = {".pyc", ".pyo"}
_INCLUDE_ROOT_FILES = {
    "SKILL.md",
    "README.md",
    "usage-guide.md",
    "commands_and_thresholds.md",
    "technical_reference.md",
    "agent_config.yaml",
}


def _parse_name(text: str, fallback: str) -> str:
    m = _FRONTMATTER.match(text)
    if not m:
        return fallback
    for line in m.group(1).splitlines():
        if line.startswith("name:"):
            return line.split(":", 1)[1].strip().strip('"').strip("'")
    return fallback


def _copy_skill_tree(src: Path, dst: Path) -> list[str]:
    copied: list[str] = []
    dst.mkdir(parents=True, exist_ok=True)
    for item in sorted(src.iterdir()):
        if item.name in _SKIP_DIRS:
            continue
        if item.is_file():
            if item.suffix in _SKIP_SUFFIXES:
                continue
            if item.name.endswith(".py") and item.name not in {"__init__.py"}:
                continue
            if item.name not in _INCLUDE_ROOT_FILES and item.suffix not in {".md", ".yaml", ".yml"}:
                continue
            shutil.copy2(item, dst / item.name)
            copied.append(str(item.relative_to(src)))
            continue
        if item.name in {"references", "assets"}:
            shutil.copytree(item, dst / item.name, dirs_exist_ok=True)
            copied.append(f"{item.name}/")
    return copied


def sync_skills(
    awesome_root: Path,
    dest_root: Path | None = None,
    *,
    dry_run: bool = False,
) -> dict:
    dest_root = dest_root or (REPO_ROOT / "skills")
    index_path = awesome_root / "bioskill_index_v3.csv"
    if not index_path.exists():
        raise FileNotFoundError(f"Missing index: {index_path}")

    existing_names = set()
    for skill_md in dest_root.glob("*/SKILL.md"):
        text = skill_md.read_text(encoding="utf-8", errors="ignore")
        existing_names.add(_parse_name(text, skill_md.parent.name).lower())

    manifest: list[dict] = []
    imported = skipped = 0

    with index_path.open(newline="", encoding="utf-8") as f:
        rows = [r for r in csv.DictReader(f) if r.get("category") == "single-cell"]

    for row in rows:
        src = awesome_root / "skills" / row["archive_path"]
        skill_md = src / "SKILL.md"
        if not skill_md.exists():
            raise FileNotFoundError(f"Missing SKILL.md: {skill_md}")

        text = skill_md.read_text(encoding="utf-8", errors="ignore")
        skill_name = _parse_name(text, row["folder_name"])
        dest_name = skill_name
        if dest_name.lower() in existing_names:
            skipped += 1
            manifest.append(
                {
                    "name": skill_name,
                    "status": "skipped_existing",
                    "source": row["source_repo"],
                    "archive_path": row["archive_path"],
                    "description": row["description"],
                }
            )
            continue
        if (dest_root / dest_name).exists():
            dest_name = f"{row['source_repo']}--{row['folder_name']}"
        dest = dest_root / dest_name
        entry = {
            "name": skill_name,
            "folder": dest_name,
            "status": "imported",
            "source": row["source_repo"],
            "archive_path": row["archive_path"],
            "description": row["description"],
            "upstream": "https://github.com/BioTender-max/awesome-bio-agent-skills",
        }
        if not dry_run:
            copied = _copy_skill_tree(src, dest)
            entry["files"] = copied
            existing_names.add(skill_name.lower())
            imported += 1
        manifest.append(entry)

    summary = {
        "category": "single-cell",
        "upstream_repo": "https://github.com/BioTender-max/awesome-bio-agent-skills",
        "index_file": "bioskill_index_v3.csv",
        "total_indexed": len(rows),
        "imported": imported,
        "skipped_existing": skipped,
        "skills": manifest,
    }
    if not dry_run:
        manifest_path = dest_root / "awesome_single_cell_manifest.json"
        manifest_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--awesome-root",
        type=Path,
        default=Path("/tmp/awesome-bio-agent-skills"),
        help="Path to cloned awesome-bio-agent-skills repository",
    )
    parser.add_argument("--dest", type=Path, default=REPO_ROOT / "skills")
    parser.add_argument("--dry-run", action="store_true")
    args = parser.parse_args()
    summary = sync_skills(args.awesome_root, args.dest, dry_run=args.dry_run)
    print(json.dumps({k: summary[k] for k in ("total_indexed", "imported", "skipped_existing")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
