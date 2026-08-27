from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

from scagent.config import REPO_ROOT, load_config, resolve_path

_TOKEN_RE = re.compile(r"[a-z0-9_]+|[\u4e00-\u9fff]", re.I)


def tokenize(text: str) -> list[str]:
    return [t.lower() for t in _TOKEN_RE.findall(text or "")]


def chunk_text(text: str, chunk_size: int, overlap: int) -> list[str]:
    text = re.sub(r"\n{3,}", "\n\n", text).strip()
    if not text:
        return []
    if len(text) <= chunk_size:
        return [text]
    chunks: list[str] = []
    start = 0
    while start < len(text):
        end = min(len(text), start + chunk_size)
        chunks.append(text[start:end].strip())
        if end >= len(text):
            break
        start = max(0, end - overlap)
    return [c for c in chunks if c]


def _read_pdf(path: Path) -> str:
    from pypdf import PdfReader

    reader = PdfReader(str(path))
    pages = []
    for page in reader.pages:
        pages.append(page.extract_text() or "")
    return "\n".join(pages)


def _iter_source_files(collection_dir: Path) -> Iterable[Path]:
    if not collection_dir.exists():
        return []
    files: list[Path] = []
    for ext in ("*.md", "*.txt", "*.pdf"):
        files.extend(collection_dir.rglob(ext))
    return [p for p in files if p.is_file() and p.name != ".DS_Store"]


def collection_dir(cfg: dict, collection: str) -> Path:
    """Directory for a RAG collection. Override with rag.collection_dirs.<name>."""
    mapping = (cfg.get("rag") or {}).get("collection_dirs") or {}
    if collection in mapping:
        p = Path(mapping[collection])
        if not p.is_absolute():
            p = Path(cfg.get("_root") or REPO_ROOT) / p
        return p
    return resolve_path(cfg, "knowledge") / collection


def ingest(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    knowledge = resolve_path(cfg, "knowledge")
    index_dir = resolve_path(cfg, "index")
    index_dir.mkdir(parents=True, exist_ok=True)
    chunk_size = int(cfg["rag"]["chunk_size"])
    overlap = int(cfg["rag"]["chunk_overlap"])

    records: list[dict] = []
    root = Path(cfg.get("_root") or REPO_ROOT)
    for collection in cfg["rag"]["collections"]:
        col_dir = collection_dir(cfg, collection)
        for path in _iter_source_files(col_dir):
            if path.name == "README.md":
                continue
            if path.suffix.lower() == ".pdf":
                text = _read_pdf(path)
            else:
                text = path.read_text(encoding="utf-8", errors="replace")
            try:
                rel = str(path.relative_to(root))
            except ValueError:
                rel = str(path.relative_to(knowledge))
            for i, chunk in enumerate(chunk_text(text, chunk_size, overlap)):
                records.append(
                    {
                        "id": f"{rel}::{i}",
                        "collection": collection,
                        "source": rel,
                        "chunk_index": i,
                        "text": chunk,
                    }
                )

    out = index_dir / "chunks.jsonl"
    with out.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    (index_dir / "meta.json").write_text(
        json.dumps({"n_chunks": len(records)}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return out
