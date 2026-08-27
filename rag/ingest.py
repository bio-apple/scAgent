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


def _units(text: str) -> list[str]:
    text = re.sub(r"\n{3,}", "\n\n", text or "").strip()
    if not text:
        return []
    marks = [0]
    for m in re.finditer(
        r"(?m)^(?:#{1,6}\s+|第[一二三四五六七八九十0-9]+[章节篇].*|\d+\.\d+\s+\S)",
        text,
    ):
        if m.start() > 0:
            marks.append(m.start())
    marks.append(len(text))
    sections: list[str] = []
    for a, b in zip(marks, marks[1:]):
        bit = text[a:b].strip()
        if bit:
            sections.append(bit)
    units: list[str] = []
    for sec in sections:
        paras = [p.strip() for p in re.split(r"\n\s*\n", sec) if p.strip()]
        units.extend(paras or [sec])
    return units


def _sentences(text: str) -> list[str]:
    parts = re.split(r"(?<=[。．.!?！？])\s+", text)
    return [p.strip() for p in parts if p.strip()]


def chunk_semantic(text: str, chunk_size: int, overlap: int) -> list[str]:
    """Pack heading/paragraph units; do not cut a paragraph unless it exceeds chunk_size."""
    units = _units(text)
    if not units:
        return []
    chunks: list[str] = []
    buf: list[str] = []
    size = 0

    def flush() -> None:
        nonlocal buf, size
        if not buf:
            return
        chunks.append("\n\n".join(buf))
        if overlap > 0:
            kept: list[str] = []
            acc = 0
            for u in reversed(buf):
                if acc + len(u) > overlap:
                    break
                kept.append(u)
                acc += len(u)
            buf = list(reversed(kept))
            size = sum(len(x) for x in buf) + 2 * max(len(buf) - 1, 0)
        else:
            buf, size = [], 0

    for unit in units:
        if len(unit) > chunk_size:
            flush()
            sent_buf: list[str] = []
            sent_size = 0
            for sent in _sentences(unit):
                if sent_size + len(sent) > chunk_size and sent_buf:
                    chunks.append(" ".join(sent_buf))
                    if overlap:
                        sent_buf, sent_size = [sent_buf[-1]], len(sent_buf[-1])
                    else:
                        sent_buf, sent_size = [], 0
                sent_buf.append(sent)
                sent_size += len(sent) + 1
            if sent_buf:
                chunks.append(" ".join(sent_buf))
            buf, size = [], 0
            continue
        if buf and size + len(unit) + 2 > chunk_size:
            flush()
        buf.append(unit)
        size += len(unit) + 2
    if buf:
        chunks.append("\n\n".join(buf))
    return [c for c in chunks if c]


def chunk_document(text: str, chunk_size: int, overlap: int, *, mode: str = "semantic") -> list[str]:
    if str(mode).lower() == "fixed":
        return chunk_text(text, chunk_size, overlap)
    return chunk_semantic(text, chunk_size, overlap)


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


def _index_fresh(cfg: dict, out: Path) -> bool:
    if not out.exists() or out.stat().st_size == 0:
        return False
    newest = 0.0
    for collection in cfg["rag"]["collections"]:
        for path in _iter_source_files(collection_dir(cfg, collection)):
            if path.name == "README.md":
                continue
            newest = max(newest, path.stat().st_mtime)
    return newest <= out.stat().st_mtime


def collection_dir(cfg: dict, collection: str) -> Path:
    """Directory for a RAG collection. Override with rag.collection_dirs.<name>."""
    mapping = (cfg.get("rag") or {}).get("collection_dirs") or {}
    if collection in mapping:
        p = Path(mapping[collection])
        if not p.is_absolute():
            p = Path(cfg.get("_root") or REPO_ROOT) / p
        return p
    return resolve_path(cfg, "knowledge") / collection


def ingest(cfg: dict | None = None, *, force: bool = False) -> Path:
    cfg = cfg or load_config()
    knowledge = resolve_path(cfg, "knowledge")
    index_dir = resolve_path(cfg, "index")
    index_dir.mkdir(parents=True, exist_ok=True)
    out = index_dir / "chunks.jsonl"
    if not force and _index_fresh(cfg, out):
        return out
    chunk_size = int(cfg["rag"]["chunk_size"])
    overlap = int(cfg["rag"]["chunk_overlap"])
    chunk_mode = str((cfg.get("rag") or {}).get("chunking") or "semantic")

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
            for i, chunk in enumerate(chunk_document(text, chunk_size, overlap, mode=chunk_mode)):
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
        json.dumps({"n_chunks": len(records), "chunking": chunk_mode}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return out
