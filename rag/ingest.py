from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Iterable

from scagent.config import REPO_ROOT, load_config, resolve_path

_TOKEN_RE = re.compile(r"[a-z0-9_]+|[\u4e00-\u9fff]", re.I)

SKIP_DIR_NAMES = {
    "_build",
    "figures",
    "datasets",
    ".git",
    "__pycache__",
    ".ipynb_checkpoints",
    "node_modules",
    ".github",
    "changelog.d",
    "scripts",
    "template",
    "dropdowns",
}
SKIP_FILE_NAMES = {
    ".DS_Store",
    "README.md",
    "CONTRIBUTING.md",
    "CHANGELOG.md",
    "CODE_OF_CONDUCT.md",
    "LICENSE.md",
    "LICENSE",
}


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


def _read_pdf(path: Path, *, cfg: dict | None = None) -> str:
    from scagent.knowledge.parser import extract_pdf_text

    cfg = cfg or load_config()
    papers_cfg = (cfg.get("rag") or {}).get("papers") or {}
    backend = str(papers_cfg.get("parser_backend") or "mineru")
    text, _ = extract_pdf_text(path, backend=backend, cfg=cfg)
    return text


def _read_ipynb(path: Path) -> str:
    data = json.loads(path.read_text(encoding="utf-8", errors="replace"))
    parts: list[str] = []
    for cell in data.get("cells") or []:
        src = cell.get("source")
        if isinstance(src, list):
            src = "".join(src)
        text = str(src or "").strip()
        if not text:
            continue
        kind = cell.get("cell_type")
        if kind == "markdown":
            parts.append(text)
        elif kind == "code":
            parts.append("```\n" + text + "\n```")
    return "\n\n".join(parts)


def _iter_papers_files(cfg: dict, col_dir: Path) -> Iterable[tuple[Path, str]]:
    """Yield (path, kind) for papers collection: curated md, parsed md; skip raw pdf when parsed exists."""
    from scagent.knowledge.parser import parsed_paths

    papers_cfg = (cfg.get("rag") or {}).get("papers") or {}
    index_parsed_only = papers_cfg.get("index_parsed_only", True)
    _, parsed_dir = parsed_paths(cfg)
    parsed_stems = {p.stem for p in parsed_dir.glob("*.md")} if parsed_dir.exists() else set()

    for path in sorted(col_dir.glob("*.md")):
        if path.name.startswith("."):
            continue
        yield path, "md"

    if parsed_dir.exists():
        for path in sorted(parsed_dir.glob("*.md")):
            yield path, "parsed"

    if not index_parsed_only:
        for path in sorted(col_dir.glob("*.pdf")):
            if path.stem in parsed_stems:
                continue
            yield path, "pdf"
    else:
        for path in sorted(col_dir.glob("*.pdf")):
            if path.stem not in parsed_stems:
                yield path, "pdf"


def _iter_source_files(collection_dir: Path, *, collection: str = "", cfg: dict | None = None) -> Iterable[Path]:
    if not collection_dir.exists():
        return []
    files: list[Path] = []
    for ext in ("*.md", "*.txt", "*.pdf", "*.ipynb", "*.json"):
        files.extend(collection_dir.rglob(ext))
    out: list[Path] = []
    for p in files:
        if not p.is_file() or p.name in SKIP_FILE_NAMES:
            continue
        if any(part in SKIP_DIR_NAMES for part in p.parts):
            continue
        out.append(p)
    return out


def _iter_collection_files(cfg: dict, collection: str, col_dir: Path) -> Iterable[tuple[Path, str]]:
    if collection == "papers":
        yield from _iter_papers_files(cfg, col_dir)
        return
    for path in _iter_source_files(col_dir):
        kind = path.suffix.lower().lstrip(".")
        yield path, kind


def _index_sources(cfg: dict) -> Iterable[Path]:
    """All source paths that affect index freshness."""
    from scagent.knowledge.parser import parsed_paths

    _, parsed_dir = parsed_paths(cfg)
    for collection in cfg["rag"]["collections"]:
        col_dir = collection_dir(cfg, collection)
        for path, _kind in _iter_collection_files(cfg, collection, col_dir):
            if path.name != "README.md":
                yield path
    if parsed_dir.exists():
        for path in parsed_dir.glob("*.md"):
            yield path


def _index_fresh(cfg: dict, out: Path) -> bool:
    if not out.exists() or out.stat().st_size == 0:
        return False
    newest = 0.0
    for path in _index_sources(cfg):
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

    from scagent.knowledge.parser import chunk_paper_markdown, ensure_papers_parsed

    papers_cfg = (cfg.get("rag") or {}).get("papers") or {}
    if papers_cfg.get("auto_parse", False):
        ensure_papers_parsed(cfg, force=force)

    chunk_size = int(cfg["rag"]["chunk_size"])
    overlap = int(cfg["rag"]["chunk_overlap"])
    chunk_mode = str((cfg.get("rag") or {}).get("chunking") or "semantic")
    paper_chunk_size = int(papers_cfg.get("chunk_size") or chunk_size)
    paper_overlap = int(papers_cfg.get("chunk_overlap") or overlap)

    records: list[dict] = []
    root = Path(cfg.get("_root") or REPO_ROOT)
    for collection in cfg["rag"]["collections"]:
        col_dir = collection_dir(cfg, collection)
        for path, kind in _iter_collection_files(cfg, collection, col_dir):
            if path.name == "README.md":
                continue
            section_chunks = None
            if kind == "pdf":
                units = [_read_pdf(path, cfg=cfg)]
            elif kind == "ipynb":
                units = [_read_ipynb(path)]
            elif kind == "json":
                from scagent.kb import flatten_json_texts

                units = flatten_json_texts(path, collection=collection) or [
                    path.read_text(encoding="utf-8", errors="replace")
                ]
            elif collection == "papers" and kind in {"md", "parsed"}:
                from scagent.knowledge.parser import sanitize_paper_text

                text = sanitize_paper_text(path.read_text(encoding="utf-8", errors="replace"))
                section_chunks = chunk_paper_markdown(
                    text,
                    source="",
                    stem=path.stem,
                    chunk_size=paper_chunk_size,
                    overlap=paper_overlap,
                )
                units = []
            else:
                units = [path.read_text(encoding="utf-8", errors="replace")]
            try:
                rel = str(path.relative_to(root))
            except ValueError:
                rel = str(path.relative_to(knowledge))
            chunks: list[str] = []
            chunk_meta: list[dict] = []
            if section_chunks is not None:
                for sc in section_chunks:
                    sc["source"] = rel
                    chunks.append(sc["text"])
                    chunk_meta.append(sc)
            else:
                for unit in units:
                    if kind == "json":
                        chunks.append(unit)
                        chunk_meta.append({})
                    else:
                        for part in chunk_document(unit, chunk_size, overlap, mode=chunk_mode):
                            chunks.append(part)
                            chunk_meta.append({})
            for i, chunk in enumerate(chunks):
                if not chunk or not str(chunk).strip():
                    continue
                meta = chunk_meta[i] if i < len(chunk_meta) else {}
                rec = {
                    "id": f"{rel}::{i}",
                    "collection": collection,
                    "source": rel,
                    "stem": path.stem,
                    "chunk_index": i,
                    "text": chunk,
                }
                if meta.get("section"):
                    rec["section"] = meta["section"]
                records.append(rec)

    out = index_dir / "chunks.jsonl"
    with out.open("w", encoding="utf-8") as f:
        for rec in records:
            f.write(json.dumps(rec, ensure_ascii=False) + "\n")
    (index_dir / "meta.json").write_text(
        json.dumps({"n_chunks": len(records), "chunking": chunk_mode}, ensure_ascii=False, indent=2),
        encoding="utf-8",
    )
    return out
