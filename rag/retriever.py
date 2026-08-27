from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

from rank_bm25 import BM25Okapi

from rag.ingest import ingest, tokenize
from scagent.config import load_config, resolve_path


def index_path(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    return resolve_path(cfg, "index") / "chunks.jsonl"


def ensure_index(cfg: dict | None = None) -> Path:
    path = index_path(cfg)
    if not path.exists() or path.stat().st_size == 0:
        ingest(cfg)
    return index_path(cfg)


def load_chunks(cfg: dict | None = None, collection: str | None = None) -> list[dict]:
    path = ensure_index(cfg)
    chunks: list[dict] = []
    if not path.exists():
        return chunks
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if collection and rec.get("collection") != collection:
            continue
        chunks.append(rec)
    return chunks


@lru_cache(maxsize=8)
def _bm25_bundle(collection: str | None, index_mtime: float) -> tuple[BM25Okapi, tuple[dict, ...]]:
    del index_mtime
    chunks = tuple(load_chunks(collection=collection))
    corpus = [tokenize(c["text"]) for c in chunks]
    if not corpus:
        corpus = [[]]
        chunks = tuple()
    return BM25Okapi(corpus), chunks


def retrieve(
    query: str,
    *,
    top_k: int | None = None,
    collection: str | None = None,
    cfg: dict | None = None,
    collections: list[str] | None = None,
) -> list[dict]:
    cfg = cfg or load_config()
    top_k = top_k or int(cfg["rag"]["top_k"])
    weights = {"papers": 1.0, "methods": 0.85, "markers": 1.25, "best_practices": 1.15}
    if collections:
        merged: list[dict] = []
        for col in collections:
            for hit in retrieve(query, top_k=top_k, collection=col, cfg=cfg):
                hit = dict(hit)
                hit["score"] = float(hit.get("score") or 0) * weights.get(col, 1.0)
                merged.append(hit)
        merged.sort(key=lambda h: h["score"], reverse=True)
        return merged[:top_k]
    default = collection or cfg["rag"].get("default_collection") or "papers"
    path = ensure_index(cfg)
    mtime = path.stat().st_mtime if path.exists() else 0.0
    bm25, chunks = _bm25_bundle(default, mtime)
    if not chunks:
        return []
    scores = bm25.get_scores(tokenize(query))
    ranked = sorted(range(len(chunks)), key=lambda i: float(scores[i]), reverse=True)
    out: list[dict] = []
    for i in ranked[:top_k]:
        rec = dict(chunks[i])
        rec["score"] = float(scores[i])
        if rec["score"] <= 0:
            continue
        out.append(rec)
    return out


def format_hits(hits: list[dict]) -> str:
    if not hits:
        return "（知识库未检索到相关段落。可将 PDF/Markdown 放入 knowledge/papers 后运行 scagent ingest。）"
    parts = []
    for i, h in enumerate(hits, 1):
        parts.append(
            f"[{i}] {h['source']} (score={h['score']:.3f})\n{h['text'].strip()}"
        )
    return "\n\n".join(parts)
