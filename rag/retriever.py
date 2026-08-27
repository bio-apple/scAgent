from __future__ import annotations

import json
import pickle
from functools import lru_cache
from pathlib import Path

from rank_bm25 import BM25Okapi

from rag.embed import cosine, embed_backend, embed_texts
from rag.ingest import ingest, tokenize
from rag.synonyms import expand_query
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
def _bm25_bundle(collection: str | None, index_mtime: float) -> tuple[BM25Okapi | None, tuple[dict, ...]]:
    del index_mtime
    chunks = tuple(load_chunks(collection=collection))
    corpus = [tokenize(c["text"]) for c in chunks]
    if not chunks or not any(corpus):
        return None, tuple()
    return BM25Okapi(corpus), chunks


@lru_cache(maxsize=8)
def _vector_bundle(collection: str | None, index_mtime: float, backend: str) -> tuple[tuple[list[float], ...], tuple[dict, ...]]:
    del index_mtime, backend
    chunks = tuple(load_chunks(collection=collection))
    if not chunks:
        return tuple(), tuple()
    vecs = _cached_embeddings(chunks)
    return tuple(vecs), chunks


def _cached_embeddings(chunks: tuple[dict, ...]) -> list[list[float]]:
    cfg = load_config()
    backend = embed_backend(cfg)
    idx = index_path(cfg)
    col = (chunks[0].get("collection") if chunks else None) or "all"
    cache = idx.parent / f"embeddings.{backend}.{col}.pkl"
    if cache.exists() and idx.exists() and cache.stat().st_mtime >= idx.stat().st_mtime:
        try:
            data = pickle.loads(cache.read_bytes())
            if data.get("n") == len(chunks) and data.get("backend") == backend and data.get("collection") == col:
                return data["vecs"]
        except Exception:
            pass
    vecs = embed_texts([c["text"] for c in chunks], cfg)
    try:
        cache.write_bytes(pickle.dumps({"n": len(chunks), "backend": backend, "collection": col, "vecs": vecs}, protocol=4))
    except OSError:
        pass
    return vecs


def clear_retrieve_cache() -> None:
    _bm25_bundle.cache_clear()
    _vector_bundle.cache_clear()


def _rrf_merge(ranked_lists: list[list[int]], k: int) -> list[tuple[int, float]]:
    scores: dict[int, float] = {}
    for ranked in ranked_lists:
        for rank, idx in enumerate(ranked, 1):
            scores[idx] = scores.get(idx, 0.0) + 1.0 / (k + rank)
    return sorted(scores.items(), key=lambda kv: kv[1], reverse=True)


def _overlap(query: str, text: str) -> float:
    q = set(tokenize(query))
    d = set(tokenize(text))
    if not q:
        return 0.0
    return len(q & d) / len(q)


def retrieve(
    query: str,
    *,
    top_k: int | None = None,
    collection: str | None = None,
    cfg: dict | None = None,
    collections: list[str] | None = None,
) -> list[dict]:
    cfg = cfg or load_config()
    rag = cfg.get("rag") or {}
    top_k = top_k or int(rag.get("top_k") or 6)
    weights = {"papers": 1.0, "methods": 0.85, "markers": 1.25, "best_practices": 1.15, "sops": 1.3}
    if collections:
        merged: list[dict] = []
        for col in collections:
            for hit in retrieve(query, top_k=top_k, collection=col, cfg=cfg):
                hit = dict(hit)
                hit["score"] = float(hit.get("score") or 0) * weights.get(col, 1.0)
                merged.append(hit)
        merged.sort(key=lambda h: h["score"], reverse=True)
        return merged[:top_k]
    default = collection or rag.get("default_collection") or "papers"
    path = ensure_index(cfg)
    mtime = path.stat().st_mtime if path.exists() else 0.0
    bm25, chunks = _bm25_bundle(default, mtime)
    if not chunks:
        return []
    expanded = expand_query(query)
    candidate_k = min(int(rag.get("candidate_k") or 24), len(chunks))
    rrf_k = int(rag.get("rrf_k") or 60)
    tokens = tokenize(expanded)
    bm_scores = bm25.get_scores(tokens)
    bm_rank = sorted(range(len(chunks)), key=lambda i: float(bm_scores[i]), reverse=True)
    mode = str(rag.get("retrieval") or "hybrid").lower()
    if mode != "hybrid":
        out: list[dict] = []
        for i in bm_rank[:top_k]:
            if float(bm_scores[i]) <= 0:
                continue
            rec = dict(chunks[i])
            rec["score"] = float(bm_scores[i])
            rec["retrieval"] = "bm25"
            out.append(rec)
        return out
    vecs, vchunks = _vector_bundle(default, mtime, embed_backend(cfg))
    qvec = embed_texts([expanded], cfg)[0] if vecs else []
    vec_scores = [cosine(qvec, v) if qvec else 0.0 for v in vecs]
    vec_rank = sorted(range(len(vchunks)), key=lambda i: vec_scores[i], reverse=True)
    fused = _rrf_merge([bm_rank[:candidate_k], vec_rank[:candidate_k]], rrf_k)
    reranked = []
    for idx, rrf in fused:
        rec = dict(chunks[idx])
        ov = _overlap(expanded, rec.get("text") or "")
        rec["score"] = float(rrf) + 0.15 * ov
        rec["retrieval"] = "hybrid"
        rec["score_bm25"] = float(bm_scores[idx]) if idx < len(bm_scores) else 0.0
        rec["score_vec"] = float(vec_scores[idx]) if idx < len(vec_scores) else 0.0
        reranked.append(rec)
    reranked.sort(key=lambda h: h["score"], reverse=True)
    out = [h for h in reranked if h["score"] > 0][:top_k]
    return out


def format_hits(hits: list[dict]) -> str:
    if not hits:
        return "（知识库未检索到相关段落。scagent update-kb 拉取 sc-best-practices；scagent add-doc <path> 加入实验室 SOP；或把 PDF/Markdown 放入 knowledge/papers 后 scagent ingest。）"
    parts = []
    for i, h in enumerate(hits, 1):
        parts.append(f"[{i}] {h['source']} (score={h['score']:.3f})\n{h['text'].strip()}")
    return "\n\n".join(parts)
