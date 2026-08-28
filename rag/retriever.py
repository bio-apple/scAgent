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

COLLECTION_WEIGHTS = {
    "sops": 1.3,
    "cell_ontology": 1.28,
    "marker_db": 1.28,
    "markers": 1.28,
    "disease_signature": 1.26,
    "best_practices": 1.2,
    "pathway": 1.22,
    "tissue_reference": 1.18,
    "papers": 1.0,
    "methods": 0.9,
    "upstream": 0.65,
}
AGENT_COLLECTIONS = ("best_practices", "papers", "methods", "sops", "upstream")
KB_COLLECTIONS = (
    "cell_ontology",
    "marker_db",
    "pathway",
    "disease_signature",
    "tissue_reference",
)
_COLLECTION_ALIASES = {"markers": "marker_db"}
_PRACTICE_BOOST = 1.4


def index_path(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    return resolve_path(cfg, "index") / "chunks.jsonl"


def ensure_index(cfg: dict | None = None) -> Path:
    path = index_path(cfg)
    if not path.exists() or path.stat().st_size == 0:
        ingest(cfg)
    return index_path(cfg)


def _canon_collection(name: str | None) -> str | None:
    if not name:
        return name
    return _COLLECTION_ALIASES.get(name, name)


def _allow_set(collection: str | None, collections: list[str] | tuple[str, ...] | None) -> set[str] | None:
    if collections:
        return {_canon_collection(c) or c for c in collections}
    if collection:
        col = _canon_collection(collection) or collection
        return {col, collection}
    return None


def load_chunks(
    cfg: dict | None = None,
    collection: str | None = None,
    collections: list[str] | tuple[str, ...] | None = None,
) -> list[dict]:
    path = ensure_index(cfg)
    chunks: list[dict] = []
    if not path.exists():
        return chunks
    allow = _allow_set(collection, collections)
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        rec = json.loads(line)
        if allow is not None and rec.get("collection") not in allow:
            continue
        chunks.append(rec)
    return chunks


def _cache_key(collection: str | None, collections: list[str] | tuple[str, ...] | None) -> str:
    if collections:
        return ",".join(sorted(collections))
    if collection:
        return collection
    return "__all__"


@lru_cache(maxsize=16)
def _bm25_bundle(cache_key: str, index_mtime: float) -> tuple[BM25Okapi | None, tuple[dict, ...]]:
    del index_mtime
    cols = None if cache_key in {"__all__", ""} else tuple(cache_key.split(","))
    chunks = tuple(load_chunks(collection=None, collections=cols))
    corpus = [tokenize(c["text"]) for c in chunks]
    if not chunks or not any(corpus):
        return None, tuple()
    return BM25Okapi(corpus), chunks


@lru_cache(maxsize=16)
def _vector_bundle(cache_key: str, index_mtime: float, backend: str) -> tuple[tuple[list[float], ...], tuple[dict, ...]]:
    del index_mtime, backend
    cols = None if cache_key in {"__all__", ""} else tuple(cache_key.split(","))
    chunks = tuple(load_chunks(collection=None, collections=cols))
    if not chunks:
        return tuple(), tuple()
    vecs = _cached_embeddings(chunks, cache_key)
    return tuple(vecs), chunks


def _cached_embeddings(chunks: tuple[dict, ...], cache_key: str) -> list[list[float]]:
    cfg = load_config()
    backend = embed_backend(cfg)
    idx = index_path(cfg)
    cache = idx.parent / f"embeddings.{backend}.{cache_key.replace(',', '_')}.pkl"
    if cache.exists() and idx.exists() and cache.stat().st_mtime >= idx.stat().st_mtime:
        try:
            data = pickle.loads(cache.read_bytes())
            if data.get("n") == len(chunks) and data.get("backend") == backend and data.get("collection") == cache_key:
                return data["vecs"]
        except Exception:
            pass
    vecs = embed_texts([c["text"] for c in chunks], cfg)
    try:
        cache.write_bytes(
            pickle.dumps({"n": len(chunks), "backend": backend, "collection": cache_key, "vecs": vecs}, protocol=4)
        )
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


def _stem(rec: dict) -> str:
    if rec.get("stem"):
        return str(rec["stem"]).lower()
    return Path(str(rec.get("source") or "")).stem.lower()


def _hit_weight(rec: dict, boost_stems: set[str]) -> float:
    w = COLLECTION_WEIGHTS.get(str(rec.get("collection") or ""), 1.0)
    if _stem(rec) in boost_stems:
        w *= _PRACTICE_BOOST
    return w


def _resolve_scope(
    rag: dict,
    collection: str | None,
    collections: list[str] | None,
) -> tuple[str | None, list[str] | None]:
    if collections:
        return None, [_canon_collection(c) or c for c in collections]
    if collection:
        return _canon_collection(collection), None
    default = rag.get("default_collection")
    if default in (None, "", "fused", "all"):
        return None, list(rag.get("collections") or [])
    return str(default), None


def retrieve(
    query: str,
    *,
    top_k: int | None = None,
    collection: str | None = None,
    cfg: dict | None = None,
    collections: list[str] | None = None,
    boost_stems: list[str] | None = None,
) -> list[dict]:
    cfg = cfg or load_config()
    rag = cfg.get("rag") or {}
    top_k = top_k or int(rag.get("top_k") or 6)
    collection, collections = _resolve_scope(rag, collection, collections)
    path = ensure_index(cfg)
    mtime = path.stat().st_mtime if path.exists() else 0.0
    key = _cache_key(collection, collections)
    bm25, chunks = _bm25_bundle(key, mtime)
    if not chunks:
        return []
    expanded = expand_query(query)
    candidate_k = min(int(rag.get("candidate_k") or 24), len(chunks))
    rrf_k = int(rag.get("rrf_k") or 60)
    tokens = tokenize(expanded)
    bm_scores = bm25.get_scores(tokens)
    bm_rank = sorted(range(len(chunks)), key=lambda i: float(bm_scores[i]), reverse=True)
    boost = {s.strip().lower().removesuffix(".md") for s in (boost_stems or []) if s}
    mode = str(rag.get("retrieval") or "hybrid").lower()

    def _finalize(idx: int, score: float, retrieval: str, extra: dict | None = None) -> dict:
        rec = dict(chunks[idx])
        rec["score"] = float(score) * _hit_weight(rec, boost)
        rec["retrieval"] = retrieval
        if extra:
            rec.update(extra)
        return rec

    if mode != "hybrid":
        out: list[dict] = []
        for i in bm_rank[:top_k]:
            if float(bm_scores[i]) <= 0:
                continue
            out.append(_finalize(i, float(bm_scores[i]), "bm25"))
        out.sort(key=lambda h: h["score"], reverse=True)
        return _dedup_hits(out, top_k)
    vecs, vchunks = _vector_bundle(key, mtime, embed_backend(cfg))
    qvec = embed_texts([expanded], cfg)[0] if vecs else []
    vec_scores = [cosine(qvec, v) if qvec else 0.0 for v in vecs]
    vec_rank = sorted(range(len(vchunks)), key=lambda i: vec_scores[i], reverse=True)
    fused = _rrf_merge([bm_rank[:candidate_k], vec_rank[:candidate_k]], rrf_k)
    reranked = []
    for idx, rrf in fused:
        ov = _overlap(expanded, chunks[idx].get("text") or "")
        rec = _finalize(
            idx,
            float(rrf) + 0.15 * ov,
            "hybrid",
            {
                "score_bm25": float(bm_scores[idx]) if idx < len(bm_scores) else 0.0,
                "score_vec": float(vec_scores[idx]) if idx < len(vec_scores) else 0.0,
            },
        )
        reranked.append(rec)
    reranked.sort(key=lambda h: h["score"], reverse=True)
    return _dedup_hits([h for h in reranked if h["score"] > 0], top_k)


def _dedup_hits(hits: list[dict], top_k: int) -> list[dict]:
    seen: set[tuple[str, str]] = set()
    out: list[dict] = []
    for h in hits:
        key = (str(h.get("source") or ""), (h.get("text") or "")[:96])
        if key in seen:
            continue
        seen.add(key)
        out.append(h)
        if len(out) >= top_k:
            break
    return out


def retrieve_fused(
    query: str,
    *,
    phase: str | None = None,
    route: list[str] | None = None,
    intents: list[str] | None = None,
    user_query: str | None = None,
    collections: list[str] | None = None,
    include_markers: bool = False,
    include_kb: bool | None = None,
    tissue: str | None = None,
    top_k: int | None = None,
    cfg: dict | None = None,
    ensure_practices: bool = True,
) -> list[dict]:
    """Hybrid retrieve across knowledge/ collections, boosting route-matched SOP files."""
    from scagent.best_practices_loader import load_practice_text, practices_for_phase, practices_for_route
    from scagent.kb import lookup_structured

    cfg = cfg or load_config()
    hint = user_query or query
    names = (
        practices_for_phase(phase, route=route, query=hint)
        if phase
        else practices_for_route(route, intents, hint)
    )
    cols = list(collections or AGENT_COLLECTIONS)
    want_kb = include_kb if include_kb is not None else True
    if include_markers and "marker_db" not in cols:
        cols.append("marker_db")
    if want_kb:
        for col in KB_COLLECTIONS:
            if col not in cols:
                cols.append(col)
    q = query if not names else f"{query} {' '.join(names)}"
    hits = retrieve(q, collections=cols, boost_stems=names, top_k=top_k, cfg=cfg)
    kb_hits = []
    if want_kb:
        if phase == "qc":
            kb_cols: list[str] | None = []
        elif phase == "annotation" or include_markers:
            kb_cols = ["marker_db", "cell_ontology", "tissue_reference"]
        elif phase == "interpret":
            kb_cols = ["pathway", "disease_signature"]
        else:
            kb_cols = None
        if kb_cols is not None and not kb_cols:
            kb_hits = []
        else:
            kb_hits = lookup_structured(hint, collections=kb_cols, tissue=tissue, top_k=4)
            for h in kb_hits:
                h.setdefault("retrieval", "structured")
    if not ensure_practices or not names:
        merged = kb_hits + hits
        want = top_k or int((cfg.get("rag") or {}).get("top_k") or 6)
        return _dedup_hits(merged, want)
    want = top_k or int((cfg.get("rag") or {}).get("top_k") or 6)
    covered = {_stem(h) for h in hits}
    extras: list[dict] = []
    for name in names:
        if name.lower() in covered:
            continue
        extras.append(
            {
                "id": f"practice::{name}",
                "collection": "best_practices",
                "source": f"knowledge/best_practices/{name}.md",
                "stem": name,
                "chunk_index": 0,
                "text": load_practice_text(name, max_chars=1600),
                "score": 1.0,
                "retrieval": "practice_ensure",
            }
        )
        if len(extras) >= 3:
            break
    return _dedup_hits(extras + kb_hits + hits, want)


PAPER_SECTION_WEIGHTS = {
    "methods": 1.35,
    "results": 1.3,
    "abstract": 1.2,
    "introduction": 1.05,
    "discussion": 1.0,
    "conclusion": 1.0,
    "other": 1.0,
}


def search_paper_knowledge(
    query: str,
    *,
    sections: list[str] | tuple[str, ...] | None = None,
    top_k: int | None = None,
    cfg: dict | None = None,
) -> list[dict]:
    """Retrieve from papers collection; boost Methods/Results/Abstract chunks."""
    cfg = cfg or load_config()
    rag = cfg.get("rag") or {}
    papers_cfg = rag.get("papers") or {}
    prefer = [s.lower() for s in (sections or papers_cfg.get("prefer_sections") or ["methods", "results", "abstract"])]
    top_k = top_k or int(rag.get("top_k") or 6)
    hits = retrieve(query, collections=["papers"], top_k=min(top_k * 3, 24), cfg=cfg)
    if prefer:
        filtered = [h for h in hits if str(h.get("section") or "other").lower() in prefer]
        if filtered:
            hits = filtered + [h for h in hits if h not in filtered]
    boosted: list[dict] = []
    for h in hits:
        rec = dict(h)
        sec = str(rec.get("section") or "other").lower()
        rec["score"] = float(rec.get("score") or 0) * PAPER_SECTION_WEIGHTS.get(sec, 1.0)
        boosted.append(rec)
    boosted.sort(key=lambda x: x["score"], reverse=True)
    return _dedup_hits(boosted, top_k)


def format_paper_hits(hits: list[dict]) -> str:
    if not hits:
        return (
            "No paper chunks matched. Run `scagent parse-papers` then `scagent ingest` "
            "after adding PDFs under knowledge/papers/."
        )
    parts = []
    for i, h in enumerate(hits, 1):
        sec = h.get("section") or "body"
        title_bit = h.get("stem") or Path(str(h.get("source") or "")).stem
        parts.append(
            f"[{i}] {title_bit} — {sec} ({h.get('source')}, score={h.get('score', 0):.3f})\n"
            f"{h.get('text', '').strip()}"
        )
    return "\n\n".join(parts)


def format_hits(hits: list[dict]) -> str:
    if not hits:
        return (
            "（知识库未检索到相关段落。scagent update-kb 拉取 sc-best-practices 到 knowledge/upstream；"
            "scagent add-doc <path> 加入实验室 SOP；或把 PDF/Markdown 放入 knowledge/papers 后 scagent ingest。）"
        )
    parts = []
    for i, h in enumerate(hits, 1):
        col = h.get("collection") or ""
        tag = f" [{col}]" if col else ""
        parts.append(f"[{i}] {h['source']}{tag} (score={h['score']:.3f})\n{h['text'].strip()}")
    return "\n\n".join(parts)
