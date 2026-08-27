"""Dense embeddings: sentence-transformers when installed, else stable hashing vectors."""

from __future__ import annotations

import hashlib
import math
from functools import lru_cache

from rag.ingest import tokenize
from scagent.config import load_config
from scagent.logutil import get_logger

log = get_logger("rag.embed")

HASH_DIM = 384


def _hashing_vec(text: str, dim: int = HASH_DIM) -> list[float]:
    acc = [0.0] * dim
    for tok in tokenize(text):
        digest = hashlib.md5(tok.encode("utf-8")).hexdigest()
        h1 = int(digest[:8], 16)
        h2 = int(digest[8:16], 16)
        acc[h1 % dim] += 1.0
        acc[h2 % dim] += 0.5
    n = math.sqrt(sum(x * x for x in acc)) or 1.0
    return [x / n for x in acc]


@lru_cache(maxsize=1)
def _sbert(model_name: str):
    from sentence_transformers import SentenceTransformer

    log.info("loading sentence-transformers %s", model_name)
    return SentenceTransformer(model_name)


def embed_backend(cfg: dict | None = None) -> str:
    cfg = cfg or load_config()
    rag = cfg.get("rag") or {}
    requested = str(rag.get("embed_backend") or "auto").lower()
    if requested == "hashing":
        return "hashing"
    if requested in {"sbert", "auto"}:
        try:
            import sentence_transformers  # noqa: F401

            return "sbert"
        except Exception:
            if requested == "sbert":
                log.warning("sentence-transformers 未安装，改用 hashing 向量")
            return "hashing"
    return "hashing"


def embed_texts(texts: list[str], cfg: dict | None = None) -> list[list[float]]:
    cfg = cfg or load_config()
    backend = embed_backend(cfg)
    if backend == "sbert":
        try:
            model = str((cfg.get("rag") or {}).get("embed_model") or "sentence-transformers/all-MiniLM-L6-v2")
            vecs = _sbert(model).encode(list(texts), normalize_embeddings=True, show_progress_bar=False)
            return [list(map(float, v)) for v in vecs]
        except Exception as exc:
            log.warning("sbert embed failed, hashing fallback: %s", exc)
    return [_hashing_vec(t) for t in texts]


def cosine(a: list[float], b: list[float]) -> float:
    n = min(len(a), len(b))
    if n == 0:
        return 0.0
    return float(sum(a[i] * b[i] for i in range(n)))
