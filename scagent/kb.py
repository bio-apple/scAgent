"""Structured knowledge base: Cell Ontology, markers, pathways, disease signatures, tissues."""

from __future__ import annotations

import json
from pathlib import Path

from scagent.config import REPO_ROOT

KB_ROOT = REPO_ROOT / "knowledge"
STRUCTURED_COLLECTIONS = (
    "cell_ontology",
    "marker_db",
    "pathway",
    "disease_signature",
    "tissue_reference",
)
COLLECTION_ALIASES = {"markers": "marker_db"}


def _dirs(collection: str) -> Path:
    name = COLLECTION_ALIASES.get(collection, collection)
    return KB_ROOT / name


def _format_kv(rec: dict, source: str) -> str:
    lines = []
    title = rec.get("name") or rec.get("id") or rec.get("term") or "record"
    rid = rec.get("id") or rec.get("cl_id") or ""
    lines.append(f"# {title}" + (f" ({rid})" if rid else ""))
    if source:
        lines.append(f"source: {source}")
    for key in ("synonyms", "aliases", "tissues", "expected_types", "lineage", "markers", "positive", "negative", "genes"):
        val = rec.get(key)
        if not val:
            continue
        if isinstance(val, list):
            lines.append(f"{key}: {', '.join(str(x) for x in val)}")
        else:
            lines.append(f"{key}: {val}")
    for key in ("parent", "atlas", "qc_note", "collection", "tissue"):
        if rec.get(key):
            lines.append(f"{key}: {rec[key]}")
    pw = rec.get("pathways") or []
    if pw:
        bits = []
        for p in pw:
            if isinstance(p, dict):
                bits.append(str(p.get("id") or p.get("name") or p))
            else:
                bits.append(str(p))
        lines.append(f"pathways: {', '.join(bits)}")
    cites = rec.get("citations") or []
    if cites:
        ids = []
        for c in cites:
            if isinstance(c, dict):
                ids.append(str(c.get("pmid") or c.get("doi") or c.get("title") or ""))
        lines.append("citations: " + "; ".join(x for x in ids if x))
    return "\n".join(lines)


def records_from_json(data: dict, *, collection: str, source_file: str) -> list[dict]:
    src = str(data.get("source") or collection)
    out: list[dict] = []

    def add(rec: dict, kind: str) -> None:
        text = _format_kv(rec, src)
        name = str(rec.get("name") or rec.get("id") or "")
        out.append(
            {
                "collection": collection,
                "source": source_file,
                "stem": Path(source_file).stem,
                "name": name,
                "id": rec.get("id") or rec.get("cl_id"),
                "kind": kind,
                "record": rec,
                "text": text,
            }
        )

    if isinstance(data.get("records"), list):
        for rec in data["records"]:
            if isinstance(rec, dict):
                add(rec, "ontology")
    if isinstance(data.get("states"), list):
        for rec in data["states"]:
            if isinstance(rec, dict):
                add(rec, "disease_signature")
    if isinstance(data.get("sets"), list):
        for rec in data["sets"]:
            if isinstance(rec, dict):
                add(rec, "pathway")
    tissues = data.get("tissues")
    if isinstance(tissues, dict):
        for tissue, rows in tissues.items():
            for rec in rows or []:
                if isinstance(rec, dict):
                    row = dict(rec)
                    row.setdefault("tissue", tissue)
                    add(row, "marker")
    elif isinstance(tissues, list):
        for rec in tissues:
            if isinstance(rec, dict):
                add(rec, "tissue")
    return out


def load_structured(collections: list[str] | None = None) -> list[dict]:
    cols = [COLLECTION_ALIASES.get(c, c) for c in (collections or STRUCTURED_COLLECTIONS)]
    seen: set[str] = set()
    out: list[dict] = []
    for col in cols:
        if col in seen:
            continue
        seen.add(col)
        d = _dirs(col)
        if not d.is_dir():
            continue
        for path in sorted(d.glob("*.json")):
            try:
                data = json.loads(path.read_text(encoding="utf-8"))
            except (OSError, json.JSONDecodeError):
                continue
            if not isinstance(data, dict):
                continue
            try:
                rel = str(path.relative_to(REPO_ROOT))
            except ValueError:
                rel = str(path)
            out.extend(records_from_json(data, collection=col, source_file=rel))
    return out


def flatten_json_texts(path: Path, *, collection: str) -> list[str]:
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return []
    if not isinstance(data, dict):
        return [path.read_text(encoding="utf-8", errors="replace")]
    rel = str(path)
    try:
        rel = str(path.relative_to(REPO_ROOT))
    except ValueError:
        pass
    recs = records_from_json(data, collection=collection, source_file=rel)
    return [r["text"] for r in recs if r.get("text")]


def lookup_structured(
    query: str,
    *,
    collections: list[str] | None = None,
    tissue: str | None = None,
    top_k: int = 8,
) -> list[dict]:
    """Keyword lookup over JSON catalogs (structured records, not prompt prose)."""
    q = (query or "").strip().lower()
    tokens = [t for t in q.replace(",", " ").replace(":", " ").split() if len(t) > 1]
    tissue_key = (tissue or "").strip().lower()
    scored: list[tuple[float, dict]] = []
    for rec in load_structured(collections):
        blob = " ".join(
            [
                rec.get("name") or "",
                str(rec.get("id") or ""),
                rec.get("text") or "",
                str((rec.get("record") or {}).get("tissue") or ""),
            ]
        ).lower()
        score = 0.0
        if q and q in blob:
            score += 3.0
        for t in tokens:
            if t in blob:
                score += 1.0
        rec_tissue = str((rec.get("record") or {}).get("tissue") or "").lower()
        tissues = (rec.get("record") or {}).get("tissues") or []
        tissue_hit = rec_tissue == tissue_key or tissue_key in {str(x).lower() for x in tissues} or tissue_key in blob
        if tissue_key and tissue_hit:
            score += 1.5
        elif tissue_key and rec.get("kind") in {"marker", "ontology", "tissue"} and not tissue_hit:
            score *= 0.35
        if score <= 0:
            continue
        hit = dict(rec)
        hit["score"] = score
        hit["retrieval"] = "structured"
        scored.append((score, hit))
    scored.sort(key=lambda kv: kv[0], reverse=True)
    return [h for _, h in scored[:top_k]]


def format_records(records: list[dict]) -> str:
    if not records:
        return "（结构化知识库无匹配记录。）"
    parts = []
    for i, h in enumerate(records, 1):
        col = h.get("collection") or ""
        src = h.get("source") or ""
        parts.append(f"[{i}] {src} [{col}]\n{(h.get('text') or '').strip()}")
    return "\n\n".join(parts)


def gene_sets_from_kb() -> dict[str, tuple[str, ...]]:
    sets: dict[str, tuple[str, ...]] = {}
    for rec in load_structured(["pathway"]):
        row = rec.get("record") or {}
        name = str(row.get("id") or rec.get("name") or "")
        genes = tuple(str(g).upper() for g in (row.get("genes") or []) if g)
        if name and genes:
            sets[name] = genes
    return sets
