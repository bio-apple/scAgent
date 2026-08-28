"""Causal evidence chain for biological claims: markers + pathway p-value + PubMed DOI."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from scagent.config import REPO_ROOT

CATALOG_PATH = REPO_ROOT / "knowledge" / "disease_signature" / "cell_states.json"
SKIP_LABELS = {"", "unknown", "unassigned", "unvalidated", "mixed", "none", "nan"}
DOI_RE = re.compile(r"10\.\d{4,9}/[-._;()/:A-Z0-9]+", re.I)


def load_state_catalog(path: Path | None = None) -> list[dict[str, Any]]:
    p = Path(path) if path else CATALOG_PATH
    if not p.is_file():
        return []
    data = json.loads(p.read_text(encoding="utf-8"))
    return list(data.get("states") or [])


def match_cell_state(label: str | None, catalog: list[dict] | None = None) -> dict[str, Any] | None:
    raw = str(label or "").strip()
    if raw.lower() in SKIP_LABELS:
        return None
    states = catalog if catalog is not None else load_state_catalog()
    key = raw.lower()
    for st in states:
        aliases = [str(st.get("name") or "").lower(), *[str(a).lower() for a in (st.get("aliases") or [])]]
        if key in aliases:
            return st
    for st in states:
        name = str(st.get("name") or "").lower()
        if name and name in key:
            return st
    return None


def _doi(value: str | None) -> str | None:
    if not value:
        return None
    m = DOI_RE.search(str(value))
    return m.group(0).rstrip(").,;") if m else None


def validate_claim(claim: dict | None) -> dict[str, Any]:
    """A state assertion is allowed only with markers + pathway p-value + DOI/PMID."""
    claim = dict(claim or {})
    issues: list[str] = []
    genes = []
    for item in claim.get("markers") or []:
        if isinstance(item, str) and item.strip():
            genes.append(item.strip().upper())
        elif isinstance(item, dict) and str(item.get("gene") or "").strip():
            genes.append(str(item.get("gene")).strip().upper())
    genes = list(dict.fromkeys(genes))
    if len(genes) < 2:
        issues.append("需要 ≥2 个 marker/checkpoint 基因（例如 PDCD1 + HAVCR2）")
    pw = claim.get("pathway") or {}
    pid = str(pw.get("id") or pw.get("term") or "").strip()
    pval = pw.get("pval")
    if pval is None:
        pval = pw.get("p_value")
    try:
        pval_f = float(pval) if pval is not None else None
    except (TypeError, ValueError):
        pval_f = None
    if not pid or pval_f is None:
        issues.append("需要通路/GO 编号以及 p-value")
    cites = claim.get("citations") or []
    has_lit = False
    for c in cites:
        if not isinstance(c, dict):
            continue
        if _doi(c.get("doi") or "") or str(c.get("pmid") or "").isdigit():
            has_lit = True
            break
    if not has_lit:
        issues.append("需要 PubMed PMID 或 DOI；禁止编造文献")
    ok = not issues
    return {"ok": ok, "issues": issues, "n_markers": len(genes), "pathway_id": pid, "pval": pval_f}


def _marker_items(observed: list[str], expected: list[str]) -> list[dict[str, str]]:
    items = []
    seen = set()
    for g in list(observed) + list(expected):
        u = str(g).strip().upper()
        if not u or u in seen:
            continue
        seen.add(u)
        kind = "observed" if u in {x.upper() for x in observed} else "canonical"
        items.append({"gene": u, "source": kind})
    return items


def _pick_pathway(state: dict | None, enrichment: dict, genes: list[str]) -> dict[str, Any] | None:
    wanted = []
    for pw in (state or {}).get("pathways") or []:
        wanted.append(str(pw.get("id") or "").upper())
    terms = list((enrichment or {}).get("terms") or [])
    for row in terms:
        tid = str(row.get("term") or row.get("id") or "").upper()
        if wanted and tid in wanted:
            return {
                "id": row.get("term") or row.get("id"),
                "name": row.get("name") or (state or {}).get("pathways", [{}])[0].get("name"),
                "pval": row.get("pval"),
                "fdr": row.get("fdr"),
                "overlap": row.get("overlap"),
                "method": row.get("method") or (enrichment or {}).get("engine"),
            }
    if genes:
        from scagent.enrich import GO_SETS, ora

        rows = ora(genes, gene_sets=GO_SETS, min_overlap=1)
        for row in rows:
            tid = str(row.get("term") or "").upper()
            if wanted and tid not in wanted:
                continue
            name = None
            for pw in (state or {}).get("pathways") or []:
                if str(pw.get("id") or "").upper() == tid:
                    name = pw.get("name")
            return {
                "id": row.get("term"),
                "name": name,
                "pval": row.get("pval"),
                "fdr": row.get("fdr"),
                "overlap": row.get("overlap"),
                "method": row.get("method"),
            }
        if rows:
            row = rows[0]
            return {
                "id": row.get("term"),
                "name": None,
                "pval": row.get("pval"),
                "fdr": row.get("fdr"),
                "overlap": row.get("overlap"),
                "method": row.get("method"),
            }
    if state and (state.get("pathways") or []):
        pw0 = state["pathways"][0]
        return {"id": pw0.get("id"), "name": pw0.get("name"), "pval": None, "fdr": None}
    return None


def assemble_claims(
    annotations: list[dict] | None,
    enrichment: dict | None = None,
    *,
    catalog: list[dict] | None = None,
) -> dict[str, Any]:
    catalog = catalog if catalog is not None else load_state_catalog()
    enrichment = enrichment or {}
    claims: list[dict[str, Any]] = []
    for row in annotations or []:
        label = str(row.get("fused") or row.get("marker_label") or row.get("label") or "").strip()
        if label.lower() in SKIP_LABELS:
            continue
        st = match_cell_state(label, catalog)
        observed = [str(g) for g in (row.get("positive") or []) if g]
        expected = list((st or {}).get("markers") or [])
        markers = _marker_items(observed, expected)
        pathway = _pick_pathway(st, enrichment, observed or expected)
        citations = [dict(c) for c in ((st or {}).get("citations") or [])]
        assertion = f"Cluster {row.get('cluster')} 处于 {label} 状态"
        claim = {
            "cluster": row.get("cluster"),
            "label": label,
            "assertion": assertion,
            "markers": markers,
            "pathway": pathway or {},
            "citations": citations,
            "note": "支持性证据链，不是干预因果。",
        }
        check = validate_claim(claim)
        claim["ok"] = check["ok"]
        claim["issues"] = check["issues"]
        if st is None and not citations:
            claim["issues"] = list(claim["issues"]) + ["本地 cell_states 无匹配状态，禁止编造 DOI"]
            claim["ok"] = False
        claims.append(claim)
    n_ok = sum(1 for c in claims if c.get("ok"))
    return {
        "format": "evidence-chain-v1",
        "n_claims": len(claims),
        "n_ok": n_ok,
        "all_ok": bool(claims) and n_ok == len(claims),
        "claims": claims,
    }


def render_evidence_markdown(payload: dict | None, *, lang: str = "zh") -> str:
    payload = payload or {}
    claims = payload.get("claims") or []
    zh = lang != "en"
    if not claims:
        return ("（无细胞状态断言，无需证据链。）" if zh else "(No cell-state assertions.)") + "\n"
    lines = [
        "每条生物学断言必须同时有 **marker**、**通路 p 值**、**PubMed DOI/PMID**。缺一不可。"
        if zh
        else "Each biological assertion needs **markers**, a **pathway p-value**, and a **PubMed DOI/PMID**.",
        "",
    ]
    for c in claims:
        status = "成立" if c.get("ok") else "证据不足，不得写成结论"
        if not zh:
            status = "supported" if c.get("ok") else "unsupported — do not report as fact"
        lines.append(f"### {c.get('assertion')}")
        lines.append("")
        lines.append(f"- {'状态' if zh else 'status'}: {status}")
        genes = []
        for m in c.get("markers") or []:
            if isinstance(m, dict):
                genes.append(f"{m.get('gene')} ({m.get('source')})")
            else:
                genes.append(str(m))
        lines.append(f"- markers: {', '.join(genes) or ('无' if zh else 'none')}")
        pw = c.get("pathway") or {}
        ptxt = "n/a"
        if pw.get("pval") is not None:
            ptxt = f"p={pw.get('pval')}"
            if pw.get("fdr") is not None:
                ptxt += f", FDR={pw.get('fdr')}"
        lines.append(f"- pathway: `{pw.get('id') or 'n/a'}` {pw.get('name') or ''} ({ptxt})")
        cites = c.get("citations") or []
        if cites:
            for cit in cites:
                doi = cit.get("doi") or ""
                pmid = cit.get("pmid") or ""
                title = cit.get("title") or ""
                url = f"https://doi.org/{doi}" if doi else (f"https://pubmed.ncbi.nlm.nih.gov/{pmid}" if pmid else "")
                lines.append(f"- citation: {title} DOI `{doi}` PMID {pmid} {url}".rstrip())
        else:
            lines.append("- citation: （无；不得编造）" if zh else "- citation: none; do not invent")
        if c.get("issues"):
            lines.append("- issues: " + "；".join(c["issues"]))
        lines.append(f"- {c.get('note') or ''}")
        lines.append("")
    return "\n".join(lines).rstrip() + "\n"


def write_evidence_chains(workspace: str | Path, *, tissue: str | None = None) -> dict[str, Any]:
    del tissue
    ws = Path(workspace)
    annotations: list[dict] = []
    ap = ws / "annotation_evidence.json"
    if ap.is_file():
        raw = json.loads(ap.read_text(encoding="utf-8"))
        annotations = raw if isinstance(raw, list) else list(raw.get("clusters") or raw.get("rows") or [])
    enrichment: dict = {}
    ep = ws / "pathway_enrichment.json"
    if ep.is_file():
        enrichment = json.loads(ep.read_text(encoding="utf-8"))
    payload = assemble_claims(annotations, enrichment)
    (ws / "evidence_chains.json").write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return payload
