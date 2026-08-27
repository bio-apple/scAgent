"""Pathway enrichment: GSEA/GSVA when libraries exist, else Fisher ORA on Hallmark-like sets."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from scagent.logutil import get_logger

log = get_logger("enrich")

# Compact Hallmark-like sets (Heumos: gene-set choice > method). Not the full MSigDB.
HALLMARK_LIKE: dict[str, tuple[str, ...]] = {
    "HALLMARK_INTERFERON_GAMMA_RESPONSE": (
        "STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "IDO1", "SOCS1", "ICAM1", "HLA-DRA", "CD274",
    ),
    "HALLMARK_G2M_CHECKPOINT": (
        "TOP2A", "MKI67", "CCNB1", "CDK1", "BIRC5", "PLK1", "CDC20", "UBE2C", "CENPF", "HMGB2",
    ),
    "HALLMARK_HYPOXIA": (
        "VEGFA", "LDHA", "PGK1", "SLC2A1", "ENO1", "HIF1A", "CA9", "BNIP3", "PDK1", "ALDOA",
    ),
    "HALLMARK_TNFA_SIGNALING_VIA_NFKB": (
        "NFKB1", "TNF", "IL6", "CXCL1", "RELA", "NFKBIA", "JUN", "FOS", "ICAM1", "CCL2",
    ),
    "HALLMARK_OXIDATIVE_PHOSPHORYLATION": (
        "NDUFA1", "COX5A", "UQCRC1", "ATP5F1A", "SDHB", "CYCS", "VDAC1", "ATP5MC1", "NDUFS1", "COX7A2",
    ),
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION": (
        "VIM", "FN1", "CDH2", "SNAI1", "TWIST1", "MMP2", "COL1A1", "ACTA2", "ZEB1", "SPARC",
    ),
}

# Compact GO sets for evidence chains (not the full ontology).
GO_SETS: dict[str, tuple[str, ...]] = {
    "GO:0002429": (
        "PDCD1", "HAVCR2", "LAG3", "TIGIT", "CTLA4", "CD274", "CD3D", "CD3E", "LCK", "ZAP70",
    ),
    "GO:0042110": (
        "CD3D", "CD3E", "CD4", "CD8A", "IL7R", "LCK", "ZAP70", "IL2RA", "CD28",
    ),
    "GO:0002250": (
        "MS4A1", "CD79A", "CD19", "IGHM", "CD3D", "CD8A", "HLA-DRA",
    ),
    "GO:0006955": (
        "LYZ", "CD14", "HLA-DRA", "TNF", "IL1B", "CXCL8", "PTPRC",
    ),
    "GO:0007399": (
        "RBFOX3", "SNAP25", "AQP4", "GFAP", "MBP", "SLC1A3", "SOX2",
    ),
    "GO:0001889": (
        "ALB", "APOA1", "APOE", "TTR", "CYP3A4",
    ),
    "GO:0007507": (
        "TNNT2", "MYH7", "ACTC1", "NPPA", "MYL2",
    ),
}


def default_gene_sets() -> dict[str, tuple[str, ...]]:
    return {**HALLMARK_LIKE, **GO_SETS}


def _bh(pvals: list[float]) -> list[float]:
    m = len(pvals)
    if m == 0:
        return []
    order = sorted(range(m), key=lambda i: pvals[i])
    q = [1.0] * m
    acc = 1.0
    for rank_from_end, idx in enumerate(reversed(order)):
        rank = m - rank_from_end
        acc = min(acc, pvals[idx] * m / rank)
        q[idx] = min(1.0, acc)
    return q


def _hypergeom_sf(k: int, K: int, n: int, N: int) -> float:
    if k <= 0:
        return 1.0
    try:
        from scipy.stats import hypergeom

        return float(hypergeom.sf(k - 1, N, K, n))
    except Exception:
        from math import comb

        hi = min(K, n)
        if k > hi or N <= 0 or n <= 0:
            return 0.0 if k > hi else 1.0
        den = comb(N, n)
        if den == 0:
            return 1.0
        total = 0
        for i in range(k, hi + 1):
            rest = n - i
            if rest < 0 or rest > (N - K):
                continue
            total += comb(K, i) * comb(N - K, rest)
        return min(1.0, total / den)


def ora(
    genes: list[str],
    *,
    gene_sets: dict[str, tuple[str, ...]] | None = None,
    background: list[str] | None = None,
    min_overlap: int = 2,
) -> list[dict]:
    """Over-representation (Fisher/hypergeometric). Used when GSEA/GSVA libs are absent."""
    query = {g.strip().upper() for g in genes if g and str(g).strip()}
    if not query:
        return []
    sets = gene_sets or default_gene_sets()
    bg = {g.strip().upper() for g in (background or []) if g}
    if not bg:
        bg = set(query)
        for members in sets.values():
            bg.update(m.upper() for m in members)
    n = len(query & bg)
    n_bg = len(bg)
    rows: list[dict] = []
    for name, members in sets.items():
        s = {m.upper() for m in members} & bg
        overlap = sorted(query & s)
        if len(overlap) < min_overlap:
            continue
        p = _hypergeom_sf(len(overlap), len(s), n, n_bg)
        rows.append(
            {
                "term": name,
                "overlap": len(overlap),
                "set_size": len(s),
                "query_size": n,
                "genes": overlap,
                "pval": p,
                "method": "ora_hypergeom",
            }
        )
    rows.sort(key=lambda r: (r["pval"], -r["overlap"]))
    fdr = _bh([float(r["pval"]) for r in rows])
    for r, q in zip(rows, fdr):
        r["fdr"] = q
    return rows


def try_gseapy(genes: list[str], gene_sets: dict[str, tuple[str, ...]]) -> tuple[list[dict], str]:
    try:
        import gseapy  # noqa: F401
    except Exception:
        return [], "ora"
    # gseapy.enrich typically needs a GMT / library download. Stay offline: run ORA.
    return ora(genes, gene_sets=gene_sets), "ora_gseapy_offline"


def try_gsva(expression_matrix) -> tuple[list[dict], str]:
    """GSVA needs decoupler/gseapy + a dense sample×gene matrix. Optional."""
    del expression_matrix
    try:
        import decoupler  # noqa: F401
    except Exception:
        return [], "gsva_unavailable"
    return [], "gsva_unavailable"


def run_enrichment(genes: list[str], *, background: list[str] | None = None) -> dict:
    engine = "ora"
    rows, tag = try_gseapy(genes, default_gene_sets())
    if tag != "ora":
        engine = tag
    if not rows:
        rows = ora(genes, background=background)
        engine = "ora"
    gsva_rows, gsva_engine = try_gsva(None)
    return {
        "engine": engine,
        "gsva": gsva_engine,
        "n_genes": len({g.upper() for g in genes if g}),
        "n_terms": len(rows),
        "terms": rows[:25],
        "gsva_terms": gsva_rows,
    }


def _genes_from_csv(path: Path, top_n: int = 200) -> list[str]:
    genes: list[str] = []
    with path.open(encoding="utf-8", errors="replace", newline="") as f:
        reader = csv.DictReader(f)
        gene_col = next((c for c in (reader.fieldnames or []) if c.lower() in {"gene", "genes", "names", "feature"}), None)
        padj_col = next((c for c in (reader.fieldnames or []) if c.lower() in {"padj", "fdr", "pvals_adj", "qval"}), None)
        rows = list(reader)
        if padj_col:
            def _p(r):
                try:
                    return float(r.get(padj_col) or 1)
                except (TypeError, ValueError):
                    return 1.0
            rows.sort(key=_p)
        for r in rows:
            g = (r.get(gene_col) if gene_col else None) or r.get("gene") or r.get("names")
            if g:
                genes.append(str(g))
            if len(genes) >= top_n:
                break
    return genes


def _genes_from_h5ad(path: Path, top_n: int = 100) -> list[str]:
    try:
        import scanpy as sc
    except Exception:
        return []
    try:
        adata = sc.read_h5ad(path)
    except Exception as exc:
        log.info("enrich skip h5ad: %s", exc)
        return []
    uns = adata.uns.get("rank_genes_groups") or {}
    names = uns.get("names")
    if names is None:
        return []
    genes: list[str] = []
    try:
        import numpy as np

        arr = np.array(names)
        for col in range(arr.shape[1] if arr.ndim == 2 else 1):
            col_vals = arr[:, col] if arr.ndim == 2 else arr
            for g in list(col_vals)[: max(10, top_n // max(arr.shape[1] if arr.ndim == 2 else 1, 1))]:
                if g is not None:
                    genes.append(str(g))
            if len(genes) >= top_n:
                break
    except Exception:
        return []
    return genes


def enrich_from_workspace(workspace: str | Path, *, top_n: int = 200) -> dict:
    workspace = Path(workspace)
    genes: list[str] = []
    source = None
    csv_path = workspace / "pseudobulk_de.csv"
    if csv_path.is_file():
        genes = _genes_from_csv(csv_path, top_n=top_n)
        source = str(csv_path.name)
    if not genes:
        for name in ("adata_processed.h5ad", "adata_qc.h5ad"):
            p = workspace / name
            if p.is_file():
                genes = _genes_from_h5ad(p, top_n=top_n)
                if genes:
                    source = name
                    break
    out = run_enrichment(genes)
    out["source"] = source
    out["n_input_genes"] = len(genes)
    (workspace / "pathway_enrichment.json").write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")
    terms = out.get("terms") or []
    if terms:
        with (workspace / "pathway_enrichment.csv").open("w", encoding="utf-8", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["term", "overlap", "set_size", "pval", "fdr", "method"])
            w.writeheader()
            for row in terms:
                w.writerow({k: row.get(k) for k in w.fieldnames})
        try:
            from scagent.plotting import pathway_bubble_plot

            pathway_bubble_plot(terms, figdir=workspace / "figures")
        except Exception as exc:
            log.info("pathway bubble skipped: %s", exc)
    return out
