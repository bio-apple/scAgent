"""Analysis steps: PCA, neighbors, Leiden, UMAP. Parameters from config.yaml unless overridden."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from scagent.config import analysis_params
from scagent.logutil import get_logger
from scagent.parallel import apply_scanpy_n_jobs, map_parallel

log = get_logger("analysis")


def pca(adata, *, n_pcs: int | None = None, cache_name: str | None = None):
    import scanpy as sc

    from scagent.inspect_data import detect_expression_layer

    if cache_name:
        from scagent.cache import load_h5ad

        hit = load_h5ad(cache_name)
        if hit is not None:
            return hit
    n = n_pcs or analysis_params()["n_pcs"]
    info = detect_expression_layer(adata)
    if info.get("layer") == "scaled":
        log.info("pca skip scale: X already scaled (%s)", info.get("reason"))
    else:
        if info.get("layer") == "counts":
            log.warning("pca: X looks like counts; scaling raw counts. Prefer log1p first.")
        sc.pp.scale(adata, max_value=analysis_params()["scale_max_value"])
    existing = adata.obsm["X_pca"] if "X_pca" in adata.obsm else None
    if existing is not None and getattr(existing, "shape", (0, 0))[1] >= n:
        log.info("pca skip: X_pca already present n_pcs=%s", existing.shape[1])
    else:
        sc.tl.pca(adata, n_comps=n, svd_solver="arpack", use_highly_variable=True)
        log.info("pca n_pcs=%s", n)
    if cache_name:
        from scagent.cache import save_h5ad

        save_h5ad(cache_name, adata)
    return adata


def neighbors(adata, *, n_neighbors: int | None = None, n_pcs: int | None = None, use_rep: str | None = None):
    import scanpy as sc

    p = analysis_params()
    n_neighbors = n_neighbors or p["n_neighbors"]
    n_pcs = n_pcs or p["n_pcs"]
    kwargs: dict = {"n_neighbors": n_neighbors}
    if use_rep:
        kwargs["use_rep"] = use_rep
    else:
        kwargs["n_pcs"] = min(n_pcs, adata.obsm["X_pca"].shape[1])
    sc.pp.neighbors(adata, **kwargs)
    log.info("neighbors %s", kwargs)
    return adata


def leiden(adata, *, resolution: float | None = None, key_added: str = "leiden"):
    import scanpy as sc

    p = analysis_params()
    res = resolution if resolution is not None else p.get("leiden_resolution") or 0.6
    sc.tl.leiden(adata, resolution=float(res), key_added=key_added)
    log.info("leiden resolution=%s n_clusters=%s", res, adata.obs[key_added].nunique())
    return adata


def umap(adata):
    import scanpy as sc

    sc.tl.umap(adata)
    return adata


def rank_genes(
    adata,
    groupby: str = "leiden",
    *,
    method: str | None = None,
    cross_validate: Any = None,
):
    """Exploratory cluster markers (cell-level). Not a sample-level condition DE.

    method: auto|wilcoxon|t-test|mast. MAST is optional R; missing package is skipped.
    cross_validate auto runs a second Scanpy test (wilcoxon ∩ t-test) for Reviewer overlap.
    """
    import json

    import scanpy as sc

    from scagent.deg_methods import (
        alt_scanpy_method,
        gene_overlap,
        resolve_cross_validate,
        resolve_marker_method,
    )

    apply_scanpy_n_jobs()
    cfg = _deg_cfg()
    primary = resolve_marker_method(method or cfg.get("marker_method") or "auto")
    do_cv = resolve_cross_validate(cfg.get("cross_validate"), explicit=_cv_explicit(cross_validate))
    use_raw = adata.raw is not None
    if not use_raw:
        log.warning("rank_genes: adata.raw missing; tests may run on scaled X")
        print("SCAGENT_WARN: rank_genes_groups without .raw; Wilcoxon/t-test may run on scaled X")
    scanpy_method = "wilcoxon" if primary == "mast" else primary
    sc.tl.rank_genes_groups(adata, groupby, method=scanpy_method, pts=True, use_raw=use_raw)
    methods = [scanpy_method]
    mast_status = None
    if primary == "mast":
        mast_status = _try_mast_markers(adata, groupby)
        if mast_status:
            methods.append("mast:" + str(mast_status))
    overlap = None
    if do_cv:
        alt = alt_scanpy_method(scanpy_method)
        try:
            sc.tl.rank_genes_groups(adata, groupby, method=alt, pts=True, use_raw=use_raw, key_added="rank_genes_groups_alt")
            methods.append(alt)
            a = _marker_gene_set(adata, "rank_genes_groups")
            b = _marker_gene_set(adata, "rank_genes_groups_alt")
            overlap = gene_overlap(a, b)
            Path("cluster_marker_overlap.json").write_text(json.dumps(overlap, indent=2), encoding="utf-8")
        except Exception as exc:
            log.info("marker cross-validate skipped: %s", exc)
    payload = {
        "method": scanpy_method,
        "requested": primary,
        "methods": methods,
        "cross_validate": bool(overlap),
        "n_overlap": None if not overlap else overlap.get("n_overlap"),
        "jaccard": None if not overlap else overlap.get("jaccard"),
        "mast": mast_status,
        "note": "exploratory cell-level markers; condition DE needs pseudobulk + FDR",
    }
    adata.uns["scagent_markers"] = payload
    print(
        "Exploratory cluster markers ("
        + ",".join(methods)
        + "); pvals_adj/FDR. Not a group-level result. Use pseudobulk + FDR for condition DE."
    )
    if overlap:
        print("marker_crossvalidate overlap=" + str(overlap.get("n_overlap")) + " jaccard=" + str(overlap.get("jaccard")))
    log.info("rank_genes_groups groupby=%s methods=%s use_raw=%s", groupby, methods, use_raw)
    return adata


def _cv_explicit(raw: Any) -> bool | None:
    if raw is None or raw == "auto":
        return None
    if isinstance(raw, bool):
        return raw
    s = str(raw).lower().strip()
    if s in {"off", "false", "0", "no", "none"}:
        return False
    if s in {"on", "true", "1", "yes", "always", "both"}:
        return True
    return None


def _marker_gene_set(adata, key: str, *, n_top: int = 25, fdr: float = 0.05) -> list[str]:
    uns = adata.uns.get(key) or {}
    names = uns.get("names")
    padj = uns.get("pvals_adj")
    if names is None:
        return []
    groups = names.dtype.names or ()
    genes: list[str] = []
    seen: set[str] = set()
    for g in groups:
        n_keep = 0
        for i, gene in enumerate(names[g]):
            gname = str(gene)
            if padj is not None:
                try:
                    if float(padj[g][i]) > fdr:
                        continue
                except (TypeError, ValueError, IndexError):
                    pass
            if gname not in seen:
                seen.add(gname)
                genes.append(gname)
            n_keep += 1
            if n_keep >= n_top:
                break
    return genes


def _try_mast_markers(adata, groupby: str) -> str | None:
    """Optional MAST on a HVG subset vs the most frequent cluster. Honest skip if R/MAST missing."""
    import shutil
    import subprocess
    import tempfile

    script = Path(__file__).resolve().parent / "r" / "mast.R"
    rscript = shutil.which("Rscript")
    if not script.is_file() or not rscript:
        return "skipped_no_rscript"
    if groupby not in adata.obs:
        return "skipped_no_groupby"
    n_obs = int(adata.n_obs)
    if n_obs > 4000:
        return "skipped_too_large"
    import numpy as np
    import pandas as pd

    expr = adata.raw.to_adata() if adata.raw is not None else adata
    genes = list(map(str, expr.var_names[: min(200, expr.n_vars)]))
    X = expr[:, genes].X
    if hasattr(X, "toarray"):
        X = X.toarray()
    X = np.asarray(X).T
    top = str(adata.obs[groupby].astype(str).value_counts().index[0])
    grp = (adata.obs[groupby].astype(str) == top).map({True: "in", False: "out"})
    with tempfile.TemporaryDirectory() as tmp:
        tmpp = Path(tmp)
        cpath = tmpp / "counts.csv"
        mpath = tmpp / "meta.csv"
        opath = tmpp / "out.csv"
        pd.DataFrame(X, index=genes, columns=list(map(str, adata.obs_names))).to_csv(cpath)
        pd.DataFrame({"cell": list(map(str, adata.obs_names)), "group": list(grp)}).to_csv(mpath, index=False)
        try:
            proc = subprocess.run(
                [rscript, str(script), str(cpath), str(mpath), str(opath)],
                capture_output=True,
                text=True,
                timeout=120,
                check=False,
            )
        except Exception as exc:
            log.info("MAST skipped: %s", exc)
            return "skipped"
        if opath.is_file():
            txt = opath.read_text(encoding="utf-8")[:200]
            if '"status":"skipped"' in txt.replace(" ", ""):
                return "skipped"
            if "gene" in txt:
                Path("mast_markers.csv").write_text(opath.read_text(encoding="utf-8"), encoding="utf-8")
                return "ok"
        log.info("MAST skipped: %s", (proc.stderr or proc.stdout or "")[-300:])
    return "skipped"


def score_markers_by_cluster(cluster_ids: list, score_fn, *, jobs: int | None = None) -> list:
    """Parallel per-cluster marker scoring. score_fn(cluster_id) -> result."""
    return map_parallel(score_fn, cluster_ids, jobs=jobs)


def needs_condition_de(query: str | None) -> bool:
    from agents.intent import parse_intent

    return bool(parse_intent(query).get("condition_comparison"))


_DEG_R = Path(__file__).resolve().parent / "r" / "deg.R"


def _deg_cfg() -> dict:
    try:
        from scagent.config import load_config

        return dict(load_config().get("deg") or {})
    except Exception:
        return {}


def _bh_fdr(pvals):
    import numpy as np

    p = np.asarray(pvals, dtype=float)
    n = max(len(p), 1)
    order = np.argsort(p)
    ranked = np.clip(p[order], 0, 1)
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(q, 0, 1)
    return out


def _aggregate_counts(adata, keys: list[str], layer: str = "counts"):
    import numpy as np
    import pandas as pd

    if layer not in adata.layers:
        if "counts_raw" in adata.layers:
            layer = "counts_raw"
            log.info("pseudobulk_de: using layers['counts_raw']")
        else:
            log.warning("pseudobulk_de: no counts layer; summing X (must be raw counts, not log1p)")
    if layer in adata.layers:
        X = adata.layers[layer]
    else:
        X = adata.X
    from scipy import sparse as sp

    if sp.issparse(X):
        X = X.tocsr()
    else:
        X = np.asarray(X)
    grp = adata.obs[keys].astype(str)
    grp["_row"] = np.arange(adata.n_obs)
    rows = []
    index = []
    for name, sub in grp.groupby(keys, observed=True):
        idx = sub["_row"].to_numpy()
        block = X[idx]
        summed = np.asarray(block.sum(axis=0)).ravel()
        rows.append(summed)
        index.append(name if isinstance(name, tuple) else (name,))
    mat = np.vstack(rows) if rows else np.zeros((0, adata.n_vars))
    cols = [str(g) for g in adata.var_names]
    idx = pd.MultiIndex.from_tuples(index, names=keys) if index else None
    return pd.DataFrame(mat, index=idx, columns=cols)


def pseudobulk_de(
    adata,
    *,
    sample_key: str,
    condition_key: str,
    groupby: str = "cell_type",
    layer: str = "counts",
    min_cells: int | None = None,
    min_replicates: int | None = None,
    engine: str | None = None,
    cross_validate: Any = None,
):
    """Sample × cell-type sum of raw counts, then edgeR/DESeq2 (rpy2) or t-test + BH-FDR."""
    import json

    import pandas as pd

    from scagent.deg_methods import alt_engine, gene_overlap, resolve_cross_validate, sig_genes

    cfg = _deg_cfg()
    min_cells = int(min_cells if min_cells is not None else cfg.get("min_cells") or 10)
    min_replicates = int(min_replicates if min_replicates is not None else cfg.get("min_replicates") or 2)
    engine = str(engine or cfg.get("engine") or "auto").lower().strip()
    if engine in {"edgeR", "edger_rpy2"}:
        engine = "edger"
    if engine in {"deseq2_rpy2"}:
        engine = "deseq2"
    if engine == "mast":
        log.warning("MAST is cell-level; condition DE still uses pseudobulk (edgeR/DESeq2/t-test)")
        engine = "auto"
    do_cv = resolve_cross_validate(cfg.get("cross_validate"), explicit=_cv_explicit(cross_validate))

    if sample_key not in adata.obs or condition_key not in adata.obs:
        log.warning("pseudobulk_de skipped: missing %s or %s", sample_key, condition_key)
        adata.uns["pseudobulk_de"] = {"ran": False, "reason": "missing_keys"}
        return adata
    if groupby not in adata.obs:
        groupby = "leiden" if "leiden" in adata.obs else None
    if not groupby:
        adata.uns["pseudobulk_de"] = {"ran": False, "reason": "no_cell_type"}
        return adata
    keys = [sample_key, condition_key, groupby]
    n_cells = adata.obs.groupby(keys, observed=True).size()
    keep = n_cells[n_cells >= min_cells]
    if keep.empty:
        adata.uns["pseudobulk_de"] = {"ran": False, "reason": "too_few_cells"}
        return adata
    pb = _aggregate_counts(adata, keys, layer=layer)
    pb = pb.loc[keep.index.intersection(pb.index)]
    results: list[dict] = []
    used = "ttest_bh"
    if engine != "ttest":
        r_out, used_r = _pseudobulk_r_backend(pb, sample_key, condition_key, groupby, engine=engine, min_replicates=min_replicates)
        if r_out:
            results = r_out
            used = used_r or "edger_rpy2"
        elif engine in {"edger", "deseq2"}:
            log.warning("requested DEG engine=%s unavailable; falling back to t-test+BH", engine)
    if not results:
        results = _pseudobulk_ttest(pb, sample_key, condition_key, groupby, min_replicates=min_replicates)
        used = "ttest_bh"
    df = pd.DataFrame(results)
    engines = [used]
    overlap = None
    if do_cv and len(df):
        alt = alt_engine(used if engine == "auto" else engine)
        alt_rows: list[dict] = []
        alt_used = None
        if alt != "ttest":
            alt_rows, alt_used = _pseudobulk_r_backend(
                pb, sample_key, condition_key, groupby, engine=alt, min_replicates=min_replicates
            )
        if not alt_rows and "ttest" not in used:
            alt_rows = _pseudobulk_ttest(pb, sample_key, condition_key, groupby, min_replicates=min_replicates)
            alt_used = "ttest_bh"
        if alt_rows and alt_used:
            engines.append(alt_used)
            df2 = pd.DataFrame(alt_rows)
            df2.to_csv("pseudobulk_de_alt.csv", index=False)
            overlap = gene_overlap(sig_genes(results), sig_genes(alt_rows))
            Path("pseudobulk_de_overlap.json").write_text(json.dumps(overlap, indent=2), encoding="utf-8")
    adata.uns["pseudobulk_de"] = {
        "ran": bool(len(df)),
        "engine": used,
        "requested_engine": engine,
        "engines": engines,
        "n_tests": int(len(df)),
        "min_replicates": min_replicates,
        "cross_validate": bool(overlap),
        "n_overlap": None if not overlap else overlap.get("n_overlap"),
        "jaccard": None if not overlap else overlap.get("jaccard"),
        "note": "sample-level pseudobulk + FDR; not cell-level Wilcoxon/MAST",
    }
    if len(df):
        df.to_csv("pseudobulk_de.csv", index=False)
        adata.uns["pseudobulk_de"]["n_sig"] = int((df.get("fdr", pd.Series(dtype=float)) < 0.05).sum())
    log.info("pseudobulk_de engine=%s n=%s cv=%s", used, len(df), bool(overlap))
    return adata


def _pb_celltype_blocks(pb, sample_key, condition_key, groupby, *, min_replicates: int):
    """Yield (cell_type, genes_x_samples matrix, sample ids, conditions, group_a, group_b)."""
    import numpy as np

    names = list(pb.index.names)
    si, ci, gi = names.index(sample_key), names.index(condition_key), names.index(groupby)
    grouped: dict = {}
    for idx, row in pb.iterrows():
        grouped.setdefault(idx[gi], []).append((str(idx[si]), str(idx[ci]), row.to_numpy()))
    genes = [str(g) for g in pb.columns]
    for ct, rows in grouped.items():
        conds = sorted({r[1] for r in rows})
        if len(conds) != 2:
            continue
        a = [r for r in rows if r[1] == conds[0]]
        b = [r for r in rows if r[1] == conds[1]]
        if len(a) < min_replicates or len(b) < min_replicates:
            log.warning("pseudobulk %s: need ≥%s replicates/condition, got %s vs %s", ct, min_replicates, len(a), len(b))
            continue
        ordered = a + b
        samples = [r[0] for r in ordered]
        conditions = [r[1] for r in ordered]
        mat = np.vstack([r[2] for r in ordered]).T
        yield str(ct), mat, samples, conditions, conds[0], conds[1], genes


def _read_deg_csv(path: Path, cell_type: str, group_a: str, group_b: str) -> list[dict]:
    import csv

    rows: list[dict] = []
    with path.open(encoding="utf-8") as f:
        for rec in csv.DictReader(f):
            gene = rec.get("gene")
            if not gene:
                continue
            pval = rec.get("pval")
            fdr = rec.get("fdr")
            lfc = rec.get("logFC")
            if any(v in (None, "", "NA", "NaN") for v in (pval, fdr, lfc)):
                continue
            try:
                pval_f = float(pval)
                fdr_f = float(fdr)
                lfc_f = float(lfc)
            except ValueError:
                continue
            rows.append(
                {
                    "cell_type": cell_type,
                    "gene": gene,
                    "logFC": lfc_f,
                    "pval": pval_f,
                    "fdr": fdr_f,
                    "group_a": group_a,
                    "group_b": group_b,
                }
            )
    return rows


def _write_deg_inputs(tmp, mat, samples, conditions, genes):
    import pandas as pd

    cpath = tmp / "counts.csv"
    mpath = tmp / "meta.csv"
    cdf = pd.DataFrame(mat, index=genes, columns=samples)
    cdf.index.name = "gene"
    cdf.to_csv(cpath)
    pd.DataFrame({"sample": samples, "condition": conditions}).to_csv(mpath, index=False)
    return cpath, mpath


def _deg_via_rpy2(cpath: Path, mpath: Path, opath: Path, engine: str) -> str | None:
    try:
        import rpy2.robjects as ro
    except Exception:
        return None
    try:
        path = str(_DEG_R).replace("\\", "/")
        ro.r(f'source("{path}")')
        used = ro.r["run_deg"](str(cpath), str(mpath), str(opath), engine)
        if not opath.is_file():
            return None
        return str(used[0]) if used is not None else engine
    except Exception as exc:
        log.info("DEG rpy2 failed (%s): %s", engine, exc)
        return None


def _deg_via_rscript(cpath: Path, mpath: Path, opath: Path, engine: str) -> str | None:
    import shutil
    import subprocess

    if not shutil.which("Rscript") or not _DEG_R.is_file():
        return None
    try:
        proc = subprocess.run(
            ["Rscript", str(_DEG_R), str(cpath), str(mpath), str(opath), engine],
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        if proc.returncode != 0 or not opath.is_file():
            log.info("DEG Rscript skipped: %s", (proc.stderr or proc.stdout or "")[-400:])
            return None
        used = (proc.stdout or "").strip().splitlines()
        return used[-1] if used else engine
    except Exception as exc:
        log.info("DEG Rscript failed: %s", exc)
        return None


def _pseudobulk_r_backend(pb, sample_key, condition_key, groupby, *, engine: str, min_replicates: int):
    """Try rpy2 (edgeR then DESeq2), then Rscript. Returns (rows, engine_name) or ([], None)."""
    import tempfile

    r_engine = "auto" if engine == "auto" else engine
    rows: list[dict] = []
    used_name = None
    for ct, mat, samples, conditions, ga, gb, genes in _pb_celltype_blocks(
        pb, sample_key, condition_key, groupby, min_replicates=min_replicates
    ):
        tmpdir = tempfile.TemporaryDirectory()
        try:
            tmp = Path(tmpdir.name)
            cpath, mpath = _write_deg_inputs(tmp, mat, samples, conditions, genes)
            opath = tmp / "out.csv"
            used = _deg_via_rpy2(cpath, mpath, opath, r_engine)
            tag = "rpy2"
            if used is None:
                used = _deg_via_rscript(cpath, mpath, opath, r_engine)
                tag = "rscript"
            if used is None or not opath.is_file():
                continue
            part = _read_deg_csv(opath, ct, ga, gb)
            if not part:
                continue
            rows.extend(part)
            used_name = f"{used}_{tag}"
        finally:
            tmpdir.cleanup()
    if not rows:
        return [], None
    return rows, used_name


def _pseudobulk_edger_rpy2(pb, sample_key, condition_key, groupby):
    """Backward-compatible alias: edgeR via rpy2/Rscript."""
    rows, _used = _pseudobulk_r_backend(
        pb, sample_key, condition_key, groupby, engine="edger", min_replicates=2
    )
    return rows or None


def _pseudobulk_ttest(pb, sample_key, condition_key, groupby, *, min_replicates: int = 2) -> list[dict]:
    import numpy as np
    from scipy.stats import ttest_ind

    results: list[dict] = []
    names = list(pb.index.names)
    si, ci, gi = names.index(sample_key), names.index(condition_key), names.index(groupby)
    grouped = {}
    for idx, row in pb.iterrows():
        key = idx[gi]
        grouped.setdefault(key, []).append((idx[si], idx[ci], row.to_numpy()))
    genes = list(pb.columns)
    for ct, rows in grouped.items():
        conds = sorted({r[1] for r in rows})
        if len(conds) != 2:
            continue
        a = [r[2] for r in rows if r[1] == conds[0]]
        b = [r[2] for r in rows if r[1] == conds[1]]
        if len(a) < min_replicates or len(b) < min_replicates:
            log.warning("pseudobulk %s: need ≥%s replicates/condition, got %s vs %s", ct, min_replicates, len(a), len(b))
            continue
        A, B = np.vstack(a), np.vstack(b)
        # CPM
        A = np.log1p(A / np.maximum(A.sum(axis=1, keepdims=True), 1) * 1e6)
        B = np.log1p(B / np.maximum(B.sum(axis=1, keepdims=True), 1) * 1e6)
        pvals = []
        lfc = []
        for j in range(A.shape[1]):
            t = ttest_ind(A[:, j], B[:, j], equal_var=False, nan_policy="omit")
            pvals.append(float(t.pvalue) if t.pvalue == t.pvalue else 1.0)
            lfc.append(float(B[:, j].mean() - A[:, j].mean()))
        fdr = _bh_fdr(np.array(pvals))
        for j, gene in enumerate(genes):
            if fdr[j] > 0.2 and abs(lfc[j]) < 0.25:
                continue
            results.append(
                {
                    "cell_type": ct,
                    "gene": gene,
                    "logFC": lfc[j],
                    "pval": pvals[j],
                    "fdr": float(fdr[j]),
                    "group_a": conds[0],
                    "group_b": conds[1],
                }
            )
    return results


def _embed_matrix(adata, embed_key: str | None = None):
    import numpy as np
    if embed_key and embed_key in adata.obsm:
        return np.asarray(adata.obsm[embed_key]), embed_key
    for key in ("X_pca_harmony", "X_scVI", "X_scanorama", "X_pca"):
        if key in adata.obsm:
            return np.asarray(adata.obsm[key]), key
    return None, None


def knn_ilisi(X, batch, *, k: int = 30) -> float:
    """Scaled iLISI in [0, 1]: neighborhood batch entropy / log(n_batches)."""
    import numpy as np
    from sklearn.neighbors import NearestNeighbors

    batch = np.asarray(batch).astype(str)
    n_b = len(set(batch))
    if n_b < 2 or X is None or len(X) < 4:
        return 1.0
    k = int(min(max(5, k), len(X) - 1))
    nn = NearestNeighbors(n_neighbors=k + 1).fit(X)
    idx = nn.kneighbors(return_distance=False)[:, 1:]
    scores = []
    for row in idx:
        labs, counts = np.unique(batch[row], return_counts=True)
        p = counts / counts.sum()
        ent = float(-(p * np.log(p + 1e-12)).sum())
        scores.append(ent / np.log(n_b))
    return float(np.mean(scores))


def pca_batch_r2(X, batch) -> float:
    """Mean R² of PCs ~ one-hot batch (lower = less batch in embedding)."""
    import numpy as np

    batch = np.asarray(batch).astype(str)
    if X is None or len(set(batch)) < 2:
        return 0.0
    dummies = np.eye(len(set(batch)))[pd_factorize(batch)]
    dummies = dummies[:, :-1]
    if dummies.size == 0:
        return 0.0
    r2s = []
    for j in range(min(X.shape[1], 20)):
        y = X[:, j]
        y = y - y.mean()
        if float(np.dot(y, y)) == 0:
            continue
        beta, *_ = np.linalg.lstsq(dummies, y, rcond=None)
        pred = dummies @ beta
        r2s.append(1.0 - float(np.dot(y - pred, y - pred)) / float(np.dot(y, y)))
    return float(np.mean(r2s)) if r2s else 0.0


def pd_factorize(batch):
    import numpy as np

    labs, inv = np.unique(np.asarray(batch).astype(str), return_inverse=True)
    _ = labs
    return inv


def integration_quality(adata, batch_key: str, *, embed_key: str | None = None) -> dict:
    """scIB iLISI/kBET if installed; else kNN-iLISI + PCA batch R²."""
    from scagent.config import load_config

    cfg = (load_config().get("integration") or {})
    ilisi_min = float(cfg.get("ilisi_min") or 0.8)
    kbet_min = float(cfg.get("kbet_min") or 0.5)
    r2_max = float(cfg.get("pca_batch_r2_max") or 0.5)
    out = {
        "ilisi": None,
        "kbet": None,
        "pca_batch_r2": None,
        "embed_key": None,
        "passed": True,
        "engine": "none",
    }
    if batch_key not in adata.obs or adata.obs[batch_key].nunique() < 2:
        out["passed"] = True
        out["engine"] = "single_batch"
        return out
    X, used = _embed_matrix(adata, embed_key)
    out["embed_key"] = used
    batch = adata.obs[batch_key]
    try:
        import scib.metrics as sm

        if used:
            adata.obsm.setdefault("_scagent_embed", X)
            try:
                out["ilisi"] = float(sm.ilisi_graph(adata, batch_key=batch_key, type_="embed", use_rep=used))
                out["engine"] = "scib"
            except Exception:
                pass
            try:
                kbet = sm.kBET(adata, batch_key=batch_key, type_="embed", embed=used)
                out["kbet"] = float(kbet) if kbet is not None else None
                out["engine"] = "scib"
            except Exception:
                pass
    except Exception:
        pass
    if X is not None:
        if out["ilisi"] is None:
            out["ilisi"] = knn_ilisi(X, batch)
            out["engine"] = out["engine"] if out["engine"] != "none" else "knn_ilisi"
        out["pca_batch_r2"] = pca_batch_r2(X, batch)
    fails = []
    if out["ilisi"] is not None and out["ilisi"] < ilisi_min:
        fails.append(f"iLISI {out['ilisi']:.3f} < {ilisi_min}")
    if out["kbet"] is not None and out["kbet"] < kbet_min:
        fails.append(f"kBET {out['kbet']:.3f} < {kbet_min}")
    if out["ilisi"] is None and out["kbet"] is None and out["pca_batch_r2"] is not None and out["pca_batch_r2"] > r2_max:
        fails.append(f"PCA batch R² {out['pca_batch_r2']:.3f} > {r2_max}")
    out["passed"] = not fails
    out["issues"] = fails
    try:
        from scagent.plotting import integration_diagnostics

        out["plots"] = integration_diagnostics(adata, batch_key)
    except Exception as exc:
        log.warning("integration diagnostics skipped: %s", exc)
        out["plots"] = []
    log.info("integration_quality %s", {k: out[k] for k in out if k != "plots"})
    return out
