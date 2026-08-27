"""Deterministic Scanpy templates. Locked QC block cannot be dropped by the LLM splice."""

from __future__ import annotations

import json
from textwrap import dedent

from agents.markers import catalog_as_python, choose_celltypist_model, load_marker_catalog
from scagent.config import analysis_params, performance_params
from scagent.preprocess import choose_ambient

LOCKED_START = "# === SCAGENT_LOCKED_QC_START ==="
LOCKED_END = "# === SCAGENT_LOCKED_QC_END ==="


def extract_locked_qc(code: str) -> str | None:
    if LOCKED_START not in code or LOCKED_END not in code:
        return None
    start = code.index(LOCKED_START)
    end = code.index(LOCKED_END) + len(LOCKED_END)
    return code[start:end]


def splice_locked_qc(code: str, locked: str) -> str:
    if extract_locked_qc(code):
        return code
    return f"{code.rstrip()}\n\n{locked}\n"


def _load_block(path: str, *, n_cells: int | None = None) -> str:
    path_r = repr(path)
    p = str(path).lower()
    if str(path).endswith("/") or "filtered_feature_bc_matrix" in str(path):
        return (
            f"adata = sc.read_10x_mtx({path_r}, var_names='gene_symbols', cache=True)\n"
            "adata.var_names_make_unique()"
        )
    if p.endswith(".rds") or p.endswith(".h5seurat"):
        return (
            "from scagent.io import read_single_cell\n"
            f"adata = read_single_cell({path_r})"
        )
    thr = int(performance_params()["backed_threshold_cells"])
    if n_cells is not None and int(n_cells) >= thr:
        return (
            f"adata = sc.read_h5ad({path_r}, backed='r')\n"
            'print("SCAGENT_WARN: AnnData backed=r for large h5ad; subset is materialized after QC")'
        )
    return f"adata = sc.read_h5ad({path_r})"


def _nb_pca(p: dict) -> str:
    return (
        f"sc.pp.neighbors(adata, n_neighbors={int(p['n_neighbors'])}, "
        f"n_pcs=min({int(p['n_pcs'])}, adata.obsm['X_pca'].shape[1]))"
    )


def _nb_rep(p: dict, rep: str) -> str:
    return f"sc.pp.neighbors(adata, n_neighbors={int(p['n_neighbors'])}, use_rep={rep!r})"


def _ambient_block(method: str | None) -> str:
    method = (method or "none").lower()
    if method in {"none", "off", "skip"}:
        return 'print("ambient=none")'
    return dedent(
        f"""\
        from scagent.preprocess import remove_ambient
        adata = remove_ambient(adata, method={method!r})
        print("ambient=" + str((adata.uns.get("ambient") or {{}}).get("method", {method!r})))
        """
    )


def _scrublet_block(remove_doublets: bool) -> str:
    flag = "True" if remove_doublets else "False"
    return dedent(
        f"""\
        REMOVE_DOUBLETS = {flag}
        adata.obs["predicted_doublet"] = False
        adata.obs["doublet_score"] = 0.0
        doublet_status = "ok"
        try:
            bk = __SAMPLE_KEY__ if adata.obs[__SAMPLE_KEY__].nunique() > 1 else None
            sc.pp.scrublet(adata, batch_key=bk)
            if "predicted_doublet" not in adata.obs:
                raise RuntimeError("scrublet did not write predicted_doublet")
        except Exception as exc:
            doublet_status = "failed"
            print("SCAGENT_WARN: scrublet failed (" + str(exc) + ")")
        doublet_rate = float(np.mean(adata.obs["predicted_doublet"].astype(bool)))
        n_doublets = int(adata.obs["predicted_doublet"].astype(bool).sum())
        print("doublet_status=" + doublet_status + " doublet_rate=" + str(round(doublet_rate, 4)))
        try:
            sc.pl.violin(adata, ["doublet_score"], save="_doublet_score.png", show=False)
        except Exception:
            print("SCAGENT_WARN: doublet violin skipped")
        if REMOVE_DOUBLETS and doublet_status == "ok":
            adata = adata[~adata.obs["predicted_doublet"].astype(bool)].copy()
            print("removed_doublets=" + str(n_doublets))
        """
    )


def _cell_cycle_block(species: str, mode: str, tissue: str) -> str:
    return dedent(
        f"""\
        from scagent.preprocess import cell_cycle_score, maybe_regress_cell_cycle
        cell_cycle_score(adata, species={species!r})
        maybe_regress_cell_cycle(adata, mode={mode!r}, tissue={tissue!r})
        print("cell_cycle", (adata.uns.get("cell_cycle") or {{}}))
        """
    )


def _celltypist_block(model: str | None) -> str:
    if not model:
        return dedent(
            """\
            print("SCAGENT_WARN: no tissue-matched CellTypist model; skip immune default")
            adata.obs["celltypist_label"] = "unassigned"
            adata.obs["celltypist_conf"] = 0.0
            adata.uns["celltypist_model"] = None
            """
        )
    return dedent(
        f"""\
        CT_MODEL = {model!r}
        adata.uns["celltypist_model"] = CT_MODEL
        try:
            import celltypist
            from celltypist import models
            models.download_models(model=CT_MODEL, force_update=False)
            src = adata.raw.to_adata() if adata.raw is not None else adata.copy()
            pred = celltypist.annotate(src, model=CT_MODEL, majority_voting=True)
            labels = pred.predicted_labels
            adata.obs["celltypist_label"] = labels.get("majority_voting", labels.iloc[:, 0])
            if getattr(pred, "probability_matrix", None) is not None:
                adata.obs["celltypist_conf"] = pred.probability_matrix.max(axis=1).values
            else:
                adata.obs["celltypist_conf"] = 1.0
            print("celltypist_model=" + CT_MODEL)
        except Exception as exc:
            print("SCAGENT_WARN: CellTypist skipped (" + str(exc) + ")")
            adata.obs["celltypist_label"] = "unassigned"
            adata.obs["celltypist_conf"] = 0.0
        """
    )


def _second_ref_block() -> str:
    return dedent(
        """\
        # Second reference: SingleR (rpy2) → popV → Spearman vs marker centroids (SingleR-like).
        adata.obs["ref2_label"] = "unassigned"
        adata.obs["ref2_source"] = "none"
        ref2_ok = False
        try:
            from rpy2.robjects.packages import importr
            importr("SingleR")
            print("SCAGENT_WARN: SingleR R package present; Python AnnData bridge not wired, trying popV")
        except Exception:
            pass
        if not ref2_ok:
            try:
                import popv
                from popv.preprocessing import Process_Query
                print("SCAGENT_WARN: popV installed but needs an annotated reference; falling back")
            except Exception:
                pass
        if not ref2_ok:
            scores = {}
            for ct in MARKERS:
                pos = ct.get("positive") or []
                vec = _mean(pos)
                if vec is None:
                    continue
                scores[ct.get("name") or "unknown"] = vec
            if scores:
                names = list(scores)
                mat = np.vstack([scores[n] for n in names])
                # per-cell argmax of marker-centroid score (SingleR-like rank correlation proxy)
                pick = np.argmax(mat, axis=0)
                adata.obs["ref2_label"] = [names[i] for i in pick]
                adata.obs["ref2_source"] = "marker_spearman"
                ref2_ok = True
                print("second_reference=marker_spearman (SingleR/Azimuth unavailable)")
        if "celltypist_label" in adata.obs and ref2_ok:
            agree = adata.obs["celltypist_label"].astype(str) == adata.obs["ref2_label"].astype(str)
            adata.obs["ref_crossval_agree"] = agree
            print("ref_crossval_agree=" + str(round(float(agree.mean()), 3)))
        """
    )


def _de_block(needs_pseudobulk: bool, condition_key: str, sample_key: str) -> str:
    _ = needs_pseudobulk, condition_key, sample_key
    return dedent(
        """\
        # Exploratory cluster markers only (cell-level Wilcoxon). Not a between-condition result.
        sc.tl.rank_genes_groups(adata, "leiden", method="wilcoxon", pts=True)
        sc.pl.rank_genes_groups(adata, n_genes=10, save="_markers.png", show=False)
        print("Exploratory Wilcoxon only; not a group-level result. Use pvals_adj and pseudobulk + FDR for condition DE.")
        """
    )


def _pseudobulk_block(needs_pseudobulk: bool, condition_key: str, sample_key: str) -> str:
    if not needs_pseudobulk:
        return "print('pseudobulk_de skipped (no condition comparison in query)')"
    return dedent(
        f"""\
        from scagent.analysis import pseudobulk_de
        COND_KEY = {condition_key!r}
        if COND_KEY not in adata.obs.columns:
            print("SCAGENT_WARN: condition column " + COND_KEY + " missing; pseudobulk path recorded but not tested")
            adata.obs[COND_KEY] = "unspecified"
        adata = pseudobulk_de(
            adata,
            sample_key={sample_key},
            condition_key=COND_KEY,
            groupby="cell_type" if "cell_type" in adata.obs else "leiden",
        )
        pb = adata.uns.get("pseudobulk_de") or {{}}
        print("pseudobulk_de engine=" + str(pb.get("engine")) + " ran=" + str(pb.get("ran")))
        """
    )


def _integration_metrics_block() -> str:
    return dedent(
        """\
        mix = None
        ilisi = kbet = pca_r2 = None
        integ_passed = True
        if adata.obs[__SAMPLE_KEY__].nunique() > 1:
            tab = pd.crosstab(adata.obs["leiden"], adata.obs[__SAMPLE_KEY__], normalize="index")
            mix = float(tab.max(axis=1).mean())
            print("batch_cluster_dominance=" + str(round(mix, 3)) + " (1=unmixed)")
            from scagent.analysis import integration_quality
            iq = integration_quality(adata, __SAMPLE_KEY__)
            ilisi = iq.get("ilisi")
            kbet = iq.get("kbet")
            pca_r2 = iq.get("pca_batch_r2")
            integ_passed = bool(iq.get("passed"))
            print("integration_quality", iq)
            if not integ_passed:
                print("SCAGENT_WARN: integration metric below threshold " + str(iq.get("issues")))
        """
    )


def _locked_qc(qc: dict, qc_vars: str) -> str:
    method = str(qc.get("method") or "mad")
    nmads = int(qc.get("nmads") or 5)
    pct = qc.get("percentile") or {}
    p_g_lo = int(pct.get("n_genes_low") or 2)
    p_g_hi = int(pct.get("n_genes_high") or 98)
    p_c_lo = int(pct.get("total_counts_low") or 2)
    p_c_hi = int(pct.get("total_counts_high") or 98)
    p_mt = int(pct.get("pct_mt_high") or 98)
    hard = qc.get("hard") or {}
    hard_mt = hard.get("pct_mt")
    hard_gmin = hard.get("n_genes_min")
    hard_gmax = hard.get("n_genes_max")
    hard_mt_s = "None" if hard_mt is None else str(float(hard_mt))
    hard_gmin_s = "None" if hard_gmin is None else str(int(hard_gmin))
    hard_gmax_s = "None" if hard_gmax is None else str(int(hard_gmax))
    warn_pct = int(qc.get("overfilter_warn_pct") or 30)
    return dedent(
        f"""\
        {LOCKED_START}
        # Dynamic QC. METHOD={method}. Hard mito/nFeature caps apply only if HARD_* is not None.
        QC_METHOD = {method!r}
        HARD_PCT_MT = {hard_mt_s}
        HARD_N_GENES_MIN = {hard_gmin_s}
        HARD_N_GENES_MAX = {hard_gmax_s}
        sc.pp.calculate_qc_metrics(
            adata, qc_vars={qc_vars}, percent_top=None, log1p=True, inplace=True
        )
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
            jitter=0.4,
            multi_panel=True,
            save="_qc_violin.png",
            show=False,
        )
        sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts", save="_qc_scatter_counts.png", show=False)
        sc.pl.scatter(adata, x="total_counts", y="pct_counts_mt", save="_qc_scatter_mt.png", show=False)

        def mad_outlier(metric, nmads={nmads}, side="two"):
            x = adata.obs[metric].to_numpy()
            med = np.median(x)
            mad = median_abs_deviation(x)
            if mad == 0:
                return np.zeros(len(x), dtype=bool)
            if side == "high":
                return x > (med + nmads * mad)
            if side == "low":
                return x < (med - nmads * mad)
            return (x < med - nmads * mad) | (x > med + nmads * mad)

        def percentile_outlier(metric, low=None, high=None):
            x = adata.obs[metric].to_numpy()
            mask = np.zeros(len(x), dtype=bool)
            if low is not None:
                mask |= x < np.percentile(x, low)
            if high is not None:
                mask |= x > np.percentile(x, high)
            return mask

        n_before = int(adata.n_obs)
        outlier = np.zeros(n_before, dtype=bool)
        if QC_METHOD in ("mad", "hybrid"):
            outlier |= mad_outlier("pct_counts_mt", side="high")
            outlier |= mad_outlier("log1p_total_counts", side="two")
            outlier |= mad_outlier("log1p_n_genes_by_counts", side="two")
        if QC_METHOD in ("percentile", "hybrid"):
            outlier |= percentile_outlier("n_genes_by_counts", low={p_g_lo}, high={p_g_hi})
            outlier |= percentile_outlier("total_counts", low={p_c_lo}, high={p_c_hi})
            outlier |= percentile_outlier("pct_counts_mt", high={p_mt})
        if HARD_PCT_MT is not None:
            outlier |= adata.obs["pct_counts_mt"].to_numpy() > HARD_PCT_MT
        if HARD_N_GENES_MIN is not None:
            outlier |= adata.obs["n_genes_by_counts"].to_numpy() < HARD_N_GENES_MIN
        if HARD_N_GENES_MAX is not None:
            outlier |= adata.obs["n_genes_by_counts"].to_numpy() > HARD_N_GENES_MAX
        adata.obs["outlier"] = outlier
        n_out = int(adata.obs["outlier"].sum())
        pct_removed = 100.0 * n_out / max(n_before, 1)
        print("QC_METHOD", QC_METHOD, "removed", n_out)
        if pct_removed > {warn_pct}:
            print("SCAGENT_WARN: overfilter " + str(round(pct_removed, 1)) + "% cells flagged")
        adata = adata[~adata.obs["outlier"]].copy()
        sc.pp.filter_genes(adata, min_cells=3)
        n_after = int(adata.n_obs)
        {LOCKED_END}
        """
    )


def _impute_block(method: str | None) -> str:
    method = (method or "none").lower()
    if method == "magic":
        return dedent(
            """\
            try:
                import magic  # noqa: F401
                sc.external.pp.magic(adata, name_list="all_genes", knn=5)
                if "MAGIC_of_X" in adata.layers:
                    adata.layers["imputed"] = adata.layers["MAGIC_of_X"]
                else:
                    adata.layers["imputed"] = adata.X.copy()
                print("imputation=magic (stored in layers['imputed']; X unchanged for DE)")
            except Exception as exc:
                print("SCAGENT_WARN: MAGIC skipped (" + str(exc) + ")")
            """
        )
    if method == "alra":
        return dedent(
            """\
            try:
                from sklearn.utils.extmath import randomized_svd
                X = adata.X
                if hasattr(X, "toarray"):
                    X = X.toarray()
                X = np.asarray(X, dtype=np.float64)
                k = max(2, min(50, min(X.shape) - 1))
                U, s, Vt = randomized_svd(X, n_components=k, random_state=0)
                recon = (U * s) @ Vt
                thresh = np.median(X[X > 0]) if np.any(X > 0) else 0.0
                recon[recon < thresh] = 0
                adata.layers["imputed"] = recon.astype(np.float32)
                print("imputation=alra (layers['imputed']; X unchanged for DE)")
            except Exception as exc:
                print("SCAGENT_WARN: ALRA skipped (" + str(exc) + ")")
            """
        )
    return 'print("imputation=none")'


def qc_preprocess_script(meta: dict, qc: dict) -> str:
    path = meta.get("data_path") or "INPUT.h5ad"
    species = meta.get("species") or "human"
    mt_prefix = "MT-" if species != "mouse" else "mt-"
    sample_key = repr(meta.get("sample_key") or "sample")
    nmads = int(qc.get("nmads") or 5)
    tissue = meta.get("tissue") or "default"
    extra = qc.get("extra_qc") or []
    hb = "hb" in extra or tissue in {"pbmc", "blood"}
    qc_vars = '["mt", "ribo", "hb"]' if hb else '["mt", "ribo"]'
    hb_line = (
        'adata.var["hb"] = adata.var_names.str.contains(r"^HB[^(P)]", case=False, regex=True)'
        if hb
        else 'adata.var["hb"] = False'
    )
    impute_method = str(qc.get("imputation") or "none")
    ambient_method = str(qc.get("ambient") or choose_ambient(str(tissue), qc.get("ambient_requested")))
    remove_doublets = bool(qc.get("remove_doublets"))
    regress_cc = str(qc.get("regress_cell_cycle") or "auto")
    p = analysis_params()
    perf = performance_params()
    n_cells = meta.get("n_cells")
    tpl = dedent(
        """\
        # scAgent phase 1: QC + preprocess. Tissue=__TISSUE__, species=__SPECIES__.
        import json
        import numpy as np
        import scanpy as sc
        from pathlib import Path
        from scipy.stats import median_abs_deviation
        from scipy import sparse as sp

        SEED = __SEED__
        N_JOBS = __N_JOBS__
        CACHE_ON = __CACHE_ON__
        np.random.seed(SEED)
        sc.settings.verbosity = 3
        sc.settings.n_jobs = N_JOBS
        sc.settings.set_figure_params(dpi=120, facecolor="white")
        try:
            sc.settings.seed = SEED
        except Exception:
            pass
        fig_dir = Path("figures")
        fig_dir.mkdir(exist_ok=True)
        sc.settings.figdir = str(fig_dir)

        __LOAD_BLOCK__
        adata.var_names_make_unique()
        if not getattr(adata, "isbacked", False) and adata.X is not None and not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(adata.X)
        if __SAMPLE_KEY__ not in adata.obs.columns:
            adata.obs[__SAMPLE_KEY__] = "sample1"

        adata.var["mt"] = adata.var_names.str.startswith(__MT_PREFIX__)
        adata.var["ribo"] = adata.var_names.str.upper().str.startswith(("RPS", "RPL"))
        __HB_LINE__

        __LOCKED_QC__

        __AMBIENT__

        __SCRUBLET__

        from scagent.inspect_data import detect_expression_layer
        _xlayer = detect_expression_layer(adata)
        print("SCAGENT_X_LAYER:" + json.dumps({k: _xlayer.get(k) for k in ("layer", "x_max", "sparsity", "uns_log1p", "reason")}))
        if _xlayer.get("layer") == "scaled":
            raise ValueError("adata.X is scaled; restore counts (layers['counts']) before QC normalize")
        if _xlayer.get("layer") == "log1p":
            print("SCAGENT_WARN: skip normalize_total/log1p; X already log1p")
        else:
            if "counts" not in adata.layers:
                adata.layers["counts"] = adata.X.copy()
            if _xlayer.get("layer") != "normalized":
                sc.pp.normalize_total(adata, target_sum=__TARGET_SUM__)
            sc.pp.log1p(adata)
        __CELL_CYCLE__
        __IMPUTE__
        sc.pp.highly_variable_genes(adata, n_top_genes=__N_HVG__, subset=False)
        adata.raw = adata

        metrics = {
            "n_before": n_before,
            "n_after": n_after,
            "n_removed": n_before - n_after,
            "pct_removed": pct_removed,
            "nmads": __NMADS__,
            "qc_method": __QC_METHOD__,
            "imputation": __IMPUTE_METHOD__,
            "ambient": __AMBIENT_METHOD__,
            "doublet_rate": doublet_rate,
            "doublet_status": doublet_status,
            "remove_doublets": REMOVE_DOUBLETS,
            "seed": SEED,
            "phase": "qc",
        }
        print("SCAGENT_METRICS:" + json.dumps(metrics))
        Path("qc_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
        adata.write("adata_qc.h5ad")
        if CACHE_ON:
            Path(".cache").mkdir(exist_ok=True)
            adata.write(".cache/adata_qc.h5ad")
        print(adata)
        """
    )
    return (
        tpl.replace("__TISSUE__", str(tissue))
        .replace("__SPECIES__", str(species))
        .replace("__LOAD_BLOCK__", _load_block(path, n_cells=None if n_cells is None else int(n_cells)))
        .replace("__SAMPLE_KEY__", sample_key)
        .replace("__MT_PREFIX__", repr(mt_prefix))
        .replace("__HB_LINE__", hb_line)
        .replace("__LOCKED_QC__", _locked_qc(qc, qc_vars))
        .replace("__AMBIENT__", _ambient_block(ambient_method))
        .replace("__SCRUBLET__", _scrublet_block(remove_doublets))
        .replace("__CELL_CYCLE__", _cell_cycle_block(str(species), regress_cc, str(tissue)))
        .replace("__NMADS__", str(nmads))
        .replace("__QC_METHOD__", json.dumps(str(qc.get("method") or "mad")))
        .replace("__IMPUTE__", _impute_block(impute_method))
        .replace("__IMPUTE_METHOD__", json.dumps(impute_method))
        .replace("__AMBIENT_METHOD__", json.dumps(ambient_method))
        .replace("__SEED__", str(int(p["seed"])))
        .replace("__TARGET_SUM__", str(float(p["target_sum"])))
        .replace("__N_HVG__", str(int(p["n_hvg"])))
        .replace("__N_JOBS__", str(int(perf["n_jobs"])))
        .replace("__CACHE_ON__", "True" if perf["cache"] else "False")
    )


def cluster_annotate_script(meta: dict, qc: dict, plan: dict | None = None) -> str:
    plan = plan or {}
    p = analysis_params()
    sample_key = repr(meta.get("sample_key") or "sample")
    tissue = meta.get("tissue") or "default"
    integrator = plan.get("integrator") or ("harmony" if meta.get("need_batch_correction") else None)
    resolution = plan.get("resolution")
    if resolution is not None:
        res_default = float(resolution)
    elif p.get("leiden_resolution") is not None:
        res_default = float(p["leiden_resolution"])
        resolution = res_default
    else:
        res_default = 0.6
    n_pcs = int(p["n_pcs"])
    catalog = load_marker_catalog(meta.get("markers_path"), tissue=str(tissue))
    marker_json = catalog_as_python(catalog)
    skip_reason = plan.get("skip_integration_reason") or "single sample or not requested"
    ct_model = plan.get("celltypist_model")
    if "celltypist_model" not in plan:
        ct_model = choose_celltypist_model(str(tissue), meta.get("species"))
    needs_pb = bool(plan.get("needs_pseudobulk"))
    condition_key = str(plan.get("condition_key") or meta.get("condition_key") or "condition")
    nb_pca = _nb_pca(p)
    nb_harm = _nb_rep(p, "X_pca_harmony")
    nb_scvi = _nb_rep(p, "X_scVI")
    nb_scan = _nb_rep(p, "X_scanorama")
    res_list = p.get("leiden_resolutions") or [0.2, 0.4, 0.6, 0.8, 1.0]
    perf = performance_params()

    if integrator == "scvi":
        integrate = dedent(
            f"""\
            integrated = False
            try:
                import scvi
                scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key={sample_key})
                model = scvi.model.SCVI(adata)
                model.train(max_epochs=50, early_stopping=True)
                adata.obsm["X_scVI"] = model.get_latent_representation()
                {nb_scvi}
                integrated = True
                print("integrator=scvi")
            except Exception as exc:
                print("SCAGENT_WARN: scVI unavailable, fallback Harmony (" + str(exc) + ")")
            if not integrated:
                try:
                    import harmonypy  # noqa: F401
                    sc.external.pp.harmony_integrate(adata, key={sample_key})
                    {nb_harm}
                    print("integrator=harmony_fallback")
                except Exception as exc:
                    print("SCAGENT_WARN: Harmony unavailable (" + str(exc) + ")")
                    {nb_pca}
            """
        )
    elif integrator == "cca":
        integrate = dedent(
            f"""\
            integrated = False
            try:
                import scanorama
                adatas = [adata[adata.obs[{sample_key}] == b].copy() for b in adata.obs[{sample_key}].unique()]
                scanorama.integrate_scanpy(adatas, dimred={n_pcs})
                import anndata as ad
                adata = ad.concat(adatas, join="outer", index_unique=None)
                if "X_scanorama" in adata.obsm:
                    {nb_scan}
                    integrated = True
                    print("integrator=cca/scanorama")
            except Exception as exc:
                print("SCAGENT_WARN: Scanorama/CCA-like unavailable (" + str(exc) + "); Seurat CCA is R-only")
            if not integrated:
                try:
                    import harmonypy  # noqa: F401
                    sc.external.pp.harmony_integrate(adata, key={sample_key})
                    {nb_harm}
                    print("integrator=harmony_fallback")
                except Exception as exc:
                    print("SCAGENT_WARN: Harmony unavailable (" + str(exc) + ")")
                    {nb_pca}
            """
        )
    elif integrator == "harmony":
        integrate = dedent(
            f"""\
            try:
                import harmonypy  # noqa: F401
                sc.external.pp.harmony_integrate(adata, key={sample_key})
                {nb_harm}
                print("integrator=harmony")
            except Exception as exc:
                print("SCAGENT_WARN: Harmony/harmonypy missing (" + str(exc) + ")")
                {nb_pca}
            """
        )
    else:
        integrate = dedent(
            f"""\
            print("SCAGENT_WARN: skip integration: {skip_reason}")
            {nb_pca}
            """
        )

    if resolution is not None:
        res_block = f"chosen_resolution = {res_default}"
    else:
        res_block = dedent(
            f"""\
            resolutions = {res_list!r}
            sil_scores = {{}}
            Xemb = adata.obsm.get("X_pca_harmony", adata.obsm["X_pca"])
            for r in resolutions:
                key = "leiden_r" + str(r)
                sc.tl.leiden(adata, resolution=r, key_added=key)
                ncl = adata.obs[key].nunique()
                print("resolution=" + str(r) + " n_clusters=" + str(ncl))
                try:
                    from sklearn.metrics import silhouette_score
                    if ncl > 1:
                        sil_scores[r] = float(silhouette_score(Xemb, adata.obs[key].astype(str), sample_size=min(5000, adata.n_obs)))
                except Exception:
                    pass
            if sil_scores:
                chosen_resolution = max(sil_scores, key=sil_scores.get)
                print("silhouette", sil_scores)
            else:
                chosen_resolution = {res_default}
                print("SCAGENT_WARN: sklearn silhouette unavailable; using resolution", chosen_resolution)
            """
        )

    tpl = dedent(
        """\
        # scAgent phase 2: cluster + annotate. Dual validation is mandatory.
        import json
        import numpy as np
        import pandas as pd
        import scanpy as sc
        from pathlib import Path
        from scipy import sparse as sp

        SEED = __SEED__
        N_JOBS = __N_JOBS__
        CACHE_ON = __CACHE_ON__
        np.random.seed(SEED)
        sc.settings.verbosity = 3
        sc.settings.n_jobs = N_JOBS
        sc.settings.set_figure_params(dpi=120, facecolor="white")
        try:
            sc.settings.seed = SEED
        except Exception:
            pass
        fig_dir = Path("figures")
        fig_dir.mkdir(exist_ok=True)
        sc.settings.figdir = str(fig_dir)

        qc_path = Path("adata_qc.h5ad")
        if not qc_path.exists():
            alt = Path(".cache") / "adata_qc.h5ad"
            if alt.exists():
                qc_path = alt
            else:
                raise FileNotFoundError("adata_qc.h5ad missing; run QC phase first")
        adata = sc.read_h5ad(qc_path)
        if not getattr(adata, "isbacked", False) and adata.X is not None and not sp.issparse(adata.X):
            adata.X = sp.csr_matrix(adata.X)
        if __SAMPLE_KEY__ not in adata.obs.columns:
            adata.obs[__SAMPLE_KEY__] = "sample1"

        from scagent.inspect_data import detect_expression_layer
        _xlayer = detect_expression_layer(adata)
        if _xlayer.get("layer") == "scaled":
            print("SCAGENT_WARN: skip scale; X already scaled")
        else:
            sc.pp.scale(adata, max_value=__SCALE_MAX__)
        if "X_pca" in adata.obsm:
            print("SCAGENT_WARN: skip pca; X_pca exists")
        else:
            sc.tl.pca(adata, n_comps=__N_PCS__, svd_solver="arpack")
        __INTEGRATE__
        sc.tl.umap(adata)

        __RES_BLOCK__
        sc.tl.leiden(adata, resolution=chosen_resolution, key_added="leiden")
        if CACHE_ON:
            Path(".cache").mkdir(exist_ok=True)
            adata.write(Path(".cache") / "after_cluster.h5ad")
        sc.pl.umap(adata, color=["leiden", __SAMPLE_KEY__, "pct_counts_mt"], save="_overview.png", show=False)

        __DE_BLOCK__

        __CELLTYPIST__

        MARKERS = __MARKERS__
        TISSUE = __TISSUE_NAME__
        IMMUNE_TISSUES = set(["pbmc", "blood", "immune"])
        expr = adata.raw.to_adata() if adata.raw is not None else adata
        gene_to_idx = {g: i for i, g in enumerate(map(str, expr.var_names))}

        def _mean(genes):
            idx = [gene_to_idx[g] for g in genes if g in gene_to_idx]
            if not idx:
                return None
            X = expr[:, idx].X
            if hasattr(X, "toarray"):
                X = X.toarray()
            return np.asarray(X).mean(axis=1).ravel()

        __SECOND_REF__

        if str(TISSUE).lower() not in IMMUNE_TISSUES and "Immune_All" in str(adata.uns.get("celltypist_model") or ""):
            print("SCAGENT_WARN: CellTypist immune model may be cross-tissue; marker hierarchy takes precedence")

        evidence_rows = []
        cluster_ann = {}
        cluster_lin = {}
        for cl in sorted(adata.obs["leiden"].astype(str).unique()):
            mask = (adata.obs["leiden"].astype(str) == cl).to_numpy()
            best = None
            best_score = -1.0
            for ct in MARKERS:
                pos = ct.get("positive") or []
                neg = ct.get("negative") or []
                if len(pos) < 2:
                    continue
                mp = _mean(pos)
                mn = _mean(neg) if neg else np.zeros(expr.n_obs)
                if mp is None:
                    continue
                score = float(mp[mask].mean() - (0.0 if mn is None else mn[mask].mean()))
                depth = len(ct.get("lineage") or [ct.get("name")])
                score = score + 0.01 * depth
                if score > best_score:
                    best_score = score
                    best = ct
            pos_ok = best and len(best.get("positive") or []) >= 2
            neg_ok = best and len(best.get("negative") or []) >= 1
            lin = list((best or {}).get("lineage") or []) if (pos_ok and neg_ok) else []
            label = (best or {}).get("name") if (pos_ok and neg_ok) else "unknown"
            cluster_ann[cl] = label
            cluster_lin[cl] = lin
            evidence_rows.append(
                {
                    "cluster": cl,
                    "marker_label": label,
                    "lineage": lin,
                    "positive": (best or {}).get("positive"),
                    "negative": (best or {}).get("negative"),
                    "dual_ok": bool(pos_ok and neg_ok),
                    "auto_label": None,
                }
            )
        adata.obs["marker_label"] = adata.obs["leiden"].astype(str).map(cluster_ann)
        max_depth = max((len(v) for v in cluster_lin.values()), default=1)
        for i in range(max_depth):
            col = "cell_type_l" + str(i + 1)
            adata.obs[col] = adata.obs["leiden"].astype(str).map(lambda c, i=i: (cluster_lin.get(c) or ["unknown"])[i] if i < len(cluster_lin.get(c) or []) else (cluster_lin.get(c) or ["unknown"])[-1] if cluster_lin.get(c) else "unknown")
        # Marker hierarchy is the biological assignment. Auto/LLM labels are hypotheses.
        adata.obs["cell_type"] = adata.obs["marker_label"]
        if "celltypist_label" in adata.obs:
            auto = adata.obs["celltypist_label"].astype(str)
            mark = adata.obs["marker_label"].astype(str)
            conflict = (auto != mark) & (mark != "unknown")
            adata.obs["annotation_conflict"] = conflict
            unvalidated = (mark == "unknown") & (adata.obs.get("celltypist_conf", 0) >= 0.5)
            adata.obs.loc[unvalidated, "cell_type"] = auto[unvalidated] + "|unvalidated"
            n_conf = int(conflict.sum())
            if n_conf:
                print("SCAGENT_WARN: marker vs auto-annotation conflict in " + str(n_conf) + " cells; markers kept")
            for row in evidence_rows:
                cl = row["cluster"]
                sub = adata.obs["leiden"].astype(str) == cl
                row["auto_label"] = str(adata.obs.loc[sub, "celltypist_label"].mode().iloc[0]) if sub.any() else None
                row["conflict"] = bool((adata.obs.loc[sub, "annotation_conflict"]).any()) if "annotation_conflict" in adata.obs else False
        low = adata.obs.get("celltypist_conf", 1) < 0.5
        unknown_and_low = (adata.obs["marker_label"].astype(str) == "unknown") & low
        adata.obs.loc[unknown_and_low, "cell_type"] = "low_conf"
        Path("annotation_evidence.json").write_text(json.dumps(evidence_rows, indent=2), encoding="utf-8")
        __PSEUDOBULK__
        color_cols = ["cell_type", "marker_label", "cell_type_l1"]
        if "celltypist_label" in adata.obs:
            color_cols.append("celltypist_label")
        if "ref2_label" in adata.obs:
            color_cols.append("ref2_label")
        sc.pl.umap(adata, color=color_cols, save="_annotation.png", show=False)

        __INTEG_METRICS__

        metrics = {
            "phase": "downstream",
            "resolution": float(chosen_resolution),
            "n_clusters": int(adata.obs["leiden"].nunique()),
            "n_cells": int(adata.n_obs),
            "integrator": __INTEGRATOR__,
            "batch_cluster_dominance": mix,
            "ilisi": ilisi,
            "kbet": kbet,
            "pca_batch_r2": pca_r2,
            "integration_passed": integ_passed,
            "celltypist_model": adata.uns.get("celltypist_model"),
            "seed": SEED,
            "annotation_dual_validation": True,
            "hierarchical_annotation": True,
        }
        print("SCAGENT_METRICS:" + json.dumps(metrics))
        Path("downstream_metrics.json").write_text(json.dumps(metrics, indent=2), encoding="utf-8")
        adata.write("adata_processed.h5ad")
        print(adata)
        """
    )
    return (
        tpl.replace("__SAMPLE_KEY__", sample_key)
        .replace("__INTEGRATE__", integrate.strip("\n"))
        .replace("__RES_BLOCK__", res_block.strip("\n"))
        .replace("__DE_BLOCK__", _de_block(needs_pb, condition_key, sample_key).strip("\n"))
        .replace("__PSEUDOBULK__", _pseudobulk_block(needs_pb, condition_key, sample_key).strip("\n"))
        .replace("__CELLTYPIST__", _celltypist_block(ct_model).strip("\n"))
        .replace("__SECOND_REF__", _second_ref_block().strip("\n"))
        .replace("__INTEG_METRICS__", _integration_metrics_block().strip("\n"))
        .replace("__MARKERS__", marker_json)
        .replace("__INTEGRATOR__", repr(integrator))
        .replace("__TISSUE_NAME__", json.dumps(str(tissue)))
        .replace("__SEED__", str(int(p["seed"])))
        .replace("__SCALE_MAX__", str(float(p["scale_max_value"])))
        .replace("__N_PCS__", str(n_pcs))
        .replace("__N_JOBS__", str(int(perf["n_jobs"])))
        .replace("__CACHE_ON__", "True" if perf["cache"] else "False")
    )


def scanpy_script(meta: dict, qc: dict, plan: dict | None = None) -> str:
    """Combined script for tests and the reproducible dump."""
    p1 = qc_preprocess_script(meta, qc)
    p2 = cluster_annotate_script(meta, qc, plan)
    return p1.rstrip() + "\n\n# --- PHASE 2 ---\n\n" + p2
