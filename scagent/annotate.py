"""Multi-evidence cell-type fusion. No single mapper (Azimuth/CellTypist) is sufficient."""

from __future__ import annotations

from collections import Counter
from typing import Any

from scagent.logutil import get_logger

log = get_logger("annotate")

CELLTYPIST_CONF_THRESHOLD = 0.8
SCANVI_LABELS_KEY = "_scanvi_supervision"
_EMPTY = {"", "unknown", "unassigned", "nan", "none", "low_conf", "mixed", "na"}


def _norm(label: Any) -> str:
    s = str(label or "").strip()
    if s.lower() in _EMPTY or s.lower() == "none":
        return ""
    return s


def labels_agree(a: Any, b: Any) -> bool:
    """True if labels are the same type under light aliasing (substring / case)."""
    x, y = _norm(a).lower(), _norm(b).lower()
    if not x or not y:
        return False
    if x == y:
        return True
    return x in y or y in x


def _canonical(group: list[str]) -> str:
    counts = Counter(group)
    return counts.most_common(1)[0][0]


def deg_catalog_labels(adata, catalog: dict, *, groupby: str = "leiden", n_top: int = 25) -> None:
    """Assign deg_label from Wilcoxon top genes ∩ catalog positive markers (independent of dual score)."""
    import numpy as np
    import pandas as pd

    types = list((catalog or {}).get("cell_types") or catalog or [])
    if isinstance(catalog, dict) and "cell_types" in catalog:
        types = catalog["cell_types"]
    adata.obs["deg_label"] = "unknown"
    rg = (adata.uns or {}).get("rank_genes_groups") or {}
    names = rg.get("names")
    if names is None or groupby not in adata.obs:
        log.info("deg_catalog_labels skipped: no rank_genes_groups")
        return
    try:
        genes_by_group = pd.DataFrame(names)
    except Exception:
        log.info("deg_catalog_labels skipped: names not tabular")
        return
    mapping: dict[str, str] = {}
    for cl in genes_by_group.columns.astype(str):
        top = [str(g) for g in genes_by_group[cl].head(n_top).tolist() if str(g) != "nan"]
        top_set = {g.upper() for g in top}
        best, best_hit = "unknown", 0.0
        for ct in types:
            pos = [str(g).upper() for g in (ct.get("positive") or [])]
            if len(pos) < 2:
                continue
            hit = sum(g in top_set for g in pos) / len(pos)
            if hit > best_hit:
                best_hit = hit
                best = str(ct.get("name") or "unknown")
        mapping[cl] = best if best_hit >= 0.25 else "unknown"
    adata.obs["deg_label"] = adata.obs[groupby].astype(str).map(lambda c: mapping.get(c, "unknown"))
    n_ok = int((np.asarray(adata.obs["deg_label"].astype(str)) != "unknown").sum())
    log.info("deg_catalog_labels assigned non-unknown=%s / %s", n_ok, adata.n_obs)


def fuse_annotation(
    adata,
    *,
    sources: tuple[str, ...] = ("marker_label", "celltypist_label", "deg_label"),
    marker_col: str = "marker_label",
    min_agree: int = 2,
):
    """Majority vote across independent evidence. Azimuth/CellTypist alone never sets cell_type."""
    import pandas as pd

    n = adata.n_obs
    cols = [c for c in sources if c in adata.obs]
    status = ["insufficient"] * n
    fused = ["unknown"] * n
    n_agree = [0] * n
    n_src = [0] * n
    mark = adata.obs[marker_col].astype(str) if marker_col in adata.obs else pd.Series(["unknown"] * n, index=adata.obs.index)

    for i, idx in enumerate(adata.obs.index):
        votes: list[str] = []
        for c in cols:
            lab = _norm(adata.obs.at[idx, c] if c in adata.obs.columns else "")
            if lab:
                votes.append(lab)
        n_src[i] = len(votes)
        if not votes:
            continue
        buckets: list[list[str]] = []
        for v in votes:
            placed = False
            for b in buckets:
                if labels_agree(v, b[0]):
                    b.append(v)
                    placed = True
                    break
            if not placed:
                buckets.append([v])
        buckets.sort(key=len, reverse=True)
        best = buckets[0]
        n_agree[i] = len(best)
        winner = _canonical(best)
        mk = _norm(mark.loc[idx] if hasattr(mark, "loc") else mark[i])
        marker_in_best = bool(mk) and any(labels_agree(mk, x) for x in best)
        autos_vs_marker = bool(mk) and n_agree[i] >= min_agree and not marker_in_best and mk not in _EMPTY

        if mk and n_agree[i] >= min_agree and marker_in_best:
            fused[i], status[i] = mk, "fused"
        elif autos_vs_marker:
            fused[i], status[i] = "mixed", "conflict_mixed"
        elif mk:
            fused[i], status[i] = mk, "marker_only"
        elif n_agree[i] >= min_agree:
            fused[i], status[i] = winner + "|unvalidated", "auto_consensus"
        else:
            fused[i], status[i] = "low_conf", "insufficient"

    adata.obs["cell_type"] = fused
    adata.obs["annotation_status"] = status
    adata.obs["annotation_n_agree"] = n_agree
    adata.obs["annotation_n_sources"] = n_src
    n_fused = sum(s == "fused" for s in status)
    n_mix = sum(s == "conflict_mixed" for s in status)
    log.info("fuse_annotation fused=%s mixed=%s / %s", n_fused, n_mix, n)
    print(
        "annotation_fusion fused="
        + str(n_fused)
        + " mixed="
        + str(n_mix)
        + " sources="
        + ",".join(cols)
    )
    return adata


def _ensure_counts_layer(adata):
    import numpy as np

    if "counts" in adata.layers:
        return "counts"
    if adata.raw is not None:
        adata.layers["counts"] = adata.raw.X.copy()
        return "counts"
    x = adata.X
    if hasattr(x, "max"):
        try:
            mx = float(x.max())
        except Exception:
            mx = None
        if mx is not None and mx > 50:
            adata.layers["counts"] = adata.X.copy()
            return "counts"
    raise RuntimeError("scANVI needs raw counts in layers['counts'] or adata.raw")


def _run_scanvi(
    adata,
    *,
    labels_key: str,
    batch_key: str | None,
    unlabeled: str = "Unknown",
    max_epochs_scvi: int = 50,
    max_epochs_scanvi: int = 20,
):
    import scvi

    layer = _ensure_counts_layer(adata)
    setup_kw: dict[str, Any] = {"layer": layer, "labels_key": labels_key, "unlabeled_category": unlabeled}
    if batch_key and batch_key in adata.obs.columns and int(adata.obs[batch_key].nunique()) > 1:
        setup_kw["batch_key"] = batch_key
    scvi.model.SCANVI.setup_anndata(adata, **setup_kw)
    scvi_model = scvi.model.SCVI(adata, n_latent=30)
    scvi_model.train(max_epochs=max_epochs_scvi, early_stopping=True)
    model = scvi.model.SCANVI.from_scvi_model(scvi_model, unlabeled_category=unlabeled)
    model.train(max_epochs=max_epochs_scanvi)
    pred = model.predict()
    soft = model.predict(soft=True)
    conf = soft.max(axis=1)
    return pred, conf, model


def ensemble_cell_annotation(
    adata,
    *,
    conf_threshold: float = CELLTYPIST_CONF_THRESHOLD,
    sample_key: str | None = "sample",
    labels_key: str = "celltypist_label",
    conf_key: str = "celltypist_conf",
    max_epochs_scvi: int = 50,
    max_epochs_scanvi: int = 20,
):
    """CellTypist first; cells with max_prob < threshold get scANVI semi-supervised labels.

    Writes ``scagent_annotation``, ``scagent_annotation_conf``, ``scagent_annotation_method``.
    """
    import numpy as np

    n = adata.n_obs
    if labels_key not in adata.obs:
        adata.obs["scagent_annotation"] = "unassigned"
        adata.obs["scagent_annotation_conf"] = 0.0
        adata.obs["scagent_annotation_method"] = "none"
        adata.uns["scagent_annotation"] = {"method": "none", "reason": "no_celltypist"}
        return adata

    labels = adata.obs[labels_key].astype(str).to_numpy()
    if conf_key in adata.obs:
        conf = np.asarray(adata.obs[conf_key], dtype=float)
    else:
        conf = np.ones(n, dtype=float)

    high = conf >= float(conf_threshold)
    low = ~high
    info: dict[str, Any] = {
        "celltypist_threshold": float(conf_threshold),
        "n_high_conf": int(high.sum()),
        "n_low_conf": int(low.sum()),
        "scanvi_ran": False,
    }

    if not low.any():
        adata.obs["scagent_annotation"] = labels
        adata.obs["scagent_annotation_conf"] = conf
        adata.obs["scagent_annotation_method"] = "celltypist"
        info["method"] = "celltypist"
        adata.uns["scagent_annotation"] = info
        print("scagent_annotation=celltypist_only n_low_conf=0")
        return adata

    supervise = np.where(high, labels, "Unknown")
    adata.obs[SCANVI_LABELS_KEY] = supervise
    batch_key = sample_key if sample_key and sample_key in adata.obs.columns else None
    try:
        scanvi_pred, scanvi_conf, _model = _run_scanvi(
            adata,
            labels_key=SCANVI_LABELS_KEY,
            batch_key=batch_key,
            max_epochs_scvi=max_epochs_scvi,
            max_epochs_scanvi=max_epochs_scanvi,
        )
        scanvi_pred = np.asarray(scanvi_pred, dtype=object)
        scanvi_conf = np.asarray(scanvi_conf, dtype=float)
        final = labels.copy()
        final[low] = scanvi_pred[low]
        final_conf = conf.copy()
        final_conf[low] = scanvi_conf[low]
        methods = np.where(high, "celltypist", "scanvi")
        adata.obs["scanvi_label"] = scanvi_pred
        adata.obs["scanvi_conf"] = scanvi_conf
        info.update({"scanvi_ran": True, "method": "celltypist+scanvi"})
        adata.obs["scagent_annotation_method"] = methods
        print(
            "scagent_annotation=celltypist+scanvi n_low_conf="
            + str(int(low.sum()))
            + " scanvi_median_conf="
            + str(round(float(np.median(scanvi_conf[low])), 3) if low.any() else 0)
        )
    except Exception as exc:
        log.warning("scANVI fallback failed: %s", exc)
        print("SCAGENT_WARN: scANVI fallback failed (" + str(exc) + "); using CellTypist labels")
        final = labels
        final_conf = conf
        adata.obs["scagent_annotation_method"] = np.where(high, "celltypist", "celltypist_low_conf")
        info.update({"scanvi_ran": False, "method": "celltypist_only", "scanvi_error": str(exc)})

    adata.obs["scagent_annotation"] = final
    adata.obs["scagent_annotation_conf"] = final_conf
    adata.uns["scagent_annotation"] = info
    return adata
