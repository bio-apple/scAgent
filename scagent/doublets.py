"""Doublet detection: Scrublet plus scDblFinder (or count-simulation) cross-check."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from scagent.logutil import get_logger

log = get_logger("doublets")

COMPLEX_TISSUES = {"tumor", "brain", "heart", "kidney", "liver", "embryo"}
_R_SCRIPT = Path(__file__).resolve().parent / "r" / "doublets.R"


def expected_doublet_rate(n_cells: int) -> float:
    """10x rule of thumb: ~0.8% per 1,000 captured cells, capped at 10%."""
    n = max(int(n_cells or 0), 1)
    return float(min(0.10, max(0.005, 0.008 * (n / 1000.0))))


def resolve_doublet_methods(
    requested: str | None,
    *,
    tissue: str | None = None,
    n_samples: int | None = None,
) -> list[str]:
    """Return engines to run. scdblfinder in the list means 'second method, prefer R'."""
    req = str(requested or "auto").lower().strip()
    n = int(n_samples or 1)
    t = str(tissue or "default").lower()
    if req == "scrublet":
        return ["scrublet"]
    if req in {"both", "scdblfinder", "cross"}:
        return ["scrublet", "scdblfinder"]
    if n > 1 or t in COMPLEX_TISSUES:
        return ["scrublet", "scdblfinder"]
    return ["scrublet"]


def _consensus(a, b):
    import numpy as np

    a = np.asarray(a, dtype=bool)
    b = np.asarray(b, dtype=bool)
    both = a & b
    disc = a ^ b
    agree = float(np.mean(a == b)) if len(a) else 1.0
    return both, disc, agree


def _simulate_doublet_scores(adata, *, seed: int = 0):
    """Independent count-simulation scorer (pair-sum doublets in SVD space)."""
    import numpy as np
    from scipy import sparse as sp
    from sklearn.decomposition import TruncatedSVD
    from sklearn.neighbors import NearestNeighbors

    X = adata.X
    if sp.issparse(X):
        X = X.tocsr().astype(np.float64)
    else:
        X = np.asarray(X, dtype=np.float64)
        X = sp.csr_matrix(X)
    n = int(X.shape[0])
    if n < 8:
        return np.zeros(n), np.zeros(n, dtype=bool)
    n_sim = int(min(n, max(64, n // 4)))
    rng = np.random.default_rng(seed)
    i = rng.integers(0, n, size=n_sim)
    j = rng.integers(0, n, size=n_sim)
    same = i == j
    if np.any(same):
        j[same] = (j[same] + 1) % n
    sim = X[i] + X[j]
    n_comp = int(min(20, max(2, n - 1), X.shape[1] - 1))
    svd = TruncatedSVD(n_components=n_comp, random_state=seed)
    z_real = svd.fit_transform(X)
    z_sim = svd.transform(sim)
    k = min(15, n_sim)
    nn = NearestNeighbors(n_neighbors=k)
    nn.fit(z_sim)
    dist, _ = nn.kneighbors(z_real)
    score = 1.0 / (1.0 + dist.mean(axis=1))
    thr = float(np.quantile(score, 1.0 - expected_doublet_rate(n)))
    pred = score >= thr
    return score, pred


def _per_sample_sim(adata, sample_key: str | None):
    import numpy as np

    n = adata.n_obs
    score = np.zeros(n, dtype=np.float64)
    pred = np.zeros(n, dtype=bool)
    if sample_key and sample_key in adata.obs.columns and int(adata.obs[sample_key].nunique()) > 1:
        col = adata.obs[sample_key].astype(str)
        for i, lab in enumerate(col.unique()):
            mask = (col == lab).to_numpy()
            if int(mask.sum()) < 8:
                continue
            s, p = _simulate_doublet_scores(adata[mask], seed=i)
            score[mask] = s
            pred[mask] = p
        return score, pred
    return _simulate_doublet_scores(adata)


def _scdblfinder_rscript(adata, sample_key: str | None) -> dict[str, Any] | None:
    import csv
    import shutil
    import subprocess
    import tempfile

    import numpy as np

    if not shutil.which("Rscript") or not _R_SCRIPT.is_file():
        return None
    sample = sample_key if sample_key and sample_key in adata.obs.columns else None
    if sample and int(adata.obs[sample].nunique()) <= 1:
        sample = None
    tmpdir = tempfile.TemporaryDirectory()
    try:
        h5ad = Path(tmpdir.name) / "dbl.h5ad"
        csv_path = Path(tmpdir.name) / "dbl.csv"
        import anndata as ad

        slim = ad.AnnData(adata.X.copy())
        slim.obs_names = adata.obs_names.to_numpy()
        slim.var_names = adata.var_names.to_numpy()
        if sample:
            slim.obs[sample] = adata.obs[sample].to_numpy()
        slim.write_h5ad(h5ad)
        proc = subprocess.run(
            ["Rscript", str(_R_SCRIPT), str(h5ad), str(csv_path), sample or "NONE"],
            capture_output=True,
            text=True,
            timeout=180,
            check=False,
        )
        if proc.returncode != 0 or not csv_path.is_file():
            log.info("scDblFinder Rscript skipped: %s", (proc.stderr or proc.stdout or "")[-400:])
            return None
        rows = list(csv.DictReader(csv_path.open(encoding="utf-8")))
        if not rows:
            return None
        by = {str(r.get("barcode") or ""): r for r in rows}
        names = [str(x) for x in adata.obs_names]
        score = np.zeros(adata.n_obs, dtype=np.float64)
        pred = np.zeros(adata.n_obs, dtype=bool)
        if all(n in by for n in names):
            for i, n in enumerate(names):
                r = by[n]
                score[i] = float(r.get("score") or 0.0)
                pred[i] = str(r.get("class") or "").lower() == "doublet"
        elif len(rows) == adata.n_obs:
            for i, r in enumerate(rows):
                score[i] = float(r.get("score") or 0.0)
                pred[i] = str(r.get("class") or "").lower() == "doublet"
        else:
            return None
        return {"score": score, "pred": pred}
    except Exception as exc:
        log.info("scDblFinder unavailable: %s", exc)
        return None
    finally:
        tmpdir.cleanup()


def _run_scrublet(adata, sample_key: str | None) -> bool:
    import scanpy as sc

    bk = None
    if sample_key and sample_key in adata.obs.columns and int(adata.obs[sample_key].nunique()) > 1:
        bk = sample_key
    sc.pp.scrublet(adata, batch_key=bk)
    if "predicted_doublet" not in adata.obs:
        raise RuntimeError("scrublet did not write predicted_doublet")
    if "doublet_score" not in adata.obs:
        adata.obs["doublet_score"] = 0.0
    return True


def detect_doublets(
    adata,
    *,
    methods: str | list[str] | None = "auto",
    sample_key: str | None = None,
    tissue: str | None = None,
    remove: bool = False,
    n_samples: int | None = None,
):
    """Scrublet, plus scDblFinder (R) or count-simulation when cross-check is on.

    Consensus `predicted_doublet` is the intersection when two methods succeed.
    Discordant cells are flagged, not auto-removed.
    """
    import numpy as np

    n_samples = n_samples if n_samples is not None else (
        int(adata.obs[sample_key].nunique()) if sample_key and sample_key in adata.obs.columns else 1
    )
    if isinstance(methods, (list, tuple)):
        want = [str(m).lower() for m in methods]
    else:
        want = resolve_doublet_methods(methods, tissue=tissue, n_samples=n_samples)
    need_second = "scdblfinder" in want or "both" in want or "sim" in want
    engines: list[str] = []
    calls: dict[str, Any] = {}
    scores: dict[str, Any] = {}

    adata.obs["predicted_doublet"] = False
    adata.obs["doublet_score"] = 0.0
    try:
        if "scrublet" in want or not want:
            _run_scrublet(adata, sample_key)
            calls["scrublet"] = adata.obs["predicted_doublet"].to_numpy().astype(bool)
            scores["scrublet"] = np.asarray(adata.obs["doublet_score"], dtype=float)
            adata.obs["doublet_scrublet"] = calls["scrublet"]
            engines.append("scrublet")
    except Exception as exc:
        log.warning("scrublet failed: %s", exc)
        print("SCAGENT_WARN: scrublet failed (" + str(exc) + ")")

    if need_second:
        r = None
        if "scdblfinder" in want or "both" in want:
            r = _scdblfinder_rscript(adata, sample_key)
        if r is not None:
            calls["scdblfinder"] = np.asarray(r["pred"], dtype=bool)
            scores["scdblfinder"] = np.asarray(r["score"], dtype=float)
            adata.obs["doublet_scdblfinder"] = calls["scdblfinder"]
            engines.append("scdblfinder")
        else:
            s, p = _per_sample_sim(adata, sample_key)
            calls["sim"] = np.asarray(p, dtype=bool)
            scores["sim"] = np.asarray(s, dtype=float)
            adata.obs["doublet_sim"] = calls["sim"]
            engines.append("sim")
            if "scdblfinder" in want or "both" in want:
                print("SCAGENT_WARN: scDblFinder unavailable; used count-simulation cross-check")

    n = adata.n_obs
    if len(calls) >= 2:
        keys = list(calls)
        both, disc, agree = _consensus(calls[keys[0]], calls[keys[1]])
        pred = both
        score = np.maximum(scores[keys[0]], scores[keys[1]])
        status = "ok"
    elif len(calls) == 1:
        k = next(iter(calls))
        pred = np.asarray(calls[k], dtype=bool)
        score = np.asarray(scores[k], dtype=float)
        disc = np.zeros(n, dtype=bool)
        agree = 1.0
        status = "partial" if need_second else "ok"
    else:
        pred = np.zeros(n, dtype=bool)
        score = np.zeros(n, dtype=float)
        disc = np.zeros(n, dtype=bool)
        agree = 0.0
        status = "failed"

    adata.obs["predicted_doublet"] = pred
    adata.obs["doublet_score"] = score
    adata.obs["doublet_discordant"] = disc
    rate = float(np.mean(pred)) if n else 0.0
    info = {
        "status": status,
        "rate": rate,
        "agreement": round(float(agree), 4),
        "methods": engines,
        "n_discordant": int(np.sum(disc)),
        "remove": bool(remove),
        "n_doublets": int(np.sum(pred)),
    }
    adata.uns["doublets"] = info
    print(
        "doublet_status="
        + status
        + " doublet_rate="
        + str(round(rate, 4))
        + " doublet_methods="
        + ",".join(engines)
        + " doublet_agreement="
        + str(info["agreement"])
    )
    try:
        import scanpy as sc

        sc.pl.violin(adata, ["doublet_score"], save="_doublet_score.png", show=False)
    except Exception:
        print("SCAGENT_WARN: doublet violin skipped")
    if remove and status != "failed":
        n_rm = int(np.sum(pred))
        adata = adata[~pred].copy()
        print("removed_doublets=" + str(n_rm))
    log.info("doublets %s", info)
    return adata
