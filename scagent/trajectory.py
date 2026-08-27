"""Trajectory / fate: judge continuity, then Palantir / DPT+PAGA / scVelo + gene trends."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from scagent.logutil import get_logger

log = get_logger("trajectory")

CONTINUOUS_TISSUES = {
    "embryo",
    "fetal",
    "development",
    "ipsc",
    "ips",
    "hematopoiesis",
    "bone_marrow",
    "marrow",
    "thymus",
    "differentiation",
    "reprogramming",
}
DISCRETE_TISSUES = {"pbmc", "blood", "immune"}
STEM_GENES = ("KIT", "PROM1", "SOX2", "POU5F1", "NANOG", "CD34", "EPCAM")
VELOCITY_LAYERS = ("spliced", "unspliced", "Ms", "Mu")


def inspect_trajectory_hints(
    *,
    tissue: str | None = None,
    layers: list[str] | None = None,
    n_obs: int | None = None,
    query: str | None = None,
) -> dict[str, Any]:
    t = str(tissue or "default").lower()
    layers = [str(x).lower() for x in (layers or [])]
    has_vel = any(x in layers for x in ("spliced", "unspliced"))
    q = str(query or "").lower()
    asked = any(
        k in q
        for k in ("轨迹", "伪时间", "命运", "分化", "palantir", "scvelo", "velocity", "monocle", "pseudotime")
    )
    candidate = bool(asked or has_vel or t in CONTINUOUS_TISSUES)
    if t in DISCRETE_TISSUES and not asked and not has_vel:
        candidate = False
    notes = []
    if has_vel:
        notes.append("检测到 spliced/unspliced，可尝试 scVelo。")
    if t in CONTINUOUS_TISSUES:
        notes.append(f"组织 {t} 常含连续分化，将评估轨迹。")
    if t in DISCRETE_TISSUES and not asked:
        notes.append("外周血/免疫更像离散群体；除非用户要求，不默认拟合命运轴。")
    return {
        "trajectory_candidate": candidate,
        "has_velocity_layers": has_vel,
        "trajectory_tissue": t,
        "trajectory_notes": notes,
        "n_cells": n_obs,
    }


def should_plan_trajectory(meta: dict | None, intent: dict | None, cfg: dict | None = None) -> bool:
    mode = str(((cfg or {}).get("modules") or {}).get("trajectory") or "auto").lower()
    if mode in {"off", "none", "skip"}:
        return False
    if mode in {"force", "on", "always"}:
        return True
    intents = (intent or {}).get("intents") or []
    if "trajectory" in intents:
        return True
    meta = meta or {}
    if meta.get("trajectory_candidate") or meta.get("has_velocity_layers"):
        return True
    t = str(meta.get("tissue") or "").lower()
    return t in CONTINUOUS_TISSUES


def _path_like_score(conn, thresh: float = 0.05) -> float:
    import numpy as np

    A = conn.toarray() if hasattr(conn, "toarray") else np.asarray(conn)
    A = np.array(A, dtype=float)
    if A.ndim != 2 or A.shape[0] < 2:
        return 0.0
    np.fill_diagonal(A, 0)
    deg = (A >= thresh).sum(axis=1)
    return float(np.mean((deg >= 1) & (deg <= 2)))


def assess_continuity(adata, *, force: bool = False) -> dict[str, Any]:
    """Decide whether the neighborhood graph looks like a continuous fate process."""
    reasons: list[str] = []
    score = 0
    n_obs = int(getattr(adata, "n_obs", 0) or 0)
    if n_obs < 40:
        out = {"verdict": "skip", "score": 0, "reasons": ["细胞数过少，伪时间不稳定"], "n_clusters": 0, "path_like": 0.0}
        if force:
            out["verdict"] = "run"
            out["confidence"] = "low"
            out["reasons"].append("用户要求轨迹，低置信仍拟合。")
        return out
    obs = adata.obs
    if "leiden" not in obs:
        return {"verdict": "skip" if not force else "run", "score": 0, "reasons": ["尚无 Leiden"], "confidence": "low" if force else "none"}
    n_cl = int(obs["leiden"].nunique())
    path_like = 0.0
    if n_cl < 3:
        reasons.append(f"簇数={n_cl}<3，更像离散群体")
        score -= 1
    else:
        score += 1
        reasons.append(f"簇数={n_cl}")
    try:
        import scanpy as sc

        if "neighbors" in adata.uns:
            sc.tl.paga(adata, groups="leiden")
            conn = (adata.uns.get("paga") or {}).get("connectivities")
            if conn is not None:
                path_like = _path_like_score(conn)
                if path_like >= 0.45:
                    score += 2
                    reasons.append(f"PAGA 连接偏路径/树 (path_like={path_like:.2f})")
                elif path_like < 0.2:
                    score -= 1
                    reasons.append(f"PAGA 更像团块 (path_like={path_like:.2f})")
                else:
                    reasons.append(f"PAGA path_like={path_like:.2f}")
    except Exception as exc:
        reasons.append(f"PAGA 未完成: {exc}")
    if "X_umap" in getattr(adata, "obsm", {}):
        try:
            import numpy as np

            umap = np.asarray(adata.obsm["X_umap"])
            cents = []
            for cl in obs["leiden"].astype(str).unique():
                cents.append(umap[obs["leiden"].astype(str).to_numpy() == cl].mean(axis=0))
            C = np.vstack(cents)
            C = C - C.mean(axis=0)
            ev = np.linalg.svd(C, compute_uv=False)
            ratio = float(ev[0] / max(ev[1], 1e-8)) if len(ev) > 1 else 1.0
            if ratio >= 2.0:
                score += 1
                reasons.append(f"UMAP 簇中心拉长 (ratio={ratio:.1f})")
        except Exception:
            pass
    layers = [str(x).lower() for x in getattr(adata, "layers", {})]
    has_vel = "spliced" in layers and "unspliced" in layers
    if has_vel:
        score += 1
        reasons.append("存在 spliced/unspliced")
    verdict = "run" if (force or score >= 2) else "skip"
    if force and score < 2:
        reasons.append("用户要求轨迹；评估偏离散，结果仅作探索。")
    if verdict == "skip":
        reasons.append("不把离散聚类强行画成命运轴。")
    return {
        "verdict": verdict,
        "score": score,
        "reasons": reasons,
        "n_clusters": n_cl,
        "path_like": path_like,
        "has_velocity_layers": has_vel,
        "confidence": "high" if score >= 3 else ("low" if force else "moderate"),
    }


def _root_index(adata) -> int:
    import numpy as np

    names = [str(g) for g in adata.var_names]
    upper = {g.upper(): i for i, g in enumerate(names)}
    idx = [upper[g] for g in STEM_GENES if g in upper]
    if idx:
        X = adata[:, [names[i] for i in idx]].X
        if hasattr(X, "toarray"):
            X = X.toarray()
        scores = np.asarray(X).mean(axis=1).ravel()
        return int(np.argmax(scores))
    vc = adata.obs["leiden"].astype(str).value_counts()
    root_cl = str(vc.index[-1])
    return int(np.flatnonzero(adata.obs["leiden"].astype(str).to_numpy() == root_cl)[0])


def _fit_dpt_paga(adata, figdir: Path) -> list[str]:
    import scanpy as sc

    from scagent.plotting import apply_figdir

    apply_figdir(figdir)
    methods = ["paga", "dpt"]
    if "paga" not in adata.uns:
        sc.tl.paga(adata, groups="leiden")
    try:
        sc.pl.paga(adata, color="leiden", threshold=0.03, save="_paga.png", show=False)
    except Exception as exc:
        log.info("paga plot skipped: %s", exc)
    sc.tl.diffmap(adata)
    adata.uns["iroot"] = _root_index(adata)
    sc.tl.dpt(adata)
    try:
        sc.pl.umap(adata, color=["dpt_pseudotime", "leiden"], save="_pseudotime.png", show=False)
    except Exception as exc:
        log.info("pseudotime umap skipped: %s", exc)
    return methods


def _try_palantir(adata) -> str | None:
    try:
        import palantir
    except Exception:
        return None
    try:
        early = adata.obs_names[_root_index(adata)]
        if hasattr(palantir.utils, "run_diffusion_maps"):
            palantir.utils.run_diffusion_maps(adata)
        if hasattr(palantir.utils, "determine_multiscale_space"):
            palantir.utils.determine_multiscale_space(adata)
        pr = palantir.core.run_palantir(adata, early, num_waypoints=min(120, max(20, adata.n_obs // 4)))
        pt = getattr(pr, "pseudotime", None)
        if pt is None and isinstance(pr, dict):
            pt = pr.get("pseudotime")
        if pt is not None:
            adata.obs["palantir_pseudotime"] = pt
        bp = getattr(pr, "branch_probs", None)
        if bp is None and isinstance(pr, dict):
            bp = pr.get("branch_probs")
        if bp is not None:
            adata.uns["palantir_branch_probs"] = "present"
        return "palantir"
    except Exception as exc:
        log.info("palantir skipped: %s", exc)
        return None


def _try_scvelo(adata, figdir: Path) -> str | None:
    layers = {str(x).lower() for x in getattr(adata, "layers", {})}
    if "spliced" not in layers or "unspliced" not in layers:
        return None
    try:
        import scvelo as scv
    except Exception:
        return "scvelo_unavailable"
    try:
        from scagent.plotting import apply_figdir

        apply_figdir(figdir)
        scv.pp.moments(adata, n_pcs=min(30, adata.n_vars), n_neighbors=min(15, max(5, adata.n_obs // 10)))
        scv.tl.velocity(adata, mode="stochastic")
        scv.tl.velocity_graph(adata)
        try:
            scv.pl.velocity_embedding_stream(adata, basis="umap", save="_velocity.png", show=False)
        except Exception:
            scv.pl.velocity_embedding(adata, basis="umap", save="_velocity.png", show=False)
        return "scvelo"
    except Exception as exc:
        log.info("scvelo skipped: %s", exc)
        return "scvelo_failed"


def _try_monocle3(workspace: Path) -> str | None:
    from scagent.config import REPO_ROOT

    h5ad = workspace / "adata_processed.h5ad"
    if not h5ad.is_file():
        h5ad = workspace / "adata_qc.h5ad"
    script = REPO_ROOT / "scagent" / "r" / "monocle3.R"
    if not script.is_file() or not h5ad.is_file():
        return None
    import shutil
    import subprocess

    rscript = shutil.which("Rscript")
    if not rscript:
        return None
    out = workspace / "monocle3_status.json"
    try:
        subprocess.run(
            [rscript, str(script), str(h5ad), str(out)],
            check=False,
            capture_output=True,
            text=True,
            timeout=120,
        )
        if out.is_file():
            payload = json.loads(out.read_text(encoding="utf-8"))
            return str(payload.get("status") or "monocle3")
    except Exception as exc:
        log.info("monocle3 skipped: %s", exc)
    return None


def plot_gene_trends(adata, *, time_key: str = "dpt_pseudotime", n_genes: int = 8, figdir: Path | None = None) -> str | None:
    if time_key not in adata.obs:
        return None
    import numpy as np

    figdir = Path(figdir or "figures")
    figdir.mkdir(parents=True, exist_ok=True)
    t = np.asarray(adata.obs[time_key], dtype=float)
    order = np.argsort(t)
    genes = []
    if "highly_variable" in adata.var:
        hv = [str(g) for g, ok in zip(adata.var_names, adata.var["highly_variable"]) if ok]
        genes = hv[:n_genes]
    if len(genes) < 3:
        genes = [str(g) for g in list(adata.var_names)[:n_genes]]
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        return None
    n = min(len(genes), n_genes)
    fig, axes = plt.subplots(n, 1, figsize=(6, 1.6 * n), sharex=True)
    if n == 1:
        axes = [axes]
    expr = adata.raw.to_adata() if adata.raw is not None else adata
    for ax, g in zip(axes, genes[:n]):
        if g not in expr.var_names:
            ax.set_visible(False)
            continue
        y = expr[:, g].X
        if hasattr(y, "toarray"):
            y = y.toarray()
        y = np.asarray(y).ravel()[order]
        tt = t[order]
        ax.scatter(tt, y, s=4, alpha=0.25, c="#0f766e", linewidths=0)
        if len(tt) > 10:
            bins = np.linspace(tt.min(), tt.max(), 12)
            idx = np.digitize(tt, bins)
            mu = [y[idx == i].mean() if np.any(idx == i) else np.nan for i in range(1, len(bins))]
            ax.plot(0.5 * (bins[:-1] + bins[1:]), mu, color="#111827", lw=1.5)
        ax.set_ylabel(g, fontsize=8)
    axes[-1].set_xlabel(time_key)
    fig.tight_layout()
    path = figdir / "gene_trends.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return str(path)


def run_trajectory_phase(adata, *, mode: str = "auto", workspace: str | Path | None = None) -> dict[str, Any]:
    """Assess continuity; if supported (or mode=force) fit PAGA/DPT, Palantir, scVelo, gene trends."""
    ws = Path(workspace or ".")
    figdir = ws / "figures"
    figdir.mkdir(parents=True, exist_ok=True)
    force = str(mode).lower() in {"force", "on", "always", "true", "1"}
    if str(mode).lower() in {"off", "none", "skip"}:
        payload = {"verdict": "skip", "reasons": ["modules.trajectory=off"], "methods": []}
        adata.uns["scagent_trajectory"] = payload
        print("SCAGENT_WARN: trajectory disabled")
        return payload
    assess = assess_continuity(adata, force=force)
    methods: list[str] = []
    if assess.get("verdict") != "run":
        payload = {**assess, "methods": []}
        adata.uns["scagent_trajectory"] = payload
        print("SCAGENT_WARN: no continuous trajectory; skipped fate axis (" + "; ".join(assess.get("reasons") or []) + ")")
        return payload
    try:
        methods.extend(_fit_dpt_paga(adata, figdir))
    except Exception as exc:
        log.info("dpt/paga failed: %s", exc)
        assess["reasons"] = list(assess.get("reasons") or []) + [f"DPT/PAGA 失败: {exc}"]
    pal = _try_palantir(adata)
    if pal:
        methods.append(pal)
    vel = _try_scvelo(adata, figdir)
    if vel:
        methods.append(vel)
    mono = _try_monocle3(ws)
    if mono:
        methods.append("monocle3:" + mono)
    trend = None
    try:
        trend = plot_gene_trends(adata, figdir=figdir)
    except Exception as exc:
        log.info("gene trends skipped: %s", exc)
    payload = {
        **assess,
        "methods": methods,
        "gene_trends": trend,
        "root_index": int(adata.uns.get("iroot") or 0),
    }
    adata.uns["scagent_trajectory"] = payload
    print("trajectory=" + ",".join(methods) + " verdict=" + str(assess.get("verdict")))
    return payload
