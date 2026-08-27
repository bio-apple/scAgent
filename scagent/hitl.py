"""Human-in-the-loop decision cards: MT filter and Leiden resolution."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from scagent.config import analysis_params, load_config, resolve_path
from scagent.logutil import get_logger

log = get_logger("hitl")

_HIGH_MT_TISSUES = {"heart", "tumor", "kidney", "liver"}


def need_hitl(state: dict) -> bool:
    return bool(state.get("interrupt_after_qc") and not state.get("auto_confirm"))


def has_mt_choice(state: dict) -> bool:
    return bool(str(state.get("qc_choice") or "").strip())


def has_resolution_choice(state: dict) -> bool:
    if str(state.get("resolution_choice") or "").strip():
        return True
    return state.get("resolution") is not None


def _median(xs: list[float]) -> float:
    if not xs:
        return 0.0
    s = sorted(xs)
    n = len(s)
    mid = n // 2
    return float(s[mid]) if n % 2 else 0.5 * (s[mid - 1] + s[mid])


def _mad(xs: list[float]) -> float:
    if not xs:
        return 0.0
    med = _median(xs)
    return _median([abs(x - med) for x in xs])


def _pct(xs: list[float], q: float) -> float:
    if not xs:
        return 0.0
    s = sorted(xs)
    i = min(len(s) - 1, max(0, int(round((q / 100.0) * (len(s) - 1)))))
    return float(s[i])


def _frac_above_mad(values: list[float], nmads: int) -> float:
    if not values:
        return 0.0
    med = _median(values)
    mad = _mad(values)
    if mad <= 0:
        return 0.0
    thr = med + nmads * mad
    return 100.0 * sum(1 for v in values if v > thr) / len(values)


def sample_pct_mt(
    data_path: str | None,
    *,
    species: str = "human",
    max_cells: int = 2000,
) -> dict[str, Any]:
    """Sample pct_counts_mt from h5ad. Missing files/deps → empty values, not fatal."""
    empty: dict[str, Any] = {"values": [], "n_sampled": 0, "n_mt_genes": 0, "note": ""}
    if not data_path or not Path(data_path).exists():
        empty["note"] = "无数据文件，选项来自组织 profile。"
        return empty
    try:
        import numpy as np
        import anndata as ad
    except ImportError:
        empty["note"] = "未安装 anndata，无法画线粒体分布。"
        return empty
    try:
        adata = ad.read_h5ad(data_path, backed="r")
        n = min(int(adata.n_obs), max_cells)
        names = [str(g) for g in adata.var_names]
        if (species or "").lower() == "mouse":
            mt = [g.startswith("mt-") or g.upper().startswith("MT-") for g in names]
        else:
            mt = [g.upper().startswith("MT-") for g in names]
        n_mt = int(sum(mt))
        if n_mt < 1:
            try:
                adata.file.close()
            except Exception:
                pass
            return {**empty, "note": "未检测到 MT 基因。"}
        view = adata[:n]
        try:
            X = view.to_memory().X
        except Exception:
            X = view.X
        try:
            adata.file.close()
        except Exception:
            pass
        if hasattr(X, "toarray"):
            arr = np.asarray(X.toarray(), dtype=np.float64)
        else:
            arr = np.asarray(X, dtype=np.float64)
        mt_idx = np.array(mt, dtype=bool)
        tot = arr.sum(axis=1)
        mt_sum = arr[:, mt_idx].sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            pct = np.where(tot > 0, 100.0 * mt_sum / tot, 0.0)
        values = [float(v) for v in np.asarray(pct).ravel().tolist()]
        return {
            "values": values,
            "n_sampled": len(values),
            "n_mt_genes": n_mt,
            "median": _median(values),
            "p98": _pct(values, 98),
            "note": "",
        }
    except Exception as exc:
        empty["note"] = f"读取 pct_mt 失败: {exc}"
        return empty


def _hist_bins(values: list[float], n_bins: int = 24) -> tuple[list[int], list[float]]:
    if not values:
        return [], []
    lo, hi = min(values), max(values)
    if hi <= lo:
        hi = lo + 1.0
    width = (hi - lo) / n_bins
    counts = [0] * n_bins
    for v in values:
        i = min(n_bins - 1, max(0, int((v - lo) / width)))
        counts[i] += 1
    edges = [lo + i * width for i in range(n_bins + 1)]
    return counts, edges


def _svg_bars(counts: list[int], edges: list[float], xlabel: str, width: int = 640, height: int = 200) -> str:
    if not counts:
        return (
            f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">'
            f'<text x="24" y="100" fill="#78716c">无分布数据</text></svg>'
        )
    pad_l, pad_b, pad_t, pad_r = 44, 36, 12, 12
    inner_w = width - pad_l - pad_r
    inner_h = height - pad_t - pad_b
    peak = max(counts) or 1
    bw = inner_w / len(counts)
    bars = []
    for i, c in enumerate(counts):
        h = inner_h * (c / peak)
        x = pad_l + i * bw
        y = pad_t + inner_h - h
        bars.append(
            f'<rect x="{x:.1f}" y="{y:.1f}" width="{max(bw - 1, 0.5):.1f}" height="{h:.1f}" fill="#0f766e"/>'
        )
    xmin, xmax = edges[0], edges[-1]
    axis = (
        f'<line x1="{pad_l}" y1="{pad_t + inner_h}" x2="{width - pad_r}" y2="{pad_t + inner_h}" stroke="#1c1917"/>'
        f'<text x="{pad_l}" y="{height - 8}" fill="#78716c" font-size="11">{xmin:.1f}</text>'
        f'<text x="{width - pad_r - 48}" y="{height - 8}" fill="#78716c" font-size="11">{xmax:.1f}</text>'
        f'<text x="{width / 2 - 40}" y="{height - 8}" fill="#44403c" font-size="11">{xlabel}</text>'
    )
    return (
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" '
        f'font-family="IBM Plex Sans, sans-serif">{"".join(bars)}{axis}</svg>'
    )


def _profile(tissue: str) -> dict:
    cfg = load_config()
    profiles = cfg.get("qc_profiles") or {}
    key = (tissue or "default").lower()
    return dict(profiles.get(key) or profiles.get("default") or {})


def pick_option(decision: dict, choice: Any) -> dict:
    options = list(decision.get("options") or [])
    rec = str(decision.get("recommended") or "")
    fallback = next((o for o in options if o.get("id") == rec), options[0] if options else {})
    if choice is None or str(choice).strip() == "":
        return fallback
    s = str(choice).strip().lower()
    for o in options:
        if str(o.get("id") or "").lower() == s:
            return o
        if str(o.get("resolution")) == s:
            return o
        if str(o.get("nmads")) == s:
            return o
    try:
        num = float(s)
    except ValueError:
        return fallback
    for o in options:
        if o.get("resolution") is not None and abs(float(o["resolution"]) - num) < 1e-9:
            return o
        if o.get("nmads") is not None and abs(float(o["nmads"]) - num) < 1e-9:
            return o
    if "resolution" in (options[0] if options else {}):
        return {
            "id": "custom",
            "label": f"自定义 resolution={num:g}",
            "resolution": num,
            "reason": "湿实验人员指定，不在预设三档内。",
        }
    return fallback


def mt_override(decision: dict, choice: Any) -> dict:
    opt = pick_option(decision, choice)
    out: dict[str, Any] = {}
    if opt.get("nmads") is not None:
        out["nmads"] = int(opt["nmads"])
    return out


def build_mt_decision(state: dict | None = None, *, sample: dict | None = None) -> dict[str, Any]:
    state = state or {}
    meta = state.get("metadata") or {}
    tissue = str(state.get("tissue") or meta.get("tissue") or "default").lower()
    prof = _profile(tissue)
    rec = int(prof.get("nmads") or (6 if tissue in _HIGH_MT_TISSUES else 5))
    sample = sample if sample is not None else sample_pct_mt(
        state.get("data_path") or meta.get("data_path"),
        species=str(meta.get("species") or "human"),
    )
    values = list(sample.get("values") or [])
    lenient_n = rec + 2
    strict_n = max(3, rec - 1)
    note = prof.get("pct_mt_note") or "先看分布再定阈值，禁止默认 mito%<5。"
    stats = ""
    if values:
        stats = (
            f"抽样 {sample.get('n_sampled')} 细胞：median={_median(values):.2f}%，"
            f"P98={_pct(values, 98):.2f}%。"
        )
    options = [
        {
            "id": "lenient",
            "label": "保守（少滤）",
            "nmads": lenient_n,
            "est_remove_pct": round(_frac_above_mad(values, lenient_n), 2),
            "reason": (
                f"MAD n={lenient_n}，只切更极端的右尾。适合心肌/肿瘤/高代谢组织，"
                f"避免把真实高线粒体细胞当死亡细胞丢掉（Yates 2025）。{note}"
            ),
        },
        {
            "id": "recommended",
            "label": "推荐（组织 MAD）",
            "nmads": rec,
            "est_remove_pct": round(_frac_above_mad(values, rec), 2),
            "reason": (
                f"按 {tissue} profile 的单侧高 MAD n={rec}。"
                f"{stats} {note}"
            ).strip(),
        },
        {
            "id": "strict",
            "label": "严格（多滤）",
            "nmads": strict_n,
            "est_remove_pct": round(_frac_above_mad(values, strict_n), 2),
            "reason": (
                f"MAD n={strict_n}，切掉更多高线粒体细胞。适合血液且右尾明显、空液滴多的情况；"
                "过过滤会丢掉稀有/应激细胞，确认前请看直方图。"
            ),
        },
    ]
    counts, edges = _hist_bins(values)
    return {
        "id": "mt",
        "title": "剔除高线粒体细胞",
        "recommended": "recommended",
        "tissue": tissue,
        "xlabel": "pct_counts_mt (%)",
        "histogram": {"counts": counts, "edges": edges},
        "sample": {k: sample.get(k) for k in ("n_sampled", "n_mt_genes", "median", "p98", "note")},
        "options": options,
        "confirm": "scagent confirm mt <lenient|recommended|strict>",
        "note": sample.get("note") or note,
    }


def build_resolution_decision(state: dict | None = None) -> dict[str, Any]:
    state = state or {}
    meta = state.get("metadata") or {}
    tissue = str(state.get("tissue") or meta.get("tissue") or "default").lower()
    n_cells = meta.get("n_cells")
    cfg = load_config()
    res_list = [float(x) for x in (analysis_params(cfg).get("leiden_resolutions") or [0.2, 0.4, 0.6, 0.8, 1.0])]
    if not res_list:
        res_list = [0.2, 0.6, 1.0]
    sil = ((state.get("artifacts") or {}).get("metrics") or {}).get("silhouette") or {}
    mid = 0.6 if 0.6 in res_list else sorted(res_list)[len(res_list) // 2]
    if isinstance(sil, dict) and sil:
        try:
            mid = float(max(sil, key=lambda k: sil[k]))
        except Exception:
            pass
    trio: list[float] = []
    for r in (min(res_list), mid, max(res_list), 0.2, 0.6, 1.0, *res_list):
        rf = float(r)
        if rf not in trio:
            trio.append(rf)
        if len(trio) >= 3:
            break
    coarse, recommended, fine = trio[0], trio[min(1, len(trio) - 1)], trio[-1]
    n_note = f"当前约 {n_cells} 细胞。" if n_cells else ""
    templates = [
        (
            "coarse",
            "粗粒度（谱系）",
            coarse,
            f"Leiden resolution={coarse:g}：T/B/髓系等谱系级分群，适合总览组成。{n_note}欠聚类时再提高。",
        ),
        (
            "recommended",
            "推荐（注释粒度）",
            recommended,
            f"Leiden resolution={recommended:g}：标准注释粒度，平衡过聚类与合并不同类型。不要默认 0.8。",
        ),
        (
            "fine",
            "细粒度（亚型）",
            fine,
            f"Leiden resolution={fine:g}：Tex/naive 等亚型。过聚类会产生无 marker 的碎片簇，仅在关心状态细分时选用。",
        ),
    ]
    options = []
    seen: set[float] = set()
    for oid, label, res, reason in templates:
        if res in seen:
            continue
        seen.add(res)
        options.append({"id": oid, "label": label, "resolution": res, "reason": reason})
    counts = [1] * len(options)
    edges = [float(o["resolution"]) for o in options]
    if len(edges) == 1:
        edges.append(edges[0] + 0.2)
    else:
        step = (edges[-1] - edges[0]) / max(len(edges) - 1, 1)
        edges = edges + [edges[-1] + step]
    xlabel = "Leiden resolution（柱=预设档）"
    if isinstance(sil, dict) and sil:
        counts = [
            max(1, int(round(10 * float(sil.get(o["resolution"], sil.get(str(o["resolution"]), 0.1))))) )
            for o in options
        ]
        xlabel = "silhouette（越高越好）"
    rec_id = "recommended" if any(o["id"] == "recommended" for o in options) else options[min(1, len(options) - 1)]["id"]
    return {
        "id": "resolution",
        "title": "确定聚类 Resolution",
        "recommended": rec_id,
        "tissue": tissue,
        "xlabel": xlabel,
        "histogram": {"counts": counts, "edges": edges},
        "options": options[:3],
        "confirm": "scagent confirm resolution <coarse|recommended|fine|0.4>",
        "note": "Resolution 不是固定 0.8。确认后才进入注释。",
    }


def render_decision_html(decision: dict) -> str:
    counts = (decision.get("histogram") or {}).get("counts") or []
    edges = (decision.get("histogram") or {}).get("edges") or []
    svg = _svg_bars(counts, edges, str(decision.get("xlabel") or ""))
    rec = str(decision.get("recommended") or "")
    cards = []
    for o in decision.get("options") or []:
        mark = " 推荐" if o.get("id") == rec else ""
        extra = ""
        if o.get("nmads") is not None:
            extra = f"MAD n={o['nmads']}"
            if o.get("est_remove_pct") is not None:
                extra += f" · 估计高 MT 剔除 {o['est_remove_pct']}%"
        if o.get("resolution") is not None:
            extra = f"resolution={o['resolution']}"
        cmd = f"scagent confirm {decision.get('id')} {o.get('id')}"
        cards.append(
            "<article class='opt'>"
            f"<h2>{o.get('label')}{mark}</h2>"
            f"<p class='meta'>{extra}</p>"
            f"<p>{o.get('reason')}</p>"
            f"<code>{cmd}</code>"
            "</article>"
        )
    note = decision.get("note") or ""
    return f"""<!DOCTYPE html>
<html lang="zh"><head><meta charset="utf-8"/>
<title>{decision.get('title') or 'HITL'}</title>
<style>
:root {{ --bg:#f4f1ea; --ink:#1c1917; --muted:#78716c; --acc:#0f766e; --card:#fffcf7; }}
body {{ margin:0; font:14px/1.45 "IBM Plex Sans","Source Han Sans SC",sans-serif; background:var(--bg); color:var(--ink); }}
header {{ padding:16px 20px; border-bottom:1px solid #e7e0d4; }}
h1 {{ font-size:18px; margin:0 0 6px; }}
.wrap {{ padding:16px 20px; max-width:960px; }}
.grid {{ display:grid; gap:12px; grid-template-columns:repeat(auto-fit,minmax(240px,1fr)); }}
.opt {{ background:var(--card); padding:12px 14px; border:1px solid #e7e0d4; }}
.opt h2 {{ font-size:14px; margin:0 0 6px; }}
.meta {{ color:var(--acc); font-variant-numeric:tabular-nums; }}
code {{ display:block; background:#1c1917; color:#fffcf7; padding:6px 8px; font-size:12px; }}
.muted {{ color:var(--muted); }}
</style></head>
<body>
<header>
  <h1>{decision.get('title')}</h1>
  <p class="muted">湿实验确认后再继续。推荐项已标出。{note}</p>
</header>
<div class="wrap">
  {svg}
  <p class="muted">{decision.get('xlabel') or ''}</p>
  <div class="grid">{"".join(cards)}</div>
  <p class="muted">{decision.get('confirm')}</p>
</div>
</body></html>
"""


def decisions_dir(cfg: dict | None = None) -> Path:
    out = resolve_path(cfg or load_config(), "outputs") / "decisions"
    out.mkdir(parents=True, exist_ok=True)
    return out


def write_decision_card(decision: dict, *, cfg: dict | None = None) -> Path:
    d = decisions_dir(cfg)
    kind = str(decision.get("id") or "decision")
    (d / f"{kind}.json").write_text(json.dumps(decision, ensure_ascii=False, indent=2), encoding="utf-8")
    html = d / f"{kind}.html"
    html.write_text(render_decision_html(decision), encoding="utf-8")
    log.info("HITL card %s → %s", kind, html)
    return html


def load_session(cfg: dict | None = None) -> dict[str, Any]:
    path = decisions_dir(cfg) / "session.json"
    if not path.is_file():
        return {}
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        return data if isinstance(data, dict) else {}
    except Exception:
        return {}


def save_session(state: dict, cfg: dict | None = None) -> Path:
    d = decisions_dir(cfg)
    prev = load_session(cfg)
    payload = {
        **prev,
        "user_query": state.get("user_query") or prev.get("user_query") or "",
        "data_path": state.get("data_path") or prev.get("data_path") or "",
        "tissue": state.get("tissue") or prev.get("tissue") or "default",
        "thread_id": state.get("thread_id") or prev.get("thread_id"),
        "language": state.get("language") or prev.get("language"),
        "markers_path": state.get("markers_path") or prev.get("markers_path"),
        "batch_key": state.get("batch_key") or prev.get("batch_key"),
        "qc_choice": state.get("qc_choice") if state.get("qc_choice") is not None else prev.get("qc_choice"),
        "resolution_choice": state.get("resolution_choice")
        if state.get("resolution_choice") is not None
        else prev.get("resolution_choice"),
        "report_lang": state.get("report_lang") or prev.get("report_lang") or "zh",
    }
    path = d / "session.json"
    path.write_text(json.dumps(payload, ensure_ascii=False, indent=2), encoding="utf-8")
    return path
