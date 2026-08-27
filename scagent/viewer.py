"""Interactive Plotly AnnData viewer: lasso/box select cells and ask the agent."""

from __future__ import annotations

import json
from collections import Counter
from pathlib import Path
from typing import Any

from scagent.config import load_config, resolve_path
from scagent.logutil import get_logger

log = get_logger("viewer")

MAX_CELLS = 15_000
_COLOR_COLS = ("leiden", "cell_type", "cell_type_l1", "annotation_status", "sample", "predicted_doublet", "doublet_call", "scagent_annotation")
_QC_COLS = ("n_genes_by_counts", "total_counts", "pct_counts_mt")

_HTML = """<!DOCTYPE html>
<html lang="zh">
<head>
<meta charset="utf-8"/>
<title>scAgent viewer</title>
<script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
<style>
  :root { --bg:#f4f1ea; --ink:#1c1917; --muted:#78716c; --acc:#0f766e; --card:#fffcf7; }
  * { box-sizing: border-box; }
  body { margin:0; font: 14px/1.45 "IBM Plex Sans", "Source Han Sans SC", sans-serif; background:var(--bg); color:var(--ink); }
  header { padding:12px 18px; border-bottom:1px solid #e7e0d4; display:flex; gap:16px; align-items:baseline; }
  header h1 { font-size:16px; margin:0; font-weight:600; }
  header span { color:var(--muted); font-size:12px; }
  #layout { display:grid; grid-template-columns: 1fr 320px; min-height: calc(100vh - 48px); }
  #plots { padding:8px; }
  #umap, #violin { width:100%; }
  #umap { height: 58vh; }
  #violin { height: 28vh; }
  aside { background:var(--card); border-left:1px solid #e7e0d4; padding:14px 16px; }
  aside h2 { font-size:12px; text-transform:uppercase; letter-spacing:.06em; color:var(--muted); margin:0 0 8px; }
  select, textarea, button, input { width:100%; font: inherit; margin: 0 0 8px; }
  textarea { min-height: 72px; resize: vertical; }
  button { background:var(--acc); color:#fff; border:0; padding:8px 10px; cursor:pointer; }
  button.secondary { background:#44403c; }
  #stats { font-variant-numeric: tabular-nums; white-space: pre-wrap; font-size:12px; color:var(--muted); }
  #reply { white-space: pre-wrap; font-size:13px; background:#f0fdf4; padding:8px; min-height:4em; }
  .warn { color:#b45309; font-size:12px; }
  .pubfigs { font-size:12px; max-height: 28vh; overflow:auto; }
  .pubfigs a { color: var(--acc); text-decoration:none; }
  .pubfigs li { margin: 0 0 6px 1.1em; }
</style>
</head>
<body>
<header>
  <h1>scAgent 交互查看器</h1>
  <span>框选 / 套索细胞后提问。静态 PNG 仍保留在 figures/。</span>
</header>
<div id="layout">
  <div id="plots">
    <div id="umap"></div>
    <div id="violin"></div>
  </div>
  <aside>
    <h2>着色</h2>
    <select id="color"></select>
    <h2>选区</h2>
    <div id="stats">未选择。用 Plotly 的 Box / Lasso 工具框细胞。</div>
    <h2>问 Agent</h2>
    <textarea id="q" placeholder="例如：分析我框选的这组细胞"></textarea>
    <button type="button" id="ask">提交选区 + 问题</button>
    <button type="button" id="dl" class="secondary">下载 selection.json</button>
    <h2>发表级主图</h2>
    <div id="pubfigs" class="pubfigs"></div>
    <p class="warn" id="hint"></p>
    <div id="reply"></div>
  </aside>
</div>
<script>
const DATA = __PAYLOAD__;
const ASK = __ASK__;
function uniq(arr){ return [...new Set(arr.map(String))]; }
function colors(col){
  const vals = DATA.obs[col] || DATA.obs[Object.keys(DATA.obs)[0]] || [];
  return vals;
}
function paint(col){
  const c = colors(col);
  const tr = { x: DATA.x, y: DATA.y, mode:'markers', type:'scattergl',
    text: DATA.ids, customdata: DATA.ids,
    marker:{ size:5, opacity:0.75, color: c, colorscale:'Portland' },
    hovertemplate: '%{text}<extra></extra>' };
  Plotly.react('umap', [tr], {
    margin:{t:24,l:40,r:10,b:40},
    dragmode:'lasso',
    xaxis:{title:'UMAP1'}, yaxis:{title:'UMAP2'},
    title: col + ' · n=' + DATA.ids.length + (DATA.sampled ? ' (抽样)' : '')
  }, {responsive:true});
  const qcKey = Object.keys(DATA.qc || {})[0];
  if (qcKey){
    Plotly.react('violin', [{ type:'violin', y: DATA.qc[qcKey], x: c, points:false, line:{color:'#0f766e'} }],
      { margin:{t:8,l:40,r:10,b:40}, title: qcKey + ' by ' + col, xaxis:{title:col} }, {responsive:true});
  }
}
const colorSel = document.getElementById('color');
(DATA.color_by || Object.keys(DATA.obs)).forEach(k => {
  const o = document.createElement('option'); o.value=k; o.textContent=k; colorSel.appendChild(o);
});
colorSel.onchange = () => paint(colorSel.value);
paint(colorSel.value || 'leiden');

(function renderPubFigs(){
  const box = document.getElementById('pubfigs');
  const items = DATA.publication_figures || [];
  if (!items.length){ box.textContent = '（执行后见 report.md 发表级清单）'; return; }
  const ul = document.createElement('ul');
  items.forEach(it => {
    const li = document.createElement('li');
    const title = it.title_zh || it.id || 'figure';
    const st = it.status || '';
    let txt = title + ' [' + st + ']';
    if (it.path && it.status === 'present') {
      const parts = String(it.path).split(';').map(s => s.trim()).filter(Boolean);
      txt += ' ';
      parts.forEach((p, i) => {
        const name = p.split('/').pop();
        const a = document.createElement('a');
        a.href = p.startsWith('figures/') ? p : ('figures/' + name);
        a.textContent = name;
        a.target = '_blank';
        li.appendChild(document.createTextNode(i ? ' · ' : ''));
        li.appendChild(a);
      });
    } else {
      li.textContent = txt;
    }
    if (!li.textContent) li.textContent = txt;
    ul.appendChild(li);
  });
  box.appendChild(ul);
})();

let selected = [];
document.getElementById('umap').on('plotly_selected', ev => {
  selected = (ev && ev.points || []).map(p => p.customdata || p.text).filter(Boolean);
  const n = selected.length;
  document.getElementById('stats').textContent = n ? ('已选 ' + n + ' 细胞\\n' + selected.slice(0,8).join('\\n') + (n>8?'\\n…':'')) : '未选择。';
});
function payload(){
  return { cell_ids: selected, n: selected.length, query: document.getElementById('q').value || '分析我框选的这组细胞' };
}
document.getElementById('dl').onclick = () => {
  if (!selected.length){ document.getElementById('hint').textContent='先框选细胞'; return; }
  const blob = new Blob([JSON.stringify(payload(), null, 2)], {type:'application/json'});
  const a = document.createElement('a'); a.href = URL.createObjectURL(blob); a.download='selection.json'; a.click();
  document.getElementById('hint').textContent = '然后: python -m scagent ask --selection selection.json';
};
document.getElementById('ask').onclick = async () => {
  if (!selected.length){ document.getElementById('hint').textContent='先框选细胞'; return; }
  const body = payload();
  document.getElementById('reply').textContent = '提交中…';
  if (ASK){
    try {
      const r = await fetch(ASK, { method:'POST', headers:{'Content-Type':'application/json'}, body: JSON.stringify(body) });
      const j = await r.json();
      document.getElementById('reply').textContent = j.summary || JSON.stringify(j);
      document.getElementById('hint').textContent = j.path ? ('已写入 ' + j.path) : '';
    } catch (e) {
      document.getElementById('reply').textContent = '本地服务不可用。已改为下载 JSON。';
      document.getElementById('dl').click();
    }
  } else {
    document.getElementById('dl').click();
    document.getElementById('reply').textContent = '无本地服务。下载 selection.json 后运行:\\npython -m scagent ask --selection selection.json';
  }
};
</script>
</body>
</html>
"""


def build_payload(adata, *, max_cells: int = MAX_CELLS, seed: int = 0) -> dict[str, Any]:
    import numpy as np

    n = int(adata.n_obs)
    idx = np.arange(n)
    sampled = n > max_cells
    if sampled:
        rng = np.random.default_rng(seed)
        idx = np.sort(rng.choice(n, size=max_cells, replace=False))
    names = np.asarray(adata.obs_names)[idx].astype(str).tolist()
    if "X_umap" in adata.obsm:
        xy = np.asarray(adata.obsm["X_umap"])[idx]
    elif "X_pca" in adata.obsm:
        xy = np.asarray(adata.obsm["X_pca"])[idx, :2]
    else:
        xy = np.zeros((len(idx), 2))
    obs: dict[str, list] = {}
    color_by: list[str] = []
    for col in _COLOR_COLS:
        if col in adata.obs:
            obs[col] = [str(v) for v in adata.obs[col].to_numpy()[idx]]
            color_by.append(col)
    qc: dict[str, list] = {}
    for col in _QC_COLS:
        if col in adata.obs:
            qc[col] = [float(v) if v == v else None for v in adata.obs[col].to_numpy()[idx]]
    return {
        "n_total": n,
        "n_shown": len(names),
        "sampled": sampled,
        "ids": names,
        "x": xy[:, 0].astype(float).tolist(),
        "y": xy[:, 1].astype(float).tolist(),
        "obs": obs,
        "qc": qc,
        "color_by": color_by or list(obs),
        "embedding": "X_umap" if "X_umap" in adata.obsm else ("X_pca" if "X_pca" in adata.obsm else "none"),
    }


def write_viewer_html(payload: dict[str, Any], dest: Path, *, ask_endpoint: str | None = None) -> Path:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    blob = json.dumps(payload, ensure_ascii=False)
    blob = blob.replace("<", "\\u003c")
    html = _HTML.replace("__PAYLOAD__", blob).replace("__ASK__", json.dumps(ask_endpoint))
    dest.write_text(html, encoding="utf-8")
    log.info("viewer %s cells=%s", dest, payload.get("n_shown"))
    return dest


def load_selection(path: str | Path) -> dict[str, Any]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    ids = data.get("cell_ids") or data.get("barcodes") or data.get("ids") or []
    ids = [str(x) for x in ids]
    return {"cell_ids": ids, "n": len(ids), "query": str(data.get("query") or "")}


def summarize_selection(adata, cell_ids: list[str], *, query: str = "") -> dict[str, Any]:
    import numpy as np

    want = set(cell_ids)
    names = np.asarray(adata.obs_names).astype(str)
    mask = np.fromiter((n in want for n in names), dtype=bool, count=len(names))
    n_hit = int(mask.sum())
    obs = adata.obs
    comp: dict[str, dict[str, int]] = {}
    for col in ("leiden", "cell_type", "cell_type_l1", "sample"):
        if col not in obs:
            continue
        vals = np.asarray(obs[col].astype(str))[mask]
        comp[col] = dict(Counter(vals.tolist()))
    qc_mean: dict[str, float] = {}
    for col in _QC_COLS:
        if col in obs:
            v = np.asarray(obs[col], dtype=float)[mask]
            v = v[np.isfinite(v)]
            if len(v):
                qc_mean[col] = float(np.mean(v))
    markers: dict[str, dict[str, float]] = {}
    for gene in ("CD3D", "MS4A1", "CD14", "NKG7", "IL7R"):
        if gene not in adata.var_names:
            continue
        x = adata[:, gene].X
        if hasattr(x, "toarray"):
            x = x.toarray()
        x = np.asarray(x).ravel()
        ins = float(np.mean(x[mask])) if n_hit else 0.0
        out = float(np.mean(x[~mask])) if int((~mask).sum()) else 0.0
        markers[gene] = {"in": round(ins, 4), "rest": round(out, 4)}
    summary = {
        "query": query or "分析我框选的这组细胞",
        "n_requested": len(cell_ids),
        "n_matched": n_hit,
        "frac": round(n_hit / max(adata.n_obs, 1), 4),
        "composition": comp,
        "qc_mean": qc_mean,
        "marker_means_in_vs_rest": markers,
    }
    lines = [
        f"选区 n={n_hit}/{adata.n_obs} ({summary['frac']:.1%})",
        f"问题: {summary['query']}",
    ]
    for col, counts in comp.items():
        top = ", ".join(f"{k}:{v}" for k, v in list(counts.items())[:8])
        lines.append(f"{col}: {top}")
    if qc_mean:
        lines.append("QC 均值: " + ", ".join(f"{k}={v:.3g}" for k, v in qc_mean.items()))
    if markers:
        lines.append("marker in vs rest: " + ", ".join(f"{g} {d['in']:.3g}/{d['rest']:.3g}" for g, d in markers.items()))
    summary["text"] = "\n".join(lines)
    return summary


def find_workspace_h5ad(workspace: Path | None = None) -> Path | None:
    cfg = load_config()
    ws = Path(workspace) if workspace else resolve_path(cfg, "workspace")
    for name in ("adata_processed.h5ad", "adata_qc.h5ad"):
        p = ws / name
        if p.is_file():
            return p
    return None


def export_workspace_viewer(
    out_dir: Path,
    *,
    h5ad: Path | None = None,
    workspace: Path | None = None,
    ask_endpoint: str | None = None,
    state: dict | None = None,
) -> Path | None:
    path = Path(h5ad) if h5ad else find_workspace_h5ad(workspace)
    if path is None or not path.is_file():
        log.info("viewer skipped: no h5ad")
        return None
    from scagent.io import read_h5ad

    adata = read_h5ad(path, backed=True)
    payload = build_payload(adata)
    if state:
        from scagent.publication_figures import build_publication_figure_inventory

        inv = build_publication_figure_inventory(state)
        payload["publication_figures"] = inv.get("items") or []
    dest = Path(out_dir) / "viewer.html"
    write_viewer_html(payload, dest, ask_endpoint=ask_endpoint)
    (Path(out_dir) / "viewer_payload.json").write_text(json.dumps({"n_total": payload["n_total"], "n_shown": payload["n_shown"], "sampled": payload["sampled"], "embedding": payload["embedding"]}, indent=2), encoding="utf-8")
    return dest


def serve_viewer(out_dir: Path, *, port: int = 8765, h5ad: Path | None = None) -> None:
    """Local server: GET viewer, POST /ask writes selection.json and returns a summary."""
    import threading
    from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
    import webbrowser

    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    html_path = export_workspace_viewer(out_dir, h5ad=h5ad, ask_endpoint="/ask")
    if html_path is None:
        raise FileNotFoundError("no adata_processed.h5ad / adata_qc.h5ad to view")
    adata_path = find_workspace_h5ad() if h5ad is None else Path(h5ad)

    class Handler(BaseHTTPRequestHandler):
        def log_message(self, fmt, *args):
            log.info("%s", fmt % args)

        def _send(self, code: int, body: bytes, ctype: str) -> None:
            self.send_response(code)
            self.send_header("Content-Type", ctype)
            self.send_header("Content-Length", str(len(body)))
            self.end_headers()
            self.wfile.write(body)

        def do_GET(self):
            if self.path in {"/", "/viewer.html"}:
                self._send(200, html_path.read_bytes(), "text/html; charset=utf-8")
                return
            self._send(404, b"not found", "text/plain")

        def do_POST(self):
            if self.path != "/ask":
                self._send(404, b"not found", "text/plain")
                return
            n = int(self.headers.get("Content-Length") or 0)
            raw = self.rfile.read(n)
            sel = json.loads(raw.decode("utf-8") or "{}")
            sel_path = out_dir / "selection.json"
            sel_path.write_text(json.dumps(sel, ensure_ascii=False, indent=2), encoding="utf-8")
            from scagent.io import read_h5ad

            adata = read_h5ad(adata_path, backed=True)
            summary = summarize_selection(adata, sel.get("cell_ids") or [], query=str(sel.get("query") or ""))
            (out_dir / "selection_report.md").write_text(summary["text"] + "\n", encoding="utf-8")
            body = json.dumps({"ok": True, "path": str(sel_path), "summary": summary["text"], "n": summary["n_matched"]}, ensure_ascii=False).encode()
            self._send(200, body, "application/json")

    httpd = ThreadingHTTPServer(("127.0.0.1", port), Handler)
    log.info("viewer http://127.0.0.1:%s", port)
    threading.Thread(target=lambda: webbrowser.open(f"http://127.0.0.1:{port}/"), daemon=True).start()
    httpd.serve_forever()
