import json

import pytest

from scagent.cli import build_parser
from scagent.viewer import build_payload, load_selection, summarize_selection, write_viewer_html


def _tiny_adata():
    np = pytest.importorskip("numpy")
    pytest.importorskip("anndata")
    from anndata import AnnData

    rng = np.random.default_rng(0)
    adata = AnnData(rng.poisson(1.0, size=(40, 8)).astype(np.float32))
    adata.obs_names = [f"c{i}" for i in range(40)]
    adata.var_names = ["MT-ND1", "CD3D", "MS4A1"] + [f"G{i}" for i in range(5)]
    adata.obs["leiden"] = ["0" if i < 15 else "1" for i in range(40)]
    adata.obs["cell_type"] = ["T" if i < 15 else "B" for i in range(40)]
    adata.obs["pct_counts_mt"] = rng.uniform(1, 8, size=40)
    adata.obsm["X_umap"] = rng.normal(size=(40, 2)).astype(np.float32)
    return adata


def test_payload_and_html(tmp_path):
    adata = _tiny_adata()
    payload = build_payload(adata, max_cells=20)
    assert payload["sampled"] is True
    assert payload["n_shown"] == 20
    assert payload["n_total"] == 40
    assert "leiden" in payload["obs"]
    html = write_viewer_html(payload, tmp_path / "viewer.html", ask_endpoint="/ask")
    text = html.read_text(encoding="utf-8")
    assert "plotly" in text.lower()
    assert "lasso" in text
    assert "plotly_selected" in text
    assert "/ask" in text
    assert "问 Agent" in text


def test_summarize_lasso_subset():
    adata = _tiny_adata()
    ids = [f"c{i}" for i in range(10)]
    s = summarize_selection(adata, ids, query="分析我框选的这组细胞")
    assert s["n_matched"] == 10
    assert s["composition"]["leiden"]["0"] == 10
    assert "CD3D" in s["marker_means_in_vs_rest"]
    assert "分析我框选" in s["text"]


def test_load_selection(tmp_path):
    p = tmp_path / "selection.json"
    p.write_text(json.dumps({"cell_ids": ["c0", "c1"], "query": "这组是什么"}), encoding="utf-8")
    sel = load_selection(p)
    assert sel["n"] == 2
    assert sel["query"] == "这组是什么"


def test_cli_view_ask_selection():
    p = build_parser()
    v = p.parse_args(["view", "--serve", "--port", "9001"])
    assert v.func.__name__ == "cmd_view"
    assert v.serve is True
    a = p.parse_args(["ask", "分析我框选的这组细胞", "--selection", "sel.json"])
    assert a.func.__name__ == "cmd_ask"
    r = p.parse_args(["run", "继续", "--selection", "sel.json"])
    assert r.selection == "sel.json"


def test_planner_records_selection():
    from agents.planner import build_plan

    plan = build_plan(
        {
            "user_query": "分析我框选的这组细胞",
            "metadata": {"tissue": "pbmc", "species": "human", "n_samples": 1},
            "language": "python",
            "selection": {"cell_ids": ["c0", "c1"], "n": 2},
        }
    )
    assert plan.get("selection_n") == 2
    assert any("框选" in x for x in plan.get("risks") or [])
