from pathlib import Path

import pytest

from agents.artifacts import collect_workspace
from agents.reviewer import publication_review
from agents.writer import render_report
from scagent.publication_figures import build_publication_figure_inventory, render_publication_figure_inventory_markdown


def _base_state(**kw):
    state = {
        "execute_code": True,
        "execution_qc": {"executed": True},
        "execution_downstream": {"executed": True},
        "execution_interpret": {"executed": True},
        "metadata": {"need_batch_correction": True, "n_samples": 2, "tissue": "pbmc"},
        "plan": {"needs_pseudobulk": True, "route": ["qc", "clustering", "enrichment"]},
        "artifacts": {
            "figure_captions": [
                {"path": "figures/qc_violin.png", "kind": "violin", "caption": "qc"},
                {"path": "figures/batch_umap_before.png", "kind": "batch_umap_before", "caption": "before"},
                {"path": "figures/batch_umap_after.png", "kind": "batch_umap_after", "caption": "after"},
                {"path": "figures/marker_heatmap.png", "kind": "marker_heatmap", "caption": "heatmap"},
                {"path": "figures/volcano.png", "kind": "volcano", "caption": "volcano"},
                {"path": "figures/pathway_bubble.png", "kind": "pathway_bubble", "caption": "bubble"},
            ],
            "figures": [],
            "metrics": {"ilisi": 0.9},
        },
        "review_qc": {"passed": True, "issues": []},
        "review_downstream": {"passed": True, "issues": [], "has_dual": True, "has_fusion": True, "has_celltypist": True},
        "code_qc": "violin",
        "code_downstream": "harmony pseudobulk_de rank_genes",
    }
    state.update(kw)
    return state


def test_inventory_complete_when_all_present():
    inv = build_publication_figure_inventory(_base_state())
    assert inv["n_required"] == 5
    assert inv["n_present"] == 5
    assert inv["complete"] is True
    assert inv["missing_ids"] == []


def test_inventory_flags_missing_volcano():
    arts = dict(_base_state()["artifacts"])
    arts["figure_captions"] = [c for c in arts["figure_captions"] if c["kind"] != "volcano"]
    inv = build_publication_figure_inventory(_base_state(artifacts=arts))
    assert "volcano" in inv["missing_ids"]
    assert inv["complete"] is False


def test_report_contains_checklist_section():
    md = render_report(_base_state(report_lang="zh"))
    assert "## 发表级图表清单" in md
    assert "QC violin" in md
    assert "pathway_bubble" in md or "通路气泡图" in md


def test_caption_kind_detects_new_figure_types(tmp_path):
    fig = tmp_path / "figures"
    fig.mkdir()
    (fig / "marker_heatmap.png").write_bytes(b"png")
    (fig / "volcano.png").write_bytes(b"png")
    (fig / "pathway_bubble.png").write_bytes(b"png")
    arts = collect_workspace(tmp_path, "downstream", {"executed": True, "ok": True})
    kinds = {c["kind"] for c in arts["figure_captions"]}
    assert "marker_heatmap" in kinds
    assert "volcano" in kinds
    assert "pathway_bubble" in kinds


def test_publication_review_uses_checklist():
    state = _base_state()
    arts = dict(state["artifacts"])
    arts["figure_captions"] = [c for c in arts["figure_captions"] if c["kind"] != "volcano"]
    state["artifacts"] = arts
    card = publication_review(state)
    fig_item = next(i for i in card["items"] if i["key"] == "figures")
    assert fig_item["status"] == "fail"
    assert "volcano" in fig_item["detail"]


def test_pathway_bubble_plot_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    from scagent.plotting import pathway_bubble_plot

    terms = [
        {"term": "HALLMARK_INTERFERON", "pval": 1e-4, "fdr": 0.01, "overlap": 5},
        {"term": "HALLMARK_HYPOXIA", "pval": 0.02, "fdr": 0.08, "overlap": 3},
    ]
    out = pathway_bubble_plot(terms, figdir=tmp_path)
    assert out is not None
    assert Path(out).is_file() and Path(out).stat().st_size > 0


def test_volcano_from_de_csv(tmp_path):
    pytest.importorskip("matplotlib")
    from scagent.plotting import volcano_from_de_csv

    csv = tmp_path / "pseudobulk_de.csv"
    csv.write_text("gene,logFC,pval,fdr\nG1,2.0,0.001,0.01\nG2,-1.0,0.2,0.3\n", encoding="utf-8")
    out = volcano_from_de_csv(csv, figdir=tmp_path)
    assert out is not None
    assert Path(out).name == "volcano.png"
