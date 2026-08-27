from pathlib import Path

import pytest

from agents.artifacts import collect_workspace
from agents.reviewer import audit_execution, publication_review
from agents.writer import render_html, render_report, stage_report_figures


def _tiny_batch_adata():
    np = pytest.importorskip("numpy")
    ad = pytest.importorskip("anndata")
    n = 40
    rng = np.random.default_rng(0)
    adata = ad.AnnData(rng.poisson(2, size=(n, 16)).astype(float))
    adata.obs["sample"] = ["a"] * 20 + ["b"] * 20
    pca = rng.normal(size=(n, 8))
    pca[:20] += 3.0
    adata.obsm["X_pca"] = pca
    adata.obsm["X_pca_harmony"] = rng.normal(size=(n, 8))
    adata.obsm["X_umap"] = rng.normal(size=(n, 2))
    return adata


def test_integration_diagnostics_writes_png(tmp_path):
    pytest.importorskip("matplotlib")
    from scagent.plotting import integration_diagnostics

    plots = integration_diagnostics(_tiny_batch_adata(), "sample", figdir=tmp_path)
    names = {Path(p).name for p in plots}
    assert "batch_pca_before.png" in names
    assert "batch_pca_after.png" in names
    assert "batch_umap_after.png" in names
    for p in plots:
        assert Path(p).is_file() and Path(p).stat().st_size > 0


def test_caption_kind_prefers_batch_umap_over_umap(tmp_path):
    fig = tmp_path / "figures"
    fig.mkdir()
    (fig / "batch_umap_before.png").write_bytes(b"png")
    (fig / "umap_overview.png").write_bytes(b"png")
    arts = collect_workspace(tmp_path, "downstream", {"executed": True, "ok": True})
    kinds = {c["kind"] for c in arts["figure_captions"]}
    assert "batch_umap_before" in kinds
    assert "umap" in kinds
    by_name = {Path(c["path"]).name: c["kind"] for c in arts["figure_captions"]}
    assert by_name["batch_umap_before.png"] == "batch_umap_before"
    assert by_name["umap_overview.png"] == "umap"


def test_report_embeds_staged_batch_figures(tmp_path):
    src = tmp_path / "ws" / "figures"
    src.mkdir(parents=True)
    png = src / "batch_pca_before.png"
    png.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 32)
    after = src / "batch_umap_after.png"
    after.write_bytes(b"\x89PNG\r\n\x1a\n" + b"\x00" * 32)
    arts = {
        "figure_captions": [
            {"path": str(png), "kind": "batch_pca_before", "caption": "校正前 PCA"},
            {"path": str(after), "kind": "batch_umap_after", "caption": "校正后 UMAP"},
        ],
        "figures": [str(png), str(after)],
        "metrics": {"ilisi": 0.91, "kbet": 0.7},
    }
    out = tmp_path / "outputs"
    staged = stage_report_figures(arts, out)
    assert staged["figure_captions"][0]["path"] == "figures/batch_pca_before.png"
    assert (out / "figures" / "batch_pca_before.png").is_file()
    state = {
        "user_query": "demo",
        "report_lang": "zh",
        "metadata": {"n_samples": 2, "need_batch_correction": True, "species": "human", "tissue": "pbmc"},
        "plan": {"integrator": "harmony", "skills": [], "route": ["qc"]},
        "artifacts": staged,
    }
    md = render_report(state)
    assert "![校正前 PCA](figures/batch_pca_before.png)" in md
    assert "![校正后 UMAP](figures/batch_umap_after.png)" in md
    html = render_html(state)
    assert '<img src="figures/batch_pca_before.png"' in html
    assert "整合诊断" in html


def test_reviewer_requires_batch_plots_only_when_executed():
    meta = {"need_batch_correction": True, "n_samples": 2}
    skip = audit_execution({}, {}, phase="downstream", execute_code=False, metadata=meta)
    assert skip["skipped"] is True
    assert not any(r["id"] == "exec.integ_plot" for r in skip["issue_records"])

    miss = audit_execution(
        {"executed": True, "ok": True},
        {
            "h5ads": {"processed": "/tmp/adata_processed.h5ad"},
            "figures": ["/tmp/umap.png"],
            "metrics": {"ilisi": 0.9},
        },
        phase="downstream",
        execute_code=True,
        metadata=meta,
    )
    assert miss["passed"] is False
    assert any(r["id"] == "exec.integ_plot" for r in miss["issue_records"])

    ok = audit_execution(
        {"executed": True, "ok": True},
        {
            "h5ads": {"processed": "/tmp/adata_processed.h5ad"},
            "figures": ["/tmp/batch_pca_before.png", "/tmp/batch_umap_after.png"],
            "metrics": {"ilisi": 0.9},
        },
        phase="downstream",
        execute_code=True,
        metadata=meta,
    )
    assert ok["passed"] is True

    single = audit_execution(
        {"executed": True, "ok": True},
        {"h5ads": {"processed": "/tmp/adata_processed.h5ad"}, "figures": ["/tmp/umap.png"], "metrics": {}},
        phase="downstream",
        execute_code=True,
        metadata={"need_batch_correction": False, "n_samples": 1},
    )
    assert single["passed"] is True


def test_publication_review_embeds_batch_diag_when_executed():
    code = "sc.external.pp.harmony_integrate(adata, 'batch')\nsc.tl.umap(adata)\n"
    base = {
        "metadata": {"need_batch_correction": True, "n_samples": 2},
        "plan": {"integrator": "harmony"},
        "code_downstream": code,
        "execute_code": True,
        "execution_downstream": {"executed": True, "ok": True},
        "review_qc": {"passed": True, "issues": []},
        "review_downstream": {"passed": True, "issues": [], "has_dual": True, "has_fusion": True, "has_celltypist": True},
        "artifacts": {"metrics": {"ilisi": 0.9}, "figures": ["workspace/figures/umap.png"]},
    }
    card = publication_review(base)
    by_key = {i["key"]: i for i in card["items"]}
    assert by_key["batch_correction"]["status"] == "fail"
    assert "诊断图" in by_key["batch_correction"]["detail"]
    assert by_key["figures"]["status"] == "fail"

    ok_state = {
        **base,
        "artifacts": {
            "metrics": {"ilisi": 0.9},
            "figures": [
                "workspace/figures/batch_pca_before.png",
                "workspace/figures/batch_umap_after.png",
                "workspace/figures/umap.png",
            ],
        },
    }
    card2 = publication_review(ok_state)
    by_key2 = {i["key"]: i for i in card2["items"]}
    assert by_key2["batch_correction"]["status"] == "pass"
    assert "已嵌入" in by_key2["batch_correction"]["detail"]
