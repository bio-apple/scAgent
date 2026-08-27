from pathlib import Path

import numpy as np
import pytest

from agents.templates import cluster_annotate_script
from scagent.analysis import _read_deg_csv, pseudobulk_de


def _condition_adata(n_genes: int = 40):
    pytest.importorskip("anndata")
    from anndata import AnnData
    from scipy import sparse

    rng = np.random.default_rng(0)
    samples = ["s1", "s2", "s3", "s4"]
    conds = ["ctrl", "ctrl", "treat", "treat"]
    n_per = 12
    n_obs = n_per * 4
    X = rng.poisson(2.0, size=(n_obs, n_genes)).astype(np.float32)
    X[2 * n_per :, 0] += 8  # gene0 up in treat samples
    adata = AnnData(sparse.csr_matrix(X))
    adata.obs_names = [f"c{i}" for i in range(n_obs)]
    adata.var_names = [f"G{i}" for i in range(n_genes)]
    adata.obs["sample"] = np.repeat(samples, n_per)
    adata.obs["condition"] = np.repeat(conds, n_per)
    adata.obs["cell_type"] = "T"
    adata.layers["counts"] = adata.X.copy()
    return adata


def test_pseudobulk_ttest_engine(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    adata = _condition_adata()
    out = pseudobulk_de(
        adata,
        sample_key="sample",
        condition_key="condition",
        groupby="cell_type",
        min_cells=5,
        min_replicates=2,
        engine="ttest",
        cross_validate=False,
    )
    info = out.uns["pseudobulk_de"]
    assert info["engine"] == "ttest_bh"
    assert info["ran"] is True
    assert Path("pseudobulk_de.csv").is_file()
    df_txt = Path("pseudobulk_de.csv").read_text(encoding="utf-8")
    assert "G0" in df_txt


def test_auto_uses_edger_when_rpy2_writes_csv(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fake_rpy2(cpath, mpath, opath, engine):
        Path(opath).write_text("gene,logFC,pval,fdr\nG0,1.5,0.001,0.01\n", encoding="utf-8")
        return "edger"

    monkeypatch.setattr("scagent.analysis._deg_via_rpy2", fake_rpy2)
    monkeypatch.setattr("scagent.analysis._deg_via_rscript", lambda *a, **k: None)
    adata = _condition_adata()
    out = pseudobulk_de(
        adata,
        sample_key="sample",
        condition_key="condition",
        groupby="cell_type",
        min_cells=5,
        engine="auto",
        cross_validate=False,
    )
    info = out.uns["pseudobulk_de"]
    assert info["engine"] == "edger_rpy2"
    assert info["n_sig"] == 1


def test_deseq2_rscript_fallback(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fake_script(cpath, mpath, opath, engine):
        assert engine in {"auto", "deseq2", "edger"}
        Path(opath).write_text("gene,logFC,pval,fdr\nG1,-0.8,0.02,0.04\n", encoding="utf-8")
        return "deseq2"

    monkeypatch.setattr("scagent.analysis._deg_via_rpy2", lambda *a, **k: None)
    monkeypatch.setattr("scagent.analysis._deg_via_rscript", fake_script)
    adata = _condition_adata()
    out = pseudobulk_de(
        adata,
        sample_key="sample",
        condition_key="condition",
        min_cells=5,
        engine="deseq2",
        cross_validate=False,
    )
    assert out.uns["pseudobulk_de"]["engine"] == "deseq2_rscript"


def test_read_deg_csv_skips_na(tmp_path):
    p = tmp_path / "out.csv"
    p.write_text("gene,logFC,pval,fdr\nA,1,0.1,0.2\nB,NA,NA,NA\n", encoding="utf-8")
    rows = _read_deg_csv(p, "T", "ctrl", "treat")
    assert [r["gene"] for r in rows] == ["A"]
    assert rows[0]["group_a"] == "ctrl"


def test_template_passes_deg_engine():
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "tissue": "pbmc", "condition_key": "condition"},
        {},
        {"integrator": None, "needs_pseudobulk": True, "condition_key": "condition", "deg_engine": "edger"},
    )
    assert "pseudobulk_de" in code
    assert "engine='edger'" in code or 'engine="edger"' in code
    compile(code, "<deg>", "exec")


def test_parse_deg_preference_from_query():
    from scagent.deg_methods import gene_overlap, parse_deg_preference

    p = parse_deg_preference("用 DESeq2 比较对照组和处理组")
    assert p["engine"] == "deseq2"
    p2 = parse_deg_preference("cluster marker 用 t-test，并与 Wilcoxon 交叉验证")
    assert p2["marker_method"] == "t-test"
    assert p2["cross_validate"] is True
    p3 = parse_deg_preference("只用 Wilcoxon 看 marker")
    assert p3["marker_method"] == "wilcoxon"
    assert p3["cross_validate"] is False
    p4 = parse_deg_preference("MAST")
    assert p4["marker_method"] == "mast"
    ov = gene_overlap(["A", "B", "C"], ["B", "C", "D"])
    assert ov["n_overlap"] == 2
    assert ov["jaccard"] == 0.5


def test_template_honors_marker_method_and_cross_validate():
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "tissue": "pbmc", "condition_key": "condition"},
        {},
        {
            "integrator": None,
            "needs_pseudobulk": True,
            "condition_key": "condition",
            "deg_engine": "deseq2",
            "marker_method": "t-test",
            "deg_cross_validate": True,
        },
    )
    assert "rank_genes" in code
    assert "t-test" in code
    assert "engine='deseq2'" in code or 'engine="deseq2"' in code
    assert "cross_validate=True" in code
    compile(code, "<deg_cv>", "exec")


def test_pseudobulk_cross_validate_overlap(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)

    def fake_rpy2(cpath, mpath, opath, engine):
        Path(opath).write_text("gene,logFC,pval,fdr\nG0,1.5,0.001,0.01\nG1,0.2,0.2,0.4\n", encoding="utf-8")
        return str(engine if engine != "auto" else "edger")

    monkeypatch.setattr("scagent.analysis._deg_via_rpy2", fake_rpy2)
    monkeypatch.setattr("scagent.analysis._deg_via_rscript", lambda *a, **k: None)
    adata = _condition_adata()
    out = pseudobulk_de(
        adata,
        sample_key="sample",
        condition_key="condition",
        groupby="cell_type",
        min_cells=5,
        engine="edger",
        cross_validate=True,
    )
    info = out.uns["pseudobulk_de"]
    assert info["cross_validate"] is True
    assert info["n_overlap"] == 1
    assert Path("pseudobulk_de_overlap.json").is_file()


def test_publication_review_reports_deg_overlap():
    from agents.reviewer import publication_review

    card = publication_review(
        {
            "metadata": {"tissue": "pbmc", "n_samples": 1, "need_batch_correction": False},
            "plan": {"needs_pseudobulk": True, "integrator": None},
            "code_downstream": "from scagent.analysis import rank_genes, pseudobulk_de\nrank_genes(adata)\npseudobulk_de(adata)\nprint('pvals_adj and pseudobulk + FDR')\n",
            "review_qc": {"passed": True},
            "review_downstream": {"passed": True, "has_dual": True},
            "execute_code": False,
            "artifacts": {
                "metrics": {
                    "deg_engine": "edger_rpy2",
                    "deg_n_overlap": 12,
                    "deg_jaccard": 0.4,
                    "deg_engines": ["edger", "deseq2"],
                }
            },
        }
    )
    deg = next(i for i in card["items"] if i["key"] == "deg")
    assert deg["status"] == "pass"
    assert "交叉验证" in deg["detail"]
    assert "12" in deg["detail"]

