import numpy as np
import pytest

from agents.code_schema import validate_script
from agents.templates import cluster_annotate_script
from scagent.annotate import ensemble_cell_annotation
from scagent.deg_methods import force_pseudobulk_de, resolve_forced_deg_engine


def test_force_pseudobulk_de_requires_replicates():
    assert force_pseudobulk_de({"condition_key": "group", "n_replicates": 2}) is True
    assert force_pseudobulk_de({"condition_key": "group", "n_replicates": 1}) is False
    assert force_pseudobulk_de({"n_replicates": 3}) is False
    assert force_pseudobulk_de({}, {"condition_key": "treat", "force_pseudobulk_de": True}) is True


def test_resolve_forced_deg_engine_blocks_ttest():
    assert resolve_forced_deg_engine("ttest") == "auto"
    assert resolve_forced_deg_engine("edger") == "edger"


def test_ensemble_celltypist_only_when_all_high_conf():
    pytest.importorskip("anndata")
    import anndata as ad

    adata = ad.AnnData(X=np.random.poisson(2, (20, 10)).astype(float))
    adata.obs["celltypist_label"] = ["T"] * 10 + ["B"] * 10
    adata.obs["celltypist_conf"] = 0.95
    out = ensemble_cell_annotation(adata, sample_key=None)
    assert (out.obs["scagent_annotation"] == out.obs["celltypist_label"]).all()
    assert out.uns["scagent_annotation"]["method"] == "celltypist"
    assert out.uns["scagent_annotation"]["scanvi_ran"] is False


def test_schema_rejects_condition_wilcoxon_when_forced():
    bad = """
import scanpy as sc
sc.pp.pca(adata)
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)
sc.tl.rank_genes_groups(adata, groupby='condition', method='wilcoxon')
"""
    meta = {"condition_key": "condition", "n_replicates": 3}
    plan = {"force_pseudobulk_de": True, "condition_key": "condition"}
    rep = validate_script(bad, phase="downstream", metadata=meta, plan=plan)
    assert rep["ok"] is False
    assert any("pseudobulk" in i for i in rep["issues"])


def test_template_includes_ensemble_and_forced_pseudobulk():
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "lung", "condition_key": "group", "n_replicates": 3},
        {"nmads": 5},
        {
            "integrator": None,
            "needs_pseudobulk": True,
            "force_pseudobulk_de": True,
            "condition_key": "group",
            "deg_engine": "auto",
        },
    )
    compile(code, "<dn>", "exec")
    assert "ensemble_cell_annotation" in code
    assert "scagent_annotation" in code
    assert "FORCE_PSEUDOBULK_DE" in code
    assert "pseudobulk_de" in code
