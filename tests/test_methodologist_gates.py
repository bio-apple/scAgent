"""Methodologist gates: marker fallback, pseudobulk force, annotate→pseudobulk, confound skills."""

from __future__ import annotations

from agents.code_schema import validate_script
from agents.markers import load_marker_catalog
from agents.reviewer import audit_code
from scagent.deg_methods import force_pseudobulk_de
from scagent.skills_loader import recommend_skills


def test_marker_catalog_refuses_pbmc_fallback_for_unknown_organs():
    lung = load_marker_catalog(tissue="lung")
    assert lung["cell_types"] and lung.get("warning") in (None, "")
    tumor = load_marker_catalog(tissue="tumor")
    assert tumor["cell_types"]
    skin = load_marker_catalog(tissue="skin")
    assert skin["cell_types"] == []
    assert "refusing PBMC" in (skin.get("warning") or "")
    unknown = load_marker_catalog(tissue="unknown")
    assert unknown["cell_types"] == []
    gut = load_marker_catalog(tissue="gut")
    assert gut["cell_types"]


def test_needs_pseudobulk_aligned_with_force_only():
    assert force_pseudobulk_de({"condition_key": "group", "n_replicates": 2}) is True
    assert force_pseudobulk_de({"condition_key": "group", "n_replicates": 1}) is False
    # Exploratory intent metadata alone must not trip hard gate
    code = "import scanpy as sc\nsc.tl.rank_genes_groups(adata, 'leiden')\n"
    audit = audit_code(
        "import scanpy as sc\n"
        "sc.pp.pca(adata)\nsc.pp.neighbors(adata)\nsc.tl.umap(adata)\nsc.tl.leiden(adata)\n"
        "sc.tl.rank_genes_groups(adata, groupby='leiden')\n",
        {"tissue": "pbmc", "needs_pseudobulk": True, "n_replicates": 1},
        phase="downstream",
    )
    # needs_pseudobulk without force + replicates must NOT emit down.pseudobulk hard fail
    assert not any(r.get("id") == "down.pseudobulk" for r in audit.get("issue_records") or [])


def test_schema_requires_annotate_before_pseudobulk():
    code = (
        "import scanpy as sc\n"
        "sc.pp.pca(adata)\n"
        "sc.pp.neighbors(adata)\n"
        "sc.tl.umap(adata)\n"
        "sc.tl.leiden(adata)\n"
        "pseudobulk_de(adata, sample_key='sample', condition_key='condition')\n"
    )
    r = validate_script(code, phase="downstream")
    assert r["ok"] is False
    assert any("annotate" in (iss or "").lower() or "pseudobulk" in (iss or "").lower() for iss in r["issues"])


def test_confounded_batch_still_gets_core_workflow():
    """Integration is gated at runtime when batch×condition confounded; core workflow skill stays."""
    skills = recommend_skills(
        {
            "n_samples": 4,
            "need_batch_correction": True,
            "batch_condition_confounded": True,
            "tissue": "pbmc",
        }
    )
    assert "seurat-workflow" in skills
    assert "cell-annotation" in skills
    skills2 = recommend_skills(
        {"n_samples": 4, "need_batch_correction": True, "tissue": "pbmc"}
    )
    assert "seurat-workflow" in skills2
