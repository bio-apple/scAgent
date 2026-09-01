"""Annotation evidence-chain gates (dual expression, R-continue, agree)."""

from __future__ import annotations

from agents.markers import choose_celltypist_model
from agents.reviewer import audit_code
from agents.templates import cluster_annotate_script
from scagent.annotate import dual_validate_expression, labels_agree


def test_template_continues_after_r_ref():
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5},
        {"integrator": None},
    )
    assert "SCAGENT_R_REF_OK" in code
    assert "_R_REF" in code
    assert "dual_validate_expression" in code
    assert "apply_ontology_ids" in code
    assert "annotation_dual_rate" in code
    assert "raise SystemExit(0)" not in code


def test_audit_requires_expression_dual():
    code = (
        "import scanpy as sc\n"
        'positive = ["CD3D", "CD3E"]\n'
        'negative = ["MS4A1"]\n'
        "fuse_annotation(adata)\n"
        'adata.obs["cell_type"] = "T cell"\n'
    )
    r = audit_code(code, {"tissue": "pbmc"}, phase="downstream")
    assert any(x.get("id") == "down.dual_expression" for x in r.get("issue_records") or [])


def test_dual_validate_expression_thresholds():
    ok = dual_validate_expression([0.5, 0.4, 0.01], [0.05], pos_min=0.1, neg_max=0.5)
    assert ok["dual_ok"] is True
    bad = dual_validate_expression([0.5, 0.01], [0.8], pos_min=0.1, neg_max=0.5)
    assert bad["dual_ok"] is False


def test_mouse_skips_human_celltypist_and_catalog():
    assert choose_celltypist_model("lung", "mouse") is None
    assert choose_celltypist_model("pbmc", "mouse") is None
    assert choose_celltypist_model("brain", "mouse") == "Developing_Mouse_Brain.pkl"
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "mouse", "tissue": "lung"},
        {},
        {"integrator": None},
    )
    assert "human marker catalog disabled" in code or "MARKERS = []" in code


def test_labels_agree_no_parent_subtype_substring():
    assert labels_agree("T cell", "CD8 T") is False
    assert labels_agree("B cell", "B cells") is True
