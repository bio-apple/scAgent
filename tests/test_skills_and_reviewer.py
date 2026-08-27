from agents.planner import choose_integrator
from agents.reviewer import audit_code, audit_execution
from agents.templates import LOCKED_START, cluster_annotate_script, qc_preprocess_script, scanpy_script
from scagent.inspect_data import detect_platform, detect_species_from_genes, inspect_data
from scagent.skills_loader import list_skills, recommend_skills


def test_existing_skills_preserved():
    names = {s.name for s in list_skills()}
    expected = {
        "anndata-data-structure",
        "celltypist-cell-annotation",
        "cellxgene-census",
        "harmony-batch-correction",
        "popv-cell-annotation",
        "scanpy-scrna-seq",
        "scvi-tools-single-cell",
        "single-cell-annotation-guide",
    }
    assert expected <= names


def test_recommend_harmony_for_multi_sample():
    skills = recommend_skills({"n_samples": 4, "need_batch_correction": True, "tissue": "pbmc"})
    assert "scanpy-scrna-seq" in skills
    assert "harmony-batch-correction" in skills
    assert "celltypist-cell-annotation" in skills


def test_choose_integrator():
    assert choose_integrator({"n_samples": 1, "need_batch_correction": False}) is None
    assert choose_integrator({"n_samples": 3, "n_cells": 8000, "need_batch_correction": True}) == "harmony"
    assert choose_integrator({"n_samples": 10, "n_cells": 1000, "need_batch_correction": True}) == "scvi"
    assert choose_integrator({"n_samples": 2, "n_cells": 200_000, "need_batch_correction": True}) == "scvi"


def test_inspect_parse_platform(tmp_path):
    p = tmp_path / "parse_biosciences_sample.h5ad"
    p.write_text("not-an-h5ad")
    meta = inspect_data(str(p), tissue="PBMC")
    assert meta["platform"] == "parse"
    assert meta["tissue"] == "pbmc"


def test_detect_species_and_platform():
    assert detect_species_from_genes(["MT-ND1", "GAPDH", "MT-CO1"]) == "human"
    assert detect_species_from_genes(["mt-Nd1", "Gapdh"]) == "mouse"
    assert detect_platform(__import__("pathlib").Path("filtered_feature_bc_matrix.h5")) == "10x"


def test_audit_requires_qc_trio():
    bad = "import scanpy as sc\nsc.pp.filter_cells(adata, min_genes=200)\n"
    r = audit_code(bad, {"need_batch_correction": True}, phase="qc")
    assert r["passed"] is False
    assert r["has_violin"] is False
    assert r["has_mad"] is False


def test_qc_template_passes_reviewer():
    meta = {"data_path": "x.h5ad", "species": "human", "need_batch_correction": True, "sample_key": "batch", "tissue": "pbmc"}
    code = qc_preprocess_script(meta, {"nmads": 5, "extra_qc": ["hb"]})
    compile(code, "<qc>", "exec")
    assert LOCKED_START in code
    assert "log1p=True" in code
    assert 'side="high"' in code
    assert "scrublet" in code
    r = audit_code(code, meta, phase="qc")
    assert r["has_violin"] and r["has_scatter"] and r["has_mad"] and r["has_locked"]
    assert r["passed"] is True


def test_downstream_template_dual_validation():
    meta = {"data_path": "x.h5ad", "species": "human", "need_batch_correction": True, "sample_key": "batch", "tissue": "pbmc"}
    plan = {"integrator": "harmony"}
    code = cluster_annotate_script(meta, {"nmads": 5}, plan)
    compile(code, "<down>", "exec")
    assert "celltypist" in code.lower()
    assert "positive" in code and "negative" in code
    assert "harmonypy" in code
    r = audit_code(code, meta, phase="downstream")
    assert r["has_celltypist"] and r["has_dual"]
    assert r["passed"] is True


def test_combined_script_compiles():
    code = scanpy_script(
        {"data_path": "x.h5ad", "species": "human", "need_batch_correction": True, "sample_key": "batch"},
        {"nmads": 5},
        {"integrator": "harmony"},
    )
    compile(code, "<scanpy_script>", "exec")
    assert "harmony" in code.lower()


def test_audit_execution_fail_and_skip():
    skip = audit_execution({}, {}, phase="qc", execute_code=False)
    assert skip["passed"] is True and skip["skipped"] is True
    fail = audit_execution(
        {"executed": True, "ok": False, "stderr": "Traceback"},
        {"metrics": {}, "h5ads": {}, "figures": []},
        phase="qc",
        execute_code=True,
    )
    assert fail["passed"] is False
    over = audit_execution(
        {"executed": True, "ok": True},
        {"metrics": {"pct_removed": 55}, "h5ads": {"qc": "/tmp/adata_qc.h5ad"}, "figures": ["/tmp/violin.png", "/tmp/scatter.png"]},
        phase="qc",
        execute_code=True,
    )
    assert over["passed"] is False
    assert any("过度过滤" in x for x in over["issues"])
