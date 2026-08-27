from agents.reviewer import audit_code
from agents.templates import scanpy_script
from scagent.inspect_data import inspect_data
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


def test_inspect_parse_platform(tmp_path):
    p = tmp_path / "parse_biosciences_sample.h5ad"
    p.write_text("not-an-h5ad")
    meta = inspect_data(str(p), tissue="PBMC")
    assert meta["platform"] == "parse"
    assert meta["tissue"] == "pbmc"


def test_audit_requires_qc_trio():
    bad = "import scanpy as sc\nsc.pp.filter_cells(adata, min_genes=200)\n"
    r = audit_code(bad, {"need_batch_correction": True})
    assert r["passed"] is False
    assert r["has_violin"] is False
    assert r["has_mad"] is False


def test_template_passes_reviewer():
    code = scanpy_script(
        {"data_path": "x.h5ad", "species": "human", "need_batch_correction": True, "sample_key": "batch"},
        {"nmads": 5},
    )
    compile(code, "<scanpy_script>", "exec")
    r = audit_code(code, {"need_batch_correction": True})
    assert r["has_violin"] and r["has_scatter"] and r["has_mad"]
    assert r["passed"] is True
    assert "harmony" in code.lower()
