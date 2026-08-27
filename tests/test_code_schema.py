from agents.code_schema import validate_script
from agents.reviewer import audit_code
from agents.templates import cluster_annotate_script, qc_preprocess_script
from sandbox.executor import write_and_maybe_run


def test_syntax_error_fails_schema():
    r = validate_script("def (\n", phase="qc")
    assert r["ok"] is False
    assert r["ast_ok"] is False
    assert any("语法" in x for x in r["issues"])


def test_de_before_pca_fails():
    code = (
        "import scanpy as sc\n"
        "sc.tl.rank_genes_groups(adata, 'leiden')\n"
        "sc.tl.pca(adata)\n"
        "sc.pp.neighbors(adata)\n"
        "sc.tl.umap(adata)\n"
        "sc.tl.leiden(adata)\n"
    )
    r = validate_script(code, phase="downstream")
    assert r["ok"] is False
    assert any(rec["id"] == "schema.dag_de" for rec in r["issue_records"])
    audit = audit_code(code, {"tissue": "pbmc"}, phase="downstream")
    assert audit["passed"] is False
    assert any(rec["id"] == "schema.dag_de" for rec in audit["issue_records"])


def test_dpt_before_clustering_fails():
    code = (
        "import scanpy as sc\n"
        "sc.tl.pca(adata)\n"
        "sc.pp.neighbors(adata)\n"
        "sc.tl.dpt(adata)\n"
        "sc.tl.umap(adata)\n"
        "sc.tl.leiden(adata)\n"
    )
    r = validate_script(code, phase="downstream")
    assert r["ok"] is False
    assert any(rec["id"] == "schema.dag_traj" for rec in r["issue_records"])


def test_seurat_findmarkers_before_pca_fails():
    code = "library(Seurat)\nFindMarkers(obj)\nRunPCA(obj)\nFindNeighbors(obj)\nFindClusters(obj)\nRunUMAP(obj)\n"
    r = validate_script(code, phase="downstream", language="r")
    assert r["ok"] is False
    assert any(rec["id"] == "schema.dag_de" for rec in r["issue_records"])


def test_valid_template_passes_schema():
    qc = qc_preprocess_script({"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"}, {"nmads": 5})
    assert validate_script(qc, phase="qc")["ok"] is True
    down = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5},
        {"integrator": None},
    )
    r = validate_script(down, phase="downstream")
    assert r["ok"] is True, r["issues"]
    assert "cluster_deg" in r["steps"]
    assert r["steps"].index("pca") < r["steps"].index("leiden")
    assert r["steps"].index("umap") < r["steps"].index("cluster_deg")


def test_trajectory_template_after_umap_leiden():
    down = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5},
        {"integrator": None, "route": ["qc", "pca", "neighbors", "umap", "leiden", "trajectory"]},
    )
    assert "run_trajectory_phase" in down
    assert "sc.tl.dpt" not in down
    r = validate_script(down, phase="downstream")
    assert r["ok"] is True, r["issues"]
    assert r["steps"].index("leiden") < r["steps"].index("trajectory")
    assert r["steps"].index("umap") < r["steps"].index("trajectory")
    audit = audit_code(down, {"tissue": "pbmc", "species": "human", "route": ["trajectory"]}, phase="downstream")
    assert audit["passed"] is True


def test_executor_blocks_schema_before_sandbox(tmp_path):
    r = write_and_maybe_run(
        "import scanpy as sc\nsc.tl.rank_genes_groups(adata, 'leiden')\n",
        workspace=tmp_path,
        execute=True,
        timeout=5,
        extra_manifest={"phase": "downstream"},
    )
    assert r["ok"] is False
    assert r["executed"] is False
    assert r["jail"] == "schema"
    assert r["returncode"] == 125


def test_qc_script_cannot_contain_dpt():
    code = "import scanpy as sc\nnp.random.seed(0)\nsc.tl.dpt(adata)\n"
    r = validate_script(code, phase="qc")
    assert r["ok"] is False
    assert any(rec["id"] == "schema.qc_order" for rec in r["issue_records"])
