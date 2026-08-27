import gzip
from pathlib import Path

import pytest

from agents.planner import choose_integrator, explain_integrator
from agents.templates import _load_block, cluster_annotate_script
from scagent.inspect_data import inspect_data
from scagent.io import (
    discover_samples,
    parse_data_spec,
    read_single_cell,
    resolve_10x_h5,
    resolve_10x_matrix_dir,
    sample_label,
)


def _write_10x_mtx(matrix_dir: Path, n_cells: int = 3, n_genes: int = 2) -> Path:
    matrix_dir.mkdir(parents=True, exist_ok=True)
    (matrix_dir / "matrix.mtx.gz").write_bytes(b"%%MatrixMarket\n")
    with gzip.open(matrix_dir / "barcodes.tsv.gz", "wt") as fh:
        for i in range(n_cells):
            fh.write(f"AAACCTG{i:04d}-1\n")
    with gzip.open(matrix_dir / "features.tsv.gz", "wt") as fh:
        for i in range(n_genes):
            fh.write(f"ENSG{i}\tGENE{i}\tGene Expression\n")
    return matrix_dir


def test_parse_data_spec_comma_and_glob(tmp_path):
    a = tmp_path / "a.h5ad"
    b = tmp_path / "b.h5ad"
    a.write_text("x")
    b.write_text("y")
    paths = parse_data_spec(f"{a},{b}")
    assert paths == [a, b]
    g = parse_data_spec(str(tmp_path / "*.h5ad"))
    assert set(g) == {a, b}


def test_resolve_cellranger_outs(tmp_path):
    mtx = _write_10x_mtx(tmp_path / "sampleA" / "outs" / "filtered_feature_bc_matrix")
    sample = tmp_path / "sampleA"
    outs = sample / "outs"
    assert resolve_10x_matrix_dir(sample) == mtx
    assert resolve_10x_matrix_dir(outs) == mtx
    assert resolve_10x_matrix_dir(mtx) == mtx
    assert sample_label(mtx) == "sampleA"
    h5 = outs / "filtered_feature_bc_matrix.h5"
    h5.write_bytes(b"not-h5")
    # mtx wins over h5 when both exist
    assert resolve_10x_matrix_dir(outs) == mtx
    # h5-only sample
    only = tmp_path / "sampleB" / "outs"
    only.mkdir(parents=True)
    cr_h5 = only / "filtered_feature_bc_matrix.h5"
    cr_h5.write_bytes(b"x")
    assert resolve_10x_h5(tmp_path / "sampleB") == cr_h5


def test_discover_two_cellranger_samples(tmp_path):
    root = tmp_path / "runs"
    s1 = _write_10x_mtx(root / "s1" / "outs" / "filtered_feature_bc_matrix", n_cells=2)
    s2 = _write_10x_mtx(root / "s2" / "outs" / "filtered_feature_bc_matrix", n_cells=4)
    found = discover_samples(root)
    assert set(found) == {s1, s2}
    meta = inspect_data(str(root), tissue="pbmc")
    assert meta["exists"] is True
    assert meta["n_samples"] == 2
    assert meta["need_batch_correction"] is True
    assert meta["sample_key"] == "sample"
    assert meta["n_cells"] == 6
    assert meta["platform"] in {"10x", "multi"}
    assert any("拼接" in n or "批次" in n for n in meta["notes"])
    assert choose_integrator(meta, "auto") == "harmony"


def test_inspect_cellranger_outs_counts(tmp_path):
    mtx = _write_10x_mtx(tmp_path / "s1" / "outs" / "filtered_feature_bc_matrix", n_cells=5, n_genes=4)
    meta = inspect_data(str(tmp_path / "s1"), tissue="pbmc")
    assert meta["platform"] == "10x"
    assert meta["n_cells"] == 5
    assert meta["n_genes"] == 4
    assert meta["n_samples"] == 1
    assert meta["need_batch_correction"] is False
    assert any("10x" in n or "Cell Ranger" in n for n in meta["notes"])
    assert resolve_10x_matrix_dir(tmp_path / "s1") == mtx


def test_inspect_h5ad_batch_triggers_auto(tmp_path):
    pytest.importorskip("anndata")
    pytest.importorskip("numpy")
    import numpy as np
    from anndata import AnnData

    rng = np.random.default_rng(0)
    adata = AnnData(rng.poisson(1.0, size=(30, 8)).astype(np.float32))
    adata.obs["batch"] = ["b1"] * 10 + ["b2"] * 10 + ["b3"] * 10
    path = tmp_path / "batched.h5ad"
    adata.write_h5ad(path)
    meta = inspect_data(str(path), tissue="pbmc")
    assert meta["sample_key"] == "batch"
    assert meta["n_samples"] == 3
    assert meta["need_batch_correction"] is True
    assert any("将触发批次校正" in n for n in meta["notes"])
    assert choose_integrator(meta, "auto") == "harmony"
    reason = explain_integrator(meta, "auto", "harmony")
    assert "Harmony" in reason
    assert "batch" in reason


def test_concat_two_h5ads(tmp_path):
    pytest.importorskip("anndata")
    pytest.importorskip("numpy")
    import numpy as np
    from anndata import AnnData

    rng = np.random.default_rng(1)

    def _one(name, n):
        a = AnnData(rng.poisson(1.0, size=(n, 6)).astype(np.float32))
        a.var_names = [f"g{i}" for i in range(6)]
        a.obs_names = [f"{name}_{i}" for i in range(n)]
        p = tmp_path / f"{name}.h5ad"
        a.write_h5ad(p)
        return p

    p1, p2 = _one("s1", 4), _one("s2", 5)
    meta = inspect_data(f"{p1},{p2}", tissue="pbmc")
    assert meta["n_samples"] == 2
    assert meta["need_batch_correction"] is True
    assert meta["n_cells"] == 9
    adata = read_single_cell(f"{p1},{p2}")
    assert "sample" in adata.obs
    assert set(adata.obs["sample"].astype(str)) == {"s1", "s2"}
    assert adata.n_obs == 9


def test_io_loom_suffix_is_dispatched(tmp_path):
    p = tmp_path / "x.loom"
    p.write_bytes(b"not-a-loom")
    with pytest.raises((ValueError, RuntimeError, OSError, ImportError)):
        read_single_cell(p)
    meta = inspect_data(str(p))
    assert meta["platform"] == "loom"
    assert meta["exists"] is True


def test_load_block_uses_io_for_new_formats():
    assert "read_single_cell" in _load_block("x.loom")
    assert "read_single_cell" in _load_block("outs/filtered_feature_bc_matrix")
    assert "read_single_cell" in _load_block("a.h5ad,b.h5ad")
    assert "read_single_cell" in _load_block("obj.rds")
    assert "sc.read_h5ad" in _load_block("x.h5ad")


def test_bbknn_integrator_template_and_route():
    from agents.dependencies import resolve_route
    from agents.reviewer import audit_code

    assert choose_integrator({"n_samples": 3, "need_batch_correction": True}, "bbknn") == "bbknn"
    assert choose_integrator({"n_samples": 3, "need_batch_correction": True}, "scanorama") == "cca"
    route = resolve_route(["clustering"], integrator="bbknn")
    assert "bbknn" in route
    assert route.index("bbknn") < route.index("neighbors")
    down = cluster_annotate_script(
        {"data_path": "x.h5ad", "need_batch_correction": True, "sample_key": "batch", "tissue": "pbmc"},
        {},
        {"integrator": "bbknn"},
    )
    compile(down, "<bbknn>", "exec")
    assert "bbknn" in down.lower()
    r = audit_code(down, {"need_batch_correction": True, "tissue": "pbmc"}, phase="downstream")
    assert r["passed"] is True


def test_writer_includes_integrator_reason():
    from agents.writer import render_report

    state = {
        "user_query": "demo",
        "report_lang": "zh",
        "thread_id": "t1",
        "metadata": {"species": "human", "tissue": "pbmc", "sample_key": "batch", "n_samples": 3},
        "plan": {
            "integrator": "harmony",
            "integrator_reason": "检测到批次列 batch（n_samples=3）：auto 选 Harmony",
            "sample_key": "batch",
            "n_samples": 3,
            "skills": [],
            "route": ["qc"],
        },
        "review_qc": {"passed": True, "issues": [], "issue_records": []},
        "artifacts": {"metrics": {"ilisi": 0.9, "kbet": 0.6, "pca_batch_r2": 0.1}, "figure_captions": []},
    }
    md = render_report(state)
    assert "检测到批次列 batch" in md
    assert "Harmony" in md
    assert "iLISI" in md
