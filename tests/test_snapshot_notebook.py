import json
from pathlib import Path

from scagent.config import load_config
from scagent.export_nb import build_notebook, build_rmd, export_analysis_notebook, package_versions
from scagent.snapshot import fork_branch, hardlink_or_copy, list_snapshots, record_h5ad


def _cfg(tmp_path):
    base = load_config()
    return {
        **base,
        "_root": str(tmp_path),
        "paths": {**base["paths"], "cache": "cache", "workspace": "ws", "outputs": "out"},
    }


def test_hardlink_or_copy_shares_inode(tmp_path):
    src = tmp_path / "a.bin"
    src.write_bytes(b"x" * 1000)
    dest = tmp_path / "b.bin"
    kind = hardlink_or_copy(src, dest)
    assert dest.is_file()
    if kind == "link":
        assert src.stat().st_ino == dest.stat().st_ino
        assert src.stat().st_nlink >= 2


def test_record_and_fork_does_not_duplicate_bytes(tmp_path):
    cfg = _cfg(tmp_path)
    src = tmp_path / "adata_qc.h5ad"
    src.write_bytes(b"fake-h5ad" * 200)
    e1 = record_h5ad(src, step="qc", thread_id="t-src", params={"seed": 0}, cfg=cfg)
    assert e1["kind"] in {"link", "copy"}
    listed = list_snapshots("t-src", cfg)
    assert listed and listed[0]["step"] == "qc"
    payload = fork_branch("t-src", "t-exp", from_step="qc", cfg=cfg)
    assert payload["entries"]
    p1 = Path(e1["path"])
    p2 = Path(payload["entries"][0]["path"])
    assert p2.is_file()
    if p1.exists() and payload["entries"][0]["kind"] == "link":
        assert p1.stat().st_ino == p2.stat().st_ino


def test_obs_delta_when_x_unchanged(tmp_path):
    import pytest

    pytest.importorskip("anndata")
    np = pytest.importorskip("numpy")
    from anndata import AnnData

    cfg = _cfg(tmp_path)
    X = np.arange(24, dtype=np.float32).reshape(6, 4)
    a = AnnData(X.copy())
    a.obs["sample"] = ["s1"] * 6
    parent_path = tmp_path / "parent.h5ad"
    a.write_h5ad(parent_path)
    parent = record_h5ad(parent_path, step="qc", thread_id="t-d", cfg=cfg)
    b = AnnData(X.copy())
    b.obs["sample"] = ["s1"] * 6
    b.obs["cell_type"] = ["T"] * 6
    child_path = tmp_path / "child.h5ad"
    b.write_h5ad(child_path)
    child = record_h5ad(child_path, step="downstream", thread_id="t-d", parent=parent, cfg=cfg)
    assert child["kind"] == "delta"
    assert child.get("obs")
    assert Path(child["obs"]).stat().st_size < Path(child_path).stat().st_size


def test_notebook_contains_seed_params_versions():
    state = {
        "user_query": "对 PBMC 做标准注释",
        "thread_id": "t-nb",
        "data_path": "data.h5ad",
        "metadata": {"tissue": "pbmc", "species": "human", "platform": "10x"},
        "plan": {"route": ["qc", "leiden", "annotate"], "integrator": None},
        "code_qc": "import numpy as np\nnp.random.seed(0)\nprint('qc')\n",
        "code_downstream": "print('down')\n",
        "artifacts": {"metrics": {"seed": 0}, "skills_fingerprint": "abc"},
    }
    nb = build_notebook(state)
    assert nb["nbformat"] == 4
    blob = json.dumps(nb, ensure_ascii=False)
    assert "seed" in blob.lower()
    assert "n_pcs" in blob
    assert "scanpy" in blob
    assert "print('qc')" in blob
    assert "print('down')" in blob
    assert "[结论]" in blob
    assert "[代码]" in blob or any(c.get("cell_type") == "code" for c in nb["cells"])
    md = "\n".join(c["source"] for c in nb["cells"] if c.get("cell_type") == "markdown")
    assert "[结论] QC" in md
    assert nb["metadata"]["scagent"]["seed"] == 0
    assert "python" in package_versions()
    assert package_versions()["scagent"]


def test_export_ipynb_and_rmd(tmp_path):
    py = export_analysis_notebook(
        {
            "user_query": "x",
            "code_qc": "print(1)",
            "plan": {"route": ["qc"]},
            "artifacts": {"metrics": {"seed": 0}},
        },
        tmp_path,
    )
    assert py.name == "analysis.ipynb"
    data = json.loads(py.read_text(encoding="utf-8"))
    assert data["cells"]
    rmd = export_analysis_notebook(
        {"user_query": "seurat", "r_degraded": True, "language": "r", "plan": {"narrative": "plan only", "route": ["plan_only"]}},
        tmp_path,
    )
    assert rmd.name == "analysis.Rmd"
    text = rmd.read_text(encoding="utf-8")
    assert "set.seed" in text
    assert "sessionInfo" in text
    assert "[结论]" in text
    assert "[代码]" in text
    assert "CreateSeuratObject" in text
    assert "FindClusters" in text
    assert build_rmd({"language": "r", "plan": {}})


def test_cli_snapshot_branch_parsers():
    from scagent.cli import build_parser

    p = build_parser()
    a = p.parse_args(["snapshots", "--thread-id", "t1"])
    assert a.func.__name__ == "cmd_snapshots"
    b = p.parse_args(["branch", "--from-thread", "t1", "--as", "t2", "--step", "qc"])
    assert b.as_name == "t2"
    assert b.func.__name__ == "cmd_branch"
