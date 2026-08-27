import pytest

from scagent.cache import load_json, save_json
from scagent.demo import write_tiny_h5ad
from scagent.io import ensure_csr, peek_h5ad_shape, read_h5ad
from scagent.parallel import map_parallel, n_jobs
from agents.templates import cluster_annotate_script, qc_preprocess_script


def test_map_parallel_serial_and_threaded():
    assert n_jobs(1) == 1
    out = map_parallel(lambda x: x * 2, [1, 2, 3], jobs=1)
    assert out == [2, 4, 6]
    out2 = map_parallel(lambda x: x + 1, [10, 20], jobs=2)
    assert out2 == [11, 21]


def test_json_cache_roundtrip(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    from scagent import cache as cache_mod

    monkeypatch.setattr(cache_mod, "cache_dir", lambda cfg=None: tmp_path / "c")
    monkeypatch.setattr(cache_mod, "cache_enabled", lambda cfg=None: True)
    (tmp_path / "c").mkdir()
    save_json("hello", {"a": 1})
    assert load_json("hello") == {"a": 1}
    assert load_json("missing") is None


def test_tiny_h5ad_sparse_csr(tmp_path):
    pytest.importorskip("anndata")
    pytest.importorskip("scipy")
    path = write_tiny_h5ad(tmp_path / "tiny_100cells.h5ad", n_cells=100, n_genes=60)
    assert path.exists()
    n_obs, n_vars = peek_h5ad_shape(path)
    assert n_obs == 100
    assert n_vars == 60
    adata = read_h5ad(path, backed=False)
    from scipy import sparse

    assert sparse.issparse(adata.X)
    dense = adata.copy()
    dense.X = adata.X.toarray()
    ensure_csr(dense)
    assert sparse.isspmatrix_csr(dense.X)
    assert adata.n_obs == 100


def test_templates_have_n_jobs_sparse_cache():
    qc = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc", "n_cells": 10},
        {"method": "mad", "nmads": 5},
    )
    assert "sc.settings.n_jobs" in qc
    assert "sp.csr_matrix" in qc
    assert "CACHE_ON" in qc
    assert ".cache/adata_qc.h5ad" in qc
    big = qc_preprocess_script(
        {"data_path": "big.h5ad", "species": "human", "tissue": "pbmc", "n_cells": 500000},
        {"nmads": 5},
    )
    assert "backed='r'" in big or 'backed="r"' in big
    down = cluster_annotate_script({"data_path": "x.h5ad", "tissue": "pbmc"}, {}, {"integrator": None})
    assert "after_cluster.h5ad" in down
    assert "sc.settings.n_jobs" in down
    compile(qc, "<qc_perf>", "exec")
    compile(down, "<down_perf>", "exec")
