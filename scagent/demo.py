"""Tiny sparse demo / test AnnData (100 cells)."""

from __future__ import annotations

from pathlib import Path

from scagent.config import REPO_ROOT
from scagent.logutil import get_logger

log = get_logger("demo")

DEFAULT_PATH = REPO_ROOT / "tests" / "data" / "tiny_100cells.h5ad"


def write_tiny_h5ad(
    path: str | Path | None = None,
    *,
    n_cells: int = 100,
    n_genes: int = 60,
) -> Path:
    import numpy as np
    from anndata import AnnData
    from scipy import sparse

    out = Path(path) if path else DEFAULT_PATH
    out.parent.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(0)
    dense = rng.poisson(1.5, size=(n_cells, n_genes)).astype(np.float32)
    X = sparse.csr_matrix(dense)
    var_names = [f"Gene{i}" for i in range(n_genes)]
    var_names[0] = "MT-ND1"
    var_names[1] = "MT-CO1"
    var_names[2] = "CD3D"
    var_names[3] = "MS4A1"
    adata = AnnData(X)
    adata.obs_names = [f"c{i}" for i in range(n_cells)]
    adata.var_names = var_names
    adata.obs["sample"] = ["s1" if i < n_cells // 2 else "s2" for i in range(n_cells)]
    adata.write_h5ad(out)
    log.info("wrote demo %s (%s cells, sparse CSR)", out, n_cells)
    return out
