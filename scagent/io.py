"""IO for AnnData (.h5ad / 10x) and Seurat (.rds / .h5seurat)."""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from scagent.config import performance_params
from scagent.logutil import get_logger

log = get_logger("io")

R_IO = Path(__file__).resolve().parent / "r" / "io.R"


def peek_h5ad_shape(path: str | Path) -> tuple[int | None, int | None]:
    """n_obs, n_vars without loading the matrix. None if unreadable."""
    path = Path(path)
    try:
        import h5py
    except ImportError:
        return None, None
    try:
        with h5py.File(path, "r") as f:
            x = f.get("X")
            if x is None:
                return None, None
            shape = getattr(x, "shape", None)
            if shape is not None and len(shape) == 2:
                return int(shape[0]), int(shape[1])
            attrs = getattr(x, "attrs", None)
            if attrs is not None and "shape" in attrs:
                sh = attrs["shape"]
                return int(sh[0]), int(sh[1])
    except Exception as exc:
        log.debug("peek_h5ad_shape failed: %s", exc)
    return None, None


def ensure_csr(adata):
    """Keep counts sparse. No-op for backed objects."""
    if getattr(adata, "isbacked", False):
        return adata
    try:
        from scipy import sparse
    except ImportError:
        return adata
    X = adata.X
    if X is None:
        return adata
    if sparse.issparse(X):
        if not sparse.isspmatrix_csr(X):
            adata.X = X.tocsr()
        return adata
    adata.X = sparse.csr_matrix(X)
    log.info("converted dense X to CSR nnz=%s", adata.X.nnz)
    return adata


def read_single_cell(path: str | Path, **kwargs: Any):
    """Dispatch reader by suffix/directory. Returns AnnData."""
    p = Path(path)
    if not p.exists():
        raise FileNotFoundError(p)
    suf = p.suffix.lower()
    if p.is_dir() or (p / "matrix.mtx").exists() or (p / "matrix.mtx.gz").exists():
        return read_10x(p, **kwargs)
    if suf == ".h5ad":
        return read_h5ad(p, **kwargs)
    if suf in {".h5"}:
        return read_10x_h5(p, **kwargs)
    if suf in {".rds", ".h5seurat"}:
        return read_seurat_rds(p, **kwargs)
    raise ValueError(f"unsupported single-cell format: {p}")


def read_h5ad(path: str | Path, *, backed: str | None | bool = None):
    """Load h5ad. backed=True/'r' or auto when n_obs >= performance.backed_threshold_cells."""
    import anndata as ad

    path = Path(path)
    n_obs, _n_vars = peek_h5ad_shape(path)
    mode: str | None
    if backed is True:
        mode = "r"
    elif backed is False:
        mode = None
    elif backed is None:
        thr = performance_params()["backed_threshold_cells"]
        mode = "r" if (n_obs is not None and n_obs >= thr) else None
    else:
        mode = str(backed)
    log.info("read h5ad %s backed=%s n_obs=%s", path, mode, n_obs)
    adata = ad.read_h5ad(path, backed=mode) if mode else ad.read_h5ad(path)
    if mode is None:
        ensure_csr(adata)
    elif n_obs:
        log.info("AnnData backed mode (avoid full RAM). Materialize a subset before scale/PCA.")
    return adata


def read_10x(path: str | Path, *, var_names: str = "gene_symbols"):
    import scanpy as sc

    log.info("read 10x mtx %s", path)
    adata = sc.read_10x_mtx(path, var_names=var_names, cache=True)
    adata.var_names_make_unique()
    return ensure_csr(adata)


def read_10x_h5(path: str | Path):
    import scanpy as sc

    log.info("read 10x h5 %s", path)
    adata = sc.read_10x_h5(path)
    adata.var_names_make_unique()
    return ensure_csr(adata)


def write_h5ad(adata, path: str | Path) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(path)
    log.info("wrote h5ad %s", path)
    return path


def read_seurat_rds(path: str | Path, *, tmp_dir: str | Path | None = None):
    """Convert Seurat .rds/.h5seurat to AnnData. Tries rpy2, then Rscript + zellkonverter."""
    path = Path(path)
    log.info("read Seurat object %s", path)
    try:
        return _read_rds_rpy2(path)
    except Exception as exc:
        log.debug("rpy2 path failed: %s", exc)
    return _read_rds_rscript(path, tmp_dir=tmp_dir)


def write_seurat_rds(adata, path: str | Path) -> Path:
    """Write AnnData to .rds via R (zellkonverter). Optional; requires R."""
    path = Path(path)
    tmp = path.with_suffix(".h5ad")
    write_h5ad(adata, tmp)
    if not shutil.which("Rscript"):
        raise RuntimeError("Rscript not found; cannot write .rds. Keep the .h5ad at " + str(tmp))
    subprocess.run(
        ["Rscript", str(R_IO), "write", str(tmp), str(path)],
        check=True,
        capture_output=True,
        text=True,
    )
    log.info("wrote rds %s", path)
    return path


def _read_rds_rpy2(path: Path):
    import anndata2ri
    import rpy2.robjects as ro
    from rpy2.robjects.packages import importr

    anndata2ri.activate()
    importr("Seurat")
    ro.r(f"obj <- readRDS('{path.as_posix()}')")
    # Prefer converting via SCE if available
    try:
        ro.r("sce <- Seurat::as.SingleCellExperiment(obj)")
        adata = ro.r("sce")
    except Exception:
        adata = ro.r("obj")
    anndata2ri.deactivate()
    return adata


def _read_rds_rscript(path: Path, tmp_dir: str | Path | None = None):
    if not shutil.which("Rscript"):
        raise RuntimeError(
            f"Cannot read {path}: install rpy2+anndata2ri, or R with zellkonverter, "
            "or convert with SeuratDisk to .h5ad first."
        )
    tmp = Path(tmp_dir) if tmp_dir else Path(tempfile.mkdtemp(prefix="scagent_rds_"))
    tmp.mkdir(parents=True, exist_ok=True)
    out = tmp / (path.stem + ".h5ad")
    proc = subprocess.run(
        ["Rscript", str(R_IO), "read", str(path), str(out)],
        capture_output=True,
        text=True,
    )
    if proc.returncode != 0 or not out.exists():
        raise RuntimeError(proc.stderr or proc.stdout or f"R conversion failed for {path}")
    return read_h5ad(out)
