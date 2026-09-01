"""IO for AnnData (.h5ad / .loom / 10x / Cell Ranger outs) and Seurat (.rds / .h5seurat)."""

from __future__ import annotations

import gzip
import shutil
import subprocess
import tempfile
from glob import glob
from pathlib import Path
from typing import Any, Iterable

from scagent.config import performance_params
from scagent.logutil import get_logger

log = get_logger("io")

R_IO = Path(__file__).resolve().parent / "r" / "io.R"

_MATRIX_NAMES = ("matrix.mtx.gz", "matrix.mtx")
_CR_H5_NAMES = ("filtered_feature_bc_matrix.h5", "raw_feature_bc_matrix.h5")
_MATRIX_DIR_NAMES = ("filtered_feature_bc_matrix", "raw_feature_bc_matrix")
_FILE_SUFFIXES = {".h5ad", ".loom", ".h5", ".rds", ".h5seurat"}


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


def sanitize_sparse_x(adata):
    """Drop out-of-bounds sparse indices (some h5ad exports have col idx == n_vars)."""
    if getattr(adata, "isbacked", False):
        return adata
    try:
        from scipy import sparse
    except ImportError:
        return adata
    X = adata.X
    if X is None or not sparse.issparse(X):
        return adata
    n_vars = int(getattr(adata, "n_vars", 0) or 0)
    if n_vars <= 0:
        return adata
    X = X.tocsr(copy=True)
    bad = X.indices >= n_vars
    if not bad.any():
        adata.X = X
        return adata
    n_bad = int(bad.sum())
    log.warning(
        "sparse X has %s out-of-bounds column indices (max=%s, n_vars=%s); dropping entries",
        n_bad,
        int(X.indices.max()),
        n_vars,
    )
    X.data = X.data.copy()
    X.indices = X.indices.copy()
    X.data[bad] = 0
    X.indices[bad] = 0
    X.eliminate_zeros()
    adata.X = X
    return adata


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
    else:
        adata.X = sparse.csr_matrix(X)
        log.info("converted dense X to CSR nnz=%s", adata.X.nnz)
    return sanitize_sparse_x(adata)


def parse_data_spec(raw: str | Path | Iterable[str | Path] | None) -> list[Path]:
    """Split comma-separated paths and expand globs. Order is preserved."""
    if raw is None:
        return []
    if isinstance(raw, (list, tuple)):
        items = [str(x).strip() for x in raw if str(x).strip()]
    else:
        s = str(raw).strip()
        if not s:
            return []
        items = [p.strip() for p in s.split(",") if p.strip()]
    out: list[Path] = []
    for item in items:
        if any(ch in item for ch in "*?["):
            matches = sorted(glob(item, recursive=True))
            if matches:
                out.extend(Path(m).expanduser() for m in matches)
            else:
                out.append(Path(item).expanduser())
        else:
            out.append(Path(item).expanduser())
    seen: set[str] = set()
    uniq: list[Path] = []
    for p in out:
        key = str(p.resolve()) if p.exists() else str(p)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(p)
    return uniq


def _has_10x_mtx(path: Path) -> bool:
    if not path.is_dir():
        return False
    return any((path / n).exists() for n in _MATRIX_NAMES)


def resolve_10x_matrix_dir(path: str | Path) -> Path | None:
    """Cell Ranger outs/ or a 10x mtx folder. Prefers filtered_feature_bc_matrix."""
    path = Path(path)
    if path.is_file():
        return None
    candidates = [
        path,
        path / "filtered_feature_bc_matrix",
        path / "raw_feature_bc_matrix",
        path / "outs" / "filtered_feature_bc_matrix",
        path / "outs" / "raw_feature_bc_matrix",
        path / "outs",
    ]
    for c in candidates:
        if _has_10x_mtx(c):
            return c
    return None


def resolve_10x_h5(path: str | Path) -> Path | None:
    """Cell Ranger filtered_feature_bc_matrix.h5 under the given path or outs/."""
    path = Path(path)
    if path.is_file() and path.suffix.lower() == ".h5":
        return path
    if not path.is_dir():
        return None
    for folder in (path, path / "outs"):
        for name in _CR_H5_NAMES:
            hit = folder / name
            if hit.is_file():
                return hit
    return None


def sample_label(path: str | Path) -> str:
    """Folder name of a Cell Ranger sample, else file stem."""
    p = Path(path)
    name = p.name
    if name in _MATRIX_DIR_NAMES:
        parent = p.parent
        if parent.name == "outs":
            return parent.parent.name or p.stem
        return parent.name or p.stem
    if name == "outs":
        return p.parent.name or p.stem
    if name in _CR_H5_NAMES:
        parent = p.parent
        if parent.name == "outs":
            return parent.parent.name or p.stem
        return parent.name or p.stem
    if p.is_file():
        return p.stem
    return name or p.stem


def discover_samples(path: str | Path) -> list[Path]:
    """One 10x/Cell Ranger sample, or children that each look like a sample."""
    path = Path(path)
    if not path.exists():
        return [path]
    if path.is_file():
        return [path]
    mtx = resolve_10x_matrix_dir(path)
    if mtx is not None:
        return [mtx]
    h5 = resolve_10x_h5(path)
    if h5 is not None:
        return [h5]
    kids: list[Path] = []
    try:
        children = sorted(path.iterdir())
    except OSError:
        return [path]
    for child in children:
        if child.name.startswith("."):
            continue
        if child.is_file() and child.suffix.lower() in _FILE_SUFFIXES:
            kids.append(child)
            continue
        if not child.is_dir():
            continue
        mtx = resolve_10x_matrix_dir(child)
        if mtx is not None:
            kids.append(mtx)
            continue
        h5 = resolve_10x_h5(child)
        if h5 is not None:
            kids.append(h5)
    return kids


def count_tsv_rows(directory: str | Path, names: tuple[str, ...]) -> int | None:
    directory = Path(directory)
    for name in names:
        f = directory / name
        if not f.is_file():
            continue
        try:
            if name.endswith(".gz"):
                with gzip.open(f, "rt") as fh:
                    return sum(1 for line in fh if line.strip())
            with f.open(encoding="utf-8", errors="replace") as fh:
                return sum(1 for line in fh if line.strip())
        except OSError:
            return None
    return None


def peek_10x_h5_shape(path: str | Path) -> tuple[int | None, int | None]:
    """n_cells, n_genes from Cell Ranger HDF5."""
    path = Path(path)
    try:
        import h5py
    except ImportError:
        return None, None
    try:
        with h5py.File(path, "r") as f:
            matrix = f.get("matrix")
            if matrix is None:
                return None, None
            if "shape" in matrix:
                sh = matrix["shape"][()]
                n_genes, n_cells = int(sh[0]), int(sh[1])
                return n_cells, n_genes
            if "barcodes" in matrix:
                n_cells = int(len(matrix["barcodes"]))
                n_genes = int(len(matrix["features"]["id"])) if "features" in matrix else None
                return n_cells, n_genes
    except Exception as exc:
        log.debug("peek_10x_h5_shape failed: %s", exc)
    return None, None


def peek_loom(path: str | Path) -> dict[str, Any]:
    """Cell/gene counts and col_attrs keys without loading the sparse matrix."""
    path = Path(path)
    out: dict[str, Any] = {"n_cells": None, "n_genes": None, "obs_columns": [], "obs_nunique": {}}
    try:
        import h5py
    except ImportError:
        return out
    try:
        with h5py.File(path, "r") as f:
            matrix = f.get("matrix")
            if matrix is not None and getattr(matrix, "shape", None) is not None and len(matrix.shape) == 2:
                out["n_genes"] = int(matrix.shape[0])
                out["n_cells"] = int(matrix.shape[1])
            cols = f.get("col_attrs")
            if cols is not None:
                names = [str(k) for k in cols.keys()]
                out["obs_columns"] = names
                nunique: dict[str, int] = {}
                for key in names:
                    try:
                        vals = cols[key][()]
                        nunique[key] = len({str(v) for v in (vals.tolist() if hasattr(vals, "tolist") else vals)})
                    except Exception:
                        continue
                out["obs_nunique"] = nunique
            rows = f.get("row_attrs")
            genes: list[str] = []
            if rows is not None:
                for gkey in ("Gene", "gene_symbols", "var_names", "Accession"):
                    if gkey in rows:
                        raw = rows[gkey][:8000]
                        genes = [x.decode("utf-8", "replace") if isinstance(x, (bytes, bytearray)) else str(x) for x in raw]
                        break
            out["genes"] = genes
    except Exception as exc:
        log.debug("peek_loom failed: %s", exc)
    return out


def read_single_cell(path: str | Path | Iterable[str | Path], **kwargs: Any):
    """Dispatch by suffix/directory. Comma-separated or a sample folder concatenates with obs[sample]."""
    sample_key = str(kwargs.pop("sample_key", None) or "sample")
    specs = parse_data_spec(path)
    if not specs:
        raise FileNotFoundError("empty data path")
    expanded: list[Path] = []
    for spec in specs:
        if not spec.exists():
            raise FileNotFoundError(spec)
        found = discover_samples(spec)
        if spec.is_dir() and not found:
            raise ValueError(
                f"unsupported directory (not 10x mtx, Cell Ranger outs/, h5ad, or loom): {spec}"
            )
        expanded.extend(found)
    seen: set[str] = set()
    uniq: list[Path] = []
    for p in expanded:
        key = str(p.resolve()) if p.exists() else str(p)
        if key in seen:
            continue
        seen.add(key)
        uniq.append(p)
    if not uniq:
        raise ValueError(f"unsupported single-cell format: {path}")
    if len(uniq) == 1:
        return _read_one(uniq[0], **kwargs)
    log.info("concat %s samples sample_key=%s", len(uniq), sample_key)
    return concat_samples(uniq, sample_key=sample_key, **kwargs)


def _read_one(path: Path, **kwargs: Any):
    suf = path.suffix.lower()
    if path.is_dir() or _has_10x_mtx(path):
        mtx = resolve_10x_matrix_dir(path)
        if mtx is not None:
            vn = kwargs.get("var_names", "gene_symbols")
            return read_10x(mtx, var_names=str(vn))
        h5 = resolve_10x_h5(path)
        if h5 is not None:
            return read_10x_h5(h5)
        raise ValueError(f"unsupported directory: {path}")
    if suf == ".h5ad":
        return read_h5ad(path, **{k: v for k, v in kwargs.items() if k in {"backed"}})
    if suf == ".loom":
        return read_loom(path)
    if suf == ".h5":
        return read_10x_h5(path)
    if suf in {".rds", ".h5seurat"}:
        return read_seurat_rds(path, **{k: v for k, v in kwargs.items() if k in {"tmp_dir"}})
    raise ValueError(f"unsupported single-cell format: {path}")


def concat_samples(paths: list[Path], *, sample_key: str = "sample", **kwargs: Any):
    import anndata as ad
    import numpy as np

    ads = []
    for p in paths:
        a = _read_one(Path(p), **kwargs)
        label = sample_label(p)
        a.obs[sample_key] = str(label)
        a.obs_names = [f"{label}_{n}" for n in map(str, a.obs_names)]
        ads.append(a)
    try:
        out = ad.concat(ads, join="inner", index_unique=None)
    except Exception as exc:
        log.debug("concat inner failed (%s); trying outer", exc)
        out = ad.concat(ads, join="outer", index_unique=None)
    if int(getattr(out, "n_vars", 0) or 0) == 0:
        out = ad.concat(ads, join="outer", index_unique=None)
    try:
        from scipy import sparse

        X = out.X
        if X is not None and not sparse.issparse(X):
            out.X = np.nan_to_num(np.asarray(X), nan=0.0)
        elif X is not None:
            X = X.tocsr()
            if getattr(X, "data", None) is not None and np.isnan(X.data).any():
                X.data = np.nan_to_num(X.data, nan=0.0)
            out.X = X
    except Exception:
        pass
    out.obs_names_make_unique()
    out.uns["scagent_concat"] = {"n_samples": len(paths), "sample_key": sample_key, "paths": [str(p) for p in paths]}
    return ensure_csr(out)


def _drop_reserved_obs_columns(adata):
    """Some exports store index as obs['_index'], which AnnData cannot write back."""
    cols = getattr(getattr(adata, "obs", None), "columns", None)
    if cols is None or "_index" not in cols:
        return adata
    log.warning("obs column '_index' is reserved; renaming to 'orig_index'")
    adata.obs.rename(columns={"_index": "orig_index"}, inplace=True)
    return adata


def read_h5ad(path: str | Path, *, backed: str | None | bool = None, use_dask: bool | None = None):
    """Load h5ad. backed=True/'r' or auto when n_obs >= performance.backed_threshold_cells.
    use_dask=True (or performance.dask.enabled + threshold) tags experimental Dask/out-of-core path."""
    import anndata as ad

    from scagent.performance import configure_scanpy_dask, dask_params

    path = Path(path)
    n_obs, _n_vars = peek_h5ad_shape(path)
    dp = dask_params()
    want_dask = use_dask if use_dask is not None else (dp["enabled"] and n_obs is not None and n_obs >= dp["threshold_cells"])
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
    if want_dask and mode is None and n_obs is not None:
        mode = "r"
    log.info("read h5ad %s backed=%s n_obs=%s dask=%s", path, mode, n_obs, want_dask)
    adata = ad.read_h5ad(path, backed=mode) if mode else ad.read_h5ad(path)
    _drop_reserved_obs_columns(adata)
    if want_dask:
        configure_scanpy_dask(adata)
    if mode is None:
        ensure_csr(adata)
    elif n_obs:
        log.info("AnnData backed mode (avoid full RAM). Materialize a subset before scale/PCA.")
    return adata


def read_loom(path: str | Path):
    """Read .loom via anndata/scanpy. Requires loompy."""
    path = Path(path)
    log.info("read loom %s", path)
    try:
        try:
            from anndata.io import read_loom as _read_loom
        except ImportError:
            import anndata as ad

            _read_loom = ad.read_loom
        try:
            adata = _read_loom(path, sparse=True, cleanup=True)
        except TypeError:
            adata = _read_loom(path)
    except ImportError as exc:
        raise RuntimeError(
            f"Cannot read {path}: install loompy (pip install loompy) for .loom support."
        ) from exc
    except Exception as exc:
        msg = str(exc).lower()
        if "loompy" in msg:
            raise RuntimeError(
                f"Cannot read {path}: install loompy (pip install loompy) for .loom support. ({exc})"
            ) from exc
        try:
            import scanpy as sc

            adata = sc.read_loom(path)
        except Exception as exc2:
            raise RuntimeError(
                f"Cannot read {path}: {exc2}. .loom needs loompy (`pip install loompy`)."
            ) from exc2
    adata.var_names_make_unique()
    return ensure_csr(adata)


def read_10x(path: str | Path, *, var_names: str = "gene_symbols"):
    import scanpy as sc

    path = Path(path)
    resolved = resolve_10x_matrix_dir(path) or path
    log.info("read 10x mtx %s", resolved)
    adata = sc.read_10x_mtx(resolved, var_names=var_names, cache=True)
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
