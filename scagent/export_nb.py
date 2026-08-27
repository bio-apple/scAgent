"""Export a dual-format Jupyter notebook (or Seurat Rmd) and execute Python via nbclient."""

from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
from datetime import date
from importlib import metadata
from pathlib import Path
from typing import Any

from scagent.config import analysis_params

_PKGS = (
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "celltypist",
    "harmonypy",
    "scvi-tools",
    "leidenalg",
    "igraph",
    "squidpy",
    "nbformat",
    "nbclient",
)


def package_versions(names: tuple[str, ...] = _PKGS) -> dict[str, str]:
    from scagent.compat import scagent_version

    out: dict[str, str] = {"python": sys.version.split()[0], "scagent": scagent_version()}
    for name in names:
        try:
            out[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            out[name] = "not-installed"
    return out


def _md_cell(source: str) -> dict[str, Any]:
    text = source if source.endswith("\n") else source + "\n"
    return {"cell_type": "markdown", "metadata": {}, "source": text}


def _as_nb_text(text: str) -> list[str]:
    if not text:
        return []
    if not text.endswith("\n"):
        text += "\n"
    return [text]


def _stream_outputs(stdout: str | None = None, stderr: str | None = None) -> list[dict[str, Any]]:
    outs: list[dict[str, Any]] = []
    if stdout:
        outs.append({"output_type": "stream", "name": "stdout", "text": _as_nb_text(stdout)})
    if stderr:
        outs.append({"output_type": "stream", "name": "stderr", "text": _as_nb_text(stderr)})
    return outs


def _code_cell(
    source: str,
    *,
    stdout: str | None = None,
    stderr: str | None = None,
    execution_count: int | None = None,
) -> dict[str, Any]:
    text = source if source.endswith("\n") else source + "\n"
    outputs = _stream_outputs(stdout, stderr)
    count = execution_count if outputs else None
    if outputs and count is None:
        count = 1
    return {
        "cell_type": "code",
        "metadata": {},
        "execution_count": count,
        "outputs": outputs,
        "source": text,
    }


def _exe(state: dict, key: str) -> dict:
    return state.get(key) or {}


def _provenance_markdown(state: dict) -> str:
    meta = state.get("metadata") or {}
    plan = state.get("plan") or {}
    params = analysis_params()
    versions = package_versions()
    seed = (state.get("artifacts") or {}).get("metrics", {}).get("seed")
    if seed is None:
        seed = params.get("seed", 0)
    ver_lines = "\n".join(f"- `{k}=={v}`" for k, v in versions.items())
    param_lines = "\n".join(f"- `{k}`: {v}" for k, v in params.items())
    snaps = (state.get("execution_qc") or {}).get("snapshot_manifests") or []
    snap_txt = ", ".join(f"{s.get('step')}={s.get('kind')}" for s in snaps) or "none"
    spatial = "spatial" in (plan.get("route") or []) or "squidpy" in (plan.get("skills") or [])
    squidpy_note = (
        "\nSpatial steps use Squidpy (`import squidpy as sq`). This route has no spatial task, so Squidpy is not imported.\n"
        if not spatial
        else "\nSpatial steps: use Squidpy (`import squidpy as sq`) on the AnnData with spatial coordinates.\n"
    )
    return (
        f"# scAgent analysis notebook\n\n"
        f"Strict **code–result** layout: each phase is a markdown `[结论]` cell followed by a runnable code cell "
        f"(Scanpy; Squidpy only if spatial). Do not paste report prose into code cells.\n"
        f"{squidpy_note}"
        f"- date: {date.today().isoformat()}\n"
        f"- query: {state.get('user_query') or ''}\n"
        f"- thread_id: `{state.get('thread_id') or ''}`\n"
        f"- tissue: {meta.get('tissue')} | species: {meta.get('species')} | platform: {meta.get('platform')}\n"
        f"- data: `{state.get('data_path') or meta.get('data_path') or ''}`\n"
        f"- route: {' → '.join(plan.get('route') or [])}\n"
        f"- integrator: {plan.get('integrator')}\n"
        f"- seed: **{seed}**\n"
        f"- skills fingerprint: `{(state.get('artifacts') or {}).get('skills_fingerprint') or ''}`\n"
        f"- snapshots: {snap_txt}\n\n"
        f"Run cells from `workspace/` (scripts expect `adata_qc.h5ad`). "
        f"Do not treat exploratory Wilcoxon as a between-condition result.\n\n"
        f"## Parameters (`config.yaml`)\n\n{param_lines}\n\n"
        f"## Tool versions\n\n{ver_lines}\n"
    )


def build_notebook(state: dict) -> dict[str, Any]:
    from scagent.dual import PHASES, phase_conclusion, strip_code_fences

    cells = [_md_cell(_provenance_markdown(state))]
    exe_keys = {
        "qc": "execution_qc",
        "downstream": "execution_downstream",
        "interpret": "execution_interpret",
    }
    n = 0
    for phase, key, title in PHASES:
        code = strip_code_fences(state.get(key) or "")
        if not code:
            continue
        n += 1
        cells.append(_md_cell(f"## [结论] {title}\n\n{phase_conclusion(state, phase)}"))
        cells.append(_md_cell(f"## [代码] {title}"))
        exe = _exe(state, exe_keys[phase])
        cells.append(
            _code_cell(
                code,
                stdout=exe.get("stdout") or "",
                stderr=exe.get("stderr") or "",
                execution_count=n if exe.get("executed") else None,
            )
        )
    if n == 0:
        cells.append(_md_cell("No executable Python was generated (R path writes `analysis.Rmd`)."))
    return {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "language_info": {"name": "python", "pygments_lexer": "ipython3"},
            "scagent": {
                "thread_id": state.get("thread_id"),
                "seed": analysis_params().get("seed"),
                "versions": package_versions(),
                "dual": "code-result-v1",
            },
        },
        "cells": cells,
    }


def seurat_phase_chunks(state: dict) -> list[tuple[str, str, str]]:
    """Runnable Seurat chunks for dual Rmd. scAgent does not execute them."""
    params = analysis_params()
    seed = int(params.get("seed") or 0)
    n_pcs = int(params.get("n_pcs") or 40)
    n_hvg = int(params.get("n_hvg") or 2000)
    n_neighbors = int(params.get("n_neighbors") or 15)
    res = params.get("leiden_resolution")
    res_txt = "0.6" if res is None else str(res)
    data = state.get("data_path") or (state.get("metadata") or {}).get("data_path") or "path/to/data"
    tissue = (state.get("metadata") or {}).get("tissue") or "default"
    qc = (
        f"library(Seurat)\n"
        f"library(Matrix)\n"
        f"set.seed({seed})\n"
        f"# QC: MAD on nCount / nFeature / percent.mt. Do NOT hardcode mito < 5%.\n"
        f"# Edit `data_path` if this is not 10x mtx; .h5ad needs zellkonverter / sceasy.\n"
        f"data_path <- {json.dumps(str(data))}\n"
        f"tissue <- {json.dumps(str(tissue))}\n"
        f"nmads <- 5L\n"
        f"if (dir.exists(data_path)) {{\n"
        f"  counts <- Read10X(data.dir = data_path)\n"
        f"  obj <- CreateSeuratObject(counts = counts, min.cells = 0, min.features = 0)\n"
        f"}} else if (grepl('\\\\.rds$', data_path, ignore.case = TRUE)) {{\n"
        f"  obj <- readRDS(data_path)\n"
        f"}} else {{\n"
        f"  stop('Point data_path at 10x mtx dir or Seurat .rds')\n"
        f"}}\n"
        f"obj[[\"percent.mt\"]] <- PercentageFeatureSet(obj, pattern = \"^MT-\")\n"
        f"mad_cut <- function(x, n = nmads, lo = TRUE, hi = TRUE) {{\n"
        f"  med <- stats::median(x, na.rm = TRUE)\n"
        f"  md <- stats::mad(x, na.rm = TRUE)\n"
        f"  c(if (lo) med - n * md else -Inf, if (hi) med + n * md else Inf)\n"
        f"}}\n"
        f"cut_count <- mad_cut(obj$nCount_RNA)\n"
        f"cut_feat <- mad_cut(obj$nFeature_RNA)\n"
        f"cut_mt <- mad_cut(obj$percent.mt, lo = FALSE)\n"
        f"keep <- obj$nCount_RNA >= cut_count[1] & obj$nCount_RNA <= cut_count[2] &\n"
        f"  obj$nFeature_RNA >= cut_feat[1] & obj$nFeature_RNA <= cut_feat[2] &\n"
        f"  obj$percent.mt <= cut_mt[2]\n"
        f"obj <- obj[, keep]\n"
        f"obj <- NormalizeData(obj)\n"
        f"obj <- FindVariableFeatures(obj, nfeatures = {n_hvg})\n"
        f"obj <- ScaleData(obj)\n"
        f"obj <- RunPCA(obj, npcs = {n_pcs})\n"
        f"saveRDS(obj, file = \"workspace/adata_qc.rds\")\n"
    )
    down = (
        f"library(Seurat)\n"
        f"set.seed({seed})\n"
        f"obj <- readRDS(\"workspace/adata_qc.rds\")\n"
        f"obj <- FindNeighbors(obj, dims = 1:{n_pcs}, k.param = {n_neighbors})\n"
        f"obj <- FindClusters(obj, resolution = {res_txt})\n"
        f"obj <- RunUMAP(obj, dims = 1:{n_pcs})\n"
        f"# Cluster markers are exploratory; condition DE needs pseudobulk + FDR (edgeR/DESeq2).\n"
        f"markers <- FindAllMarkers(obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)\n"
        f"write.csv(markers, \"workspace/cluster_markers.csv\", row.names = FALSE)\n"
        f"# Monocle3 fate (optional; requires monocle3). Do not run before PCA/clusters.\n"
        f"# library(monocle3); cds <- as.cell_data_set(obj); cds <- cluster_cells(cds);\n"
        f"# cds <- learn_graph(cds); cds <- order_cells(cds)\n"
        f"saveRDS(obj, file = \"workspace/adata_annotated.rds\")\n"
    )
    interp = (
        "# Over-representation / GSEA in R (clusterProfiler or fgsea). Edit gene lists from DEG.\n"
        "# library(clusterProfiler)\n"
        '# ego <- enrichGO(gene = sig_genes, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH")\n'
        'message("Interpretation: run ORA/GSEA on pseudobulk DEG, not Wilcoxon cluster markers.")\n'
        "sessionInfo()\n"
    )
    titles = {
        "qc": "QC & Preprocessing",
        "downstream": "Clustering & Differential",
        "interpret": "Biological Interpretation",
    }
    chunks = [("qc", titles["qc"], qc), ("downstream", titles["downstream"], down), ("interpret", titles["interpret"], interp)]
    return chunks


def build_rmd(state: dict) -> str:
    """Dual-format Seurat Rmd. Not executed by scAgent (no IRkernel binding)."""
    from scagent.dual import phase_conclusion, report_lang

    plan = state.get("plan") or {}
    params = analysis_params()
    versions = package_versions()
    seed = params.get("seed", 0)
    zh = report_lang(state) != "en"
    header = (
        "---\n"
        'title: "scAgent Seurat dual notebook"\n'
        f'date: "{date.today().isoformat()}"\n'
        "output: html_document\n"
        "---\n\n"
        f"**Query:** {state.get('user_query') or ''}\n\n"
        f"**Seed:** `{seed}`\n\n"
        f"**Route:** {' → '.join(plan.get('route') or [])}\n\n"
        + (
            "严格「代码-结果」格式。scAgent **不执行** 此 Rmd；在 RStudio / IRkernel 中逐块运行。\n"
            if zh
            else "Strict code–result format. scAgent **does not execute** this Rmd; run chunks in RStudio / IRkernel.\n"
        )
        + f"\n{plan.get('narrative') or ''}\n\n"
    )
    body = []
    for phase, title, code in seurat_phase_chunks(state):
        body.append(f"## [结论] {title}\n\n{phase_conclusion(state, phase)}\n")
        body.append(f"## [代码] {title}\n\n```{{r {phase}}}\n{code.rstrip()}\n```\n")
    footer = (
        "\n## Parameters\n\n"
        + "\n".join(f"- `{k}`: {v}" for k, v in params.items())
        + "\n\n## Python-side versions (if mixed)\n\n"
        + "\n".join(f"- `{k}=={v}`" for k, v in versions.items())
        + "\n"
    )
    return header + "\n".join(body) + footer


def _collect_nbclient_io(cell: Any) -> tuple[str, str, bool]:
    stdout_parts: list[str] = []
    stderr_parts: list[str] = []
    ok = True
    for out in getattr(cell, "outputs", None) or []:
        kind = getattr(out, "output_type", None) or (out.get("output_type") if isinstance(out, dict) else None)
        if kind == "stream":
            name = getattr(out, "name", None) or (out.get("name") if isinstance(out, dict) else "stdout")
            text = getattr(out, "text", None)
            if text is None and isinstance(out, dict):
                text = out.get("text")
            blob = "".join(text) if isinstance(text, list) else str(text or "")
            (stdout_parts if name != "stderr" else stderr_parts).append(blob)
        elif kind == "error":
            ok = False
            tb = getattr(out, "traceback", None)
            if tb is None and isinstance(out, dict):
                tb = out.get("traceback")
            stderr_parts.append("\n".join(tb) if isinstance(tb, list) else str(tb or "error"))
    return "".join(stdout_parts), "".join(stderr_parts), ok


def _run_unjailed_subprocess(
    script: Path,
    *,
    workspace: Path,
    env: dict[str, str],
    timeout: int,
) -> tuple[str, str, int]:
    kwargs: dict[str, Any] = {
        "cwd": str(workspace),
        "env": env,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        kwargs["start_new_session"] = True
    proc = subprocess.Popen([sys.executable, str(script)], **kwargs)
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        if os.name == "posix":
            try:
                os.killpg(proc.pid, signal.SIGKILL)
            except (ProcessLookupError, PermissionError, OSError):
                proc.kill()
        else:
            proc.kill()
        stdout, stderr = proc.communicate(timeout=5)
        return stdout or "", (stderr or "") + f"\ntimeout after {timeout}s", -9
    return stdout or "", stderr or "", proc.returncode or 0


def _kernel_name() -> str | None:
    try:
        from jupyter_client.kernelspec import KernelSpecManager

        specs = KernelSpecManager().find_kernel_specs()
    except Exception:
        return None
    if "python3" in specs:
        return "python3"
    if "python" in specs:
        return "python"
    return None


def _kernel_ready(kernel: str, timeout: float = 6.0) -> bool:
    def _probe() -> bool:
        from jupyter_client import KernelManager

        km = KernelManager(kernel_name=kernel)
        km.start_kernel()
        try:
            kc = km.client()
            kc.start_channels()
            kc.wait_for_ready(timeout=timeout)
            return True
        finally:
            try:
                km.shutdown_kernel(now=True)
            except Exception:
                pass
            try:
                km.cleanup_resources()
            except Exception:
                pass

    try:
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=1) as pool:
            return bool(pool.submit(_probe).result(timeout=timeout + 2))
    except Exception:
        return False


def _try_nbclient(code: str, *, workspace: Path, timeout: int, env: dict[str, str]) -> dict[str, Any] | None:
    try:
        import nbformat
        from nbclient import NotebookClient
    except ImportError:
        return None
    kernel = _kernel_name()
    if not kernel or not _kernel_ready(kernel):
        return None
    nb = nbformat.v4.new_notebook()
    nb["metadata"]["kernelspec"] = {"name": kernel, "language": "python", "display_name": "Python 3"}
    nb.cells.append(nbformat.v4.new_code_cell(source=code if code.endswith("\n") else code + "\n"))
    client = NotebookClient(
        nb,
        timeout=timeout,
        kernel_name=kernel,
        resources={"metadata": {"path": str(workspace)}},
    )
    patched = ("PYTHONHASHSEED", "PYTHONPATH", "MPLCONFIGDIR")
    old = {k: os.environ.get(k) for k in patched}
    try:
        for k in patched:
            if k in env and env[k] is not None:
                os.environ[k] = env[k]
        try:
            client.execute()
        except Exception as exc:
            name = type(exc).__name__
            msg = str(exc).lower()
            if "no such kernel" in msg or name in {"NoSuchKernel", "NoSuchKernelError"}:
                return None
            stdout, stderr, _ok = _collect_nbclient_io(nb.cells[0] if nb.cells else {})
            return {
                "ok": False,
                "stdout": stdout,
                "stderr": (stderr + "\n" + str(exc)).strip(),
                "returncode": 1,
                "jail": "jupyter",
                "notebook_obj": nb,
            }
        stdout, stderr, ok = _collect_nbclient_io(nb.cells[0])
        return {
            "ok": ok,
            "stdout": stdout,
            "stderr": stderr,
            "returncode": 0 if ok else 1,
            "jail": "jupyter",
            "notebook_obj": nb,
        }
    finally:
        for k, v in old.items():
            if v is None:
                os.environ.pop(k, None)
            else:
                os.environ[k] = v


def write_executed_notebook(
    path: Path,
    code: str,
    stdout: str,
    stderr: str,
    *,
    ok: bool,
    notebook_obj: Any | None = None,
) -> Path:
    path = Path(path)
    if notebook_obj is not None:
        try:
            import nbformat

            path.write_text(nbformat.writes(notebook_obj), encoding="utf-8")
            return path
        except Exception:
            pass
    outputs = _stream_outputs(stdout, stderr)
    if not ok:
        outputs.append(
            {
                "output_type": "error",
                "ename": "ExecutionError",
                "evalue": (stderr or "nonzero exit")[:500],
                "traceback": [stderr or "nonzero exit"],
            }
        )
    nb = {
        "nbformat": 4,
        "nbformat_minor": 5,
        "metadata": {
            "kernelspec": {"display_name": "Python 3", "language": "python", "name": "python3"},
            "scagent": {"executor": "jupyter"},
        },
        "cells": [
            {
                "cell_type": "code",
                "metadata": {},
                "execution_count": 1,
                "outputs": outputs,
                "source": code if code.endswith("\n") else code + "\n",
            }
        ],
    }
    path.write_text(json.dumps(nb, ensure_ascii=False, indent=1), encoding="utf-8")
    return path


def execute_via_jupyter(
    code: str,
    *,
    workspace: Path,
    timeout: int,
    filename: str = "analysis.py",
    env: dict[str, str] | None = None,
    script: Path | None = None,
) -> dict[str, Any]:
    """Run generated Python as a notebook cell. No OS jail (seatbelt/bwrap)."""
    workspace = Path(workspace)
    workspace.mkdir(parents=True, exist_ok=True)
    (workspace / "figures").mkdir(exist_ok=True)
    script = Path(script) if script is not None else workspace / filename
    if not script.exists():
        script.write_text(code or "", encoding="utf-8")
    env = env or os.environ.copy()
    nb_path = workspace / f"{script.stem}.ipynb"
    via = _try_nbclient(code, workspace=workspace, timeout=timeout, env=env)
    if via is not None:
        write_executed_notebook(
            nb_path,
            code,
            via.get("stdout") or "",
            via.get("stderr") or "",
            ok=bool(via.get("ok")),
            notebook_obj=via.get("notebook_obj"),
        )
        via["notebook"] = str(nb_path)
        via["executed"] = True
        return via
    stdout, stderr, rc = _run_unjailed_subprocess(script, workspace=workspace, env=env, timeout=timeout)
    ok = rc == 0
    write_executed_notebook(nb_path, code, stdout, stderr, ok=ok)
    return {
        "ok": ok,
        "stdout": stdout,
        "stderr": stderr,
        "returncode": rc,
        "jail": "jupyter-subprocess",
        "notebook": str(nb_path),
        "executed": True,
    }


def export_analysis_notebook(state: dict, out_dir: Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    if (state.get("language") or (state.get("plan") or {}).get("language")) == "r" or state.get("r_degraded"):
        path = out_dir / "analysis.Rmd"
        path.write_text(build_rmd(state), encoding="utf-8")
        return path
    path = out_dir / "analysis.ipynb"
    path.write_text(json.dumps(build_notebook(state), ensure_ascii=False, indent=1), encoding="utf-8")
    return path
