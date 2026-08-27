"""AST + DAG schema for generated Scanpy / Seurat scripts. No single call may skip prerequisites."""

from __future__ import annotations

import ast
import re
from typing import Any

# step -> calls that must already have appeared (any one of a dimred group is enough)
_DIMRED = ("pca", "harmony", "scvi", "bbknn")

CODE_PREREQS: dict[str, tuple[str, ...]] = {
    "neighbors": _DIMRED,
    "umap": ("neighbors",),
    "leiden": ("neighbors",),
    "cluster_deg": ("pca", "leiden", "umap"),
    "pseudobulk_deg": ("leiden", "umap"),
    "annotate": ("leiden",),
    "trajectory": ("pca", "umap", "leiden"),
}

# (attr suffix) -> step. Matched against reversed Attribute chain.
_PY_CALLS: dict[tuple[str, ...], str] = {
    ("pp", "pca"): "pca",
    ("tl", "pca"): "pca",
    ("pp", "neighbors"): "neighbors",
    ("tl", "leiden"): "leiden",
    ("tl", "louvain"): "leiden",
    ("tl", "umap"): "umap",
    ("tl", "rank_genes_groups"): "cluster_deg",
    ("tl", "dpt"): "trajectory",
    ("tl", "diffmap"): "trajectory",
    ("tl", "paga"): "trajectory",
    ("tl", "draw_graph"): "trajectory",
    ("tl", "velocity"): "trajectory",
    ("tl", "velocity_graph"): "trajectory",
    ("tl", "recover_dynamics"): "trajectory",
    ("pp", "harmony_integrate"): "harmony",
    ("harmony_integrate",): "harmony",
    ("pp", "bbknn"): "neighbors",
    ("bbknn",): "neighbors",
}

_PY_NAME_CALLS = {
    "run_trajectory_phase": "trajectory",
    "pseudobulk_de": "pseudobulk_deg",
    "fuse_annotation": "annotate",
    "rank_genes": "cluster_deg",
}

_R_PATTERNS: tuple[tuple[str, str], ...] = (
    (r"\bRunPCA\s*\(", "pca"),
    (r"\bFindNeighbors\s*\(", "neighbors"),
    (r"\bFindClusters\s*\(", "leiden"),
    (r"\bRunUMAP\s*\(", "umap"),
    (r"\bFind(?:All)?Markers\s*\(", "cluster_deg"),
    (r"\blearn_graph\s*\(|\border_cells\s*\(|\borderCells\s*\(|\bmonocle3\b|\brun_palantir\b|\bscvelo\b|\bscv\.tl\.velocity", "trajectory"),
)


def _add(records: list[dict], message: str, *, id: str) -> None:
    records.append({"id": id, "severity": "fail", "source": "schema", "message": message})


def _attr_suffix(node: ast.AST) -> tuple[str, ...]:
    parts: list[str] = []
    cur: ast.AST = node
    while isinstance(cur, ast.Attribute):
        parts.append(cur.attr)
        cur = cur.value
    return tuple(reversed(parts))


def _python_steps(tree: ast.AST) -> list[tuple[int, str]]:
    found: list[tuple[int, str]] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        lineno = getattr(node, "lineno", 0) or 0
        if isinstance(node.func, ast.Name) and node.func.id in _PY_NAME_CALLS:
            found.append((lineno, _PY_NAME_CALLS[node.func.id]))
            continue
        suf = _attr_suffix(node.func)
        if len(suf) >= 2 and suf[-2:] in _PY_CALLS:
            found.append((lineno, _PY_CALLS[suf[-2:]]))
        elif suf[-1:] in _PY_CALLS:
            found.append((lineno, _PY_CALLS[suf[-1:]]))
        elif isinstance(node.func, ast.Attribute) and node.func.attr == "annotate":
            # celltypist.annotate
            val = node.func.value
            if isinstance(val, ast.Name) and "celltypist" in val.id.lower():
                found.append((lineno, "annotate"))
        elif isinstance(node.func, ast.Attribute) and node.func.attr == "SCVI":
            found.append((lineno, "scvi"))
        elif isinstance(node.func, ast.Attribute) and node.func.attr == "get_latent_representation":
            found.append((lineno, "scvi"))
    found.sort(key=lambda x: x[0])
    return found


def _r_steps(code: str) -> list[tuple[int, str]]:
    found: list[tuple[int, str]] = []
    for i, line in enumerate(code.splitlines(), 1):
        for pat, step in _R_PATTERNS:
            if re.search(pat, line):
                found.append((i, step))
    return found


def _first_index(steps: list[tuple[int, str]]) -> dict[str, int]:
    out: dict[str, int] = {}
    for lineno, step in steps:
        out.setdefault(step, lineno)
    return out


def _prereq_ok(idx: dict[str, int], step: str, need: str) -> bool:
    here = idx[step]
    if need in _DIMRED:
        return any(d in idx and idx[d] < here for d in _DIMRED)
    if need == "pca":
        return any(d in idx and idx[d] < here for d in _DIMRED)
    return need in idx and idx[need] < here


def validate_script(code: str | None, *, phase: str = "downstream", language: str = "python") -> dict[str, Any]:
    """Return ok/issues. QC only checks syntax + no DE/DPT. Downstream enforces the analysis DAG."""
    text = code or ""
    records: list[dict] = []
    steps: list[str] = []
    ast_ok = True
    lang = (language or "python").lower()

    tree = None
    if text.strip() and lang != "r":
        try:
            tree = ast.parse(text)
        except SyntaxError as exc:
            ast_ok = False
            _add(records, f"语法错误: {exc.msg} (line {exc.lineno})", id="schema.syntax")

    ordered: list[tuple[int, str]] = []
    if tree is not None:
        ordered.extend(_python_steps(tree))
    if lang == "r" or re.search(r"\blibrary\s*\(\s*Seurat|\bFindClusters\s*\(", text):
        ordered.extend(_r_steps(text))
    ordered.sort(key=lambda x: x[0])
    steps = [s for _, s in ordered]
    idx = _first_index(ordered)

    if phase == "qc":
        for banned, label in (("cluster_deg", "差异表达"), ("trajectory", "伪时间/轨迹"), ("pseudobulk_deg", "pseudobulk DE")):
            if banned in idx:
                _add(records, f"QC 阶段出现{label}；须先完成降维聚类后再进入 downstream", id="schema.qc_order")
        ok = ast_ok and not any(r.get("severity") == "fail" for r in records)
        return {"ok": ok, "ast_ok": ast_ok, "steps": steps, "issues": [r["message"] for r in records], "issue_records": records}

    if phase == "interpret":
        ok = ast_ok and not any(r.get("severity") == "fail" for r in records)
        return {"ok": ok, "ast_ok": ast_ok, "steps": steps, "issues": [r["message"] for r in records], "issue_records": records}

    for step, prereqs in CODE_PREREQS.items():
        if step not in idx:
            continue
        missing = [p for p in prereqs if not _prereq_ok(idx, step, p)]
        if not missing:
            continue
        if step in {"cluster_deg", "pseudobulk_deg"}:
            _add(
                records,
                f"DE ({step}) 出现在降维/聚类之前或缺少前提 {missing}；须先 PCA/neighbors + Leiden/UMAP",
                id="schema.dag_de",
            )
        elif step == "trajectory":
            _add(
                records,
                f"伪时间/轨迹出现在降维聚类之前或缺少前提 {missing}；DPT/PAGA/Palantir/scVelo/Monocle3 须在 PCA + UMAP + Leiden 之后",
                id="schema.dag_traj",
            )
        else:
            _add(records, f"步骤 {step} 缺少前提 {missing}", id="schema.dag")

    ok = ast_ok and not any(r.get("severity") == "fail" for r in records)
    return {
        "ok": ok,
        "ast_ok": ast_ok,
        "steps": steps,
        "issues": [r["message"] for r in records],
        "issue_records": records,
    }
