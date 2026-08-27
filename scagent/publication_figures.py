"""Publication-grade figure checklist (QC violin, batch UMAP, marker heatmap, volcano, pathway bubble)."""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable

ContextFn = Callable[[dict[str, Any]], bool]


def _ctx(state: dict) -> dict[str, Any]:
    meta = state.get("metadata") or {}
    plan = state.get("plan") or {}
    arts = state.get("artifacts") or {}
    executed = bool(state.get("execute_code")) and any(
        (state.get(k) or {}).get("executed")
        for k in ("execution_qc", "execution_downstream", "execution_interpret")
    )
    route = list(plan.get("route") or [])
    return {
        "executed": executed,
        "qc_executed": bool((state.get("execution_qc") or {}).get("executed")),
        "downstream_executed": bool((state.get("execution_downstream") or {}).get("executed")),
        "interpret_executed": bool((state.get("execution_interpret") or {}).get("executed")),
        "need_batch": bool(meta.get("need_batch_correction") or int(meta.get("n_samples") or 1) > 1),
        "needs_pseudobulk": bool(plan.get("needs_pseudobulk") or meta.get("needs_pseudobulk")),
        "enrichment_planned": "enrichment" in route or "interpret" in route or bool(state.get("interpretation_plan")),
    }


def _caps_and_paths(artifacts: dict) -> tuple[list[dict], list[str]]:
    caps = list(artifacts.get("figure_captions") or [])
    paths = list(artifacts.get("figures") or [])
    for c in caps:
        p = str(c.get("path") or "")
        if p and p not in paths:
            paths.append(p)
    return caps, paths


def _match_kind(caps: list[dict], paths: list[str], kinds: tuple[str, ...]) -> dict | None:
    for c in caps:
        kind = str(c.get("kind") or "")
        if kind in kinds:
            return c
    for kind in sorted(kinds, key=len, reverse=True):
        for p in paths:
            low = Path(p).name.lower()
            if kind in low:
                return {"path": p, "kind": kind, "caption": ""}
    return None


def _match_pair(caps: list[dict], paths: list[str], kinds: tuple[str, str]) -> tuple[dict | None, dict | None]:
    return _match_kind(caps, paths, (kinds[0],)), _match_kind(caps, paths, (kinds[1],))


PUBLICATION_FIGURE_SPECS: list[dict[str, Any]] = [
    {
        "id": "qc_violin",
        "title_zh": "QC violin",
        "title_en": "QC violin",
        "desc_zh": "n_genes / total_counts / pct_mt 分布，用于 MAD 阈值设定。",
        "desc_en": "n_genes / total_counts / pct_mt distribution for MAD thresholding.",
        "required": lambda c: c["qc_executed"],
        "kinds": ("violin",),
    },
    {
        "id": "batch_umap",
        "title_zh": "整合前后 UMAP（批次着色）",
        "title_en": "Integration UMAP (before / after, batch-colored)",
        "desc_zh": "校正前/后 UMAP 批次着色对照；混匀不是整合成功证据。",
        "desc_en": "Batch-colored UMAP before vs after correction; mixing alone is not proof of success.",
        "required": lambda c: c["need_batch"] and c["downstream_executed"],
        "pair_kinds": ("batch_umap_before", "batch_umap_after"),
    },
    {
        "id": "marker_heatmap",
        "title_zh": "Marker heatmap",
        "title_en": "Marker heatmap",
        "desc_zh": "探索性 cluster marker 热图（Wilcoxon/t-test/MAST）；非组间结论。",
        "desc_en": "Exploratory cluster marker heatmap; not a between-group DE result.",
        "required": lambda c: c["downstream_executed"],
        "kinds": ("marker_heatmap", "heatmap"),
    },
    {
        "id": "volcano",
        "title_zh": "火山图（pseudobulk DEG）",
        "title_en": "Volcano plot (pseudobulk DEG)",
        "desc_zh": "组间 sample-level pseudobulk 差异表达（logFC vs FDR）。",
        "desc_en": "Sample-level pseudobulk group DE (logFC vs FDR).",
        "required": lambda c: c["needs_pseudobulk"] and c["downstream_executed"],
        "kinds": ("volcano",),
    },
    {
        "id": "pathway_bubble",
        "title_zh": "通路气泡图",
        "title_en": "Pathway bubble plot",
        "desc_zh": "ORA/GSEA 富集 top 通路（-log10 p 与 overlap 大小）。",
        "desc_en": "ORA/GSEA top pathways (-log10 p-value; bubble size = overlap).",
        "required": lambda c: c["interpret_executed"] or (c["enrichment_planned"] and c["downstream_executed"]),
        "kinds": ("pathway_bubble", "bubble"),
    },
]


def _status_label(status: str, lang: str) -> str:
    zh = lang != "en"
    return {
        "present": "✓ 已生成" if zh else "✓ present",
        "missing": "✗ 缺失" if zh else "✗ missing",
        "na": "— 不适用" if zh else "— n/a",
        "pending": "○ 未执行" if zh else "○ not executed",
    }.get(status, status)


def build_publication_figure_inventory(state: dict) -> dict[str, Any]:
    ctx = _ctx(state)
    arts = state.get("artifacts") or {}
    caps, paths = _caps_and_paths(arts)
    items: list[dict[str, Any]] = []

    for spec in PUBLICATION_FIGURE_SPECS:
        req_fn: ContextFn = spec["required"]
        applicable = bool(req_fn(ctx))
        if not ctx["executed"]:
            status = "pending"
            path_txt = ""
            detail = ""
        elif not applicable:
            status = "na"
            path_txt = ""
            detail = ""
        elif spec.get("pair_kinds"):
            a_kind, b_kind = spec["pair_kinds"]
            a, b = _match_pair(caps, paths, (a_kind, b_kind))
            detail = ""
            if a and b:
                status = "present"
                path_txt = f"{a.get('path')} ; {b.get('path')}"
            else:
                status = "missing"
                missing = []
                if not a:
                    missing.append(a_kind)
                if not b:
                    missing.append(b_kind)
                path_txt = ""
                detail = "missing: " + ", ".join(missing)
        else:
            cap = _match_kind(caps, paths, tuple(spec.get("kinds") or ()))
            detail = ""
            if cap:
                status = "present"
                path_txt = str(cap.get("path") or "")
            else:
                status = "missing"
                path_txt = ""

        items.append(
            {
                "id": spec["id"],
                "title_zh": spec["title_zh"],
                "title_en": spec["title_en"],
                "desc_zh": spec["desc_zh"],
                "desc_en": spec["desc_en"],
                "status": status,
                "path": path_txt,
                "detail": detail if spec.get("pair_kinds") else (detail or ""),
                "required": applicable and ctx["executed"],
            }
        )

    required = [i for i in items if i["required"]]
    present = [i for i in required if i["status"] == "present"]
    missing = [i for i in required if i["status"] == "missing"]
    return {
        "items": items,
        "n_required": len(required),
        "n_present": len(present),
        "n_missing": len(missing),
        "complete": len(missing) == 0,
        "missing_ids": [i["id"] for i in missing],
    }


def render_publication_figure_inventory_markdown(state: dict, *, lang: str = "zh") -> str:
    inv = build_publication_figure_inventory(state)
    zh = lang != "en"
    lines = [
        (
            "发表级主图清单：执行后应包含下列图；缺失项会在 Publication Reviewer 中标记。"
            if zh
            else "Publication main-figure checklist. Missing items are flagged in the publication review."
        ),
        "",
        "| " + ("图" if zh else "Figure") + " | " + ("说明" if zh else "Description") + " | " + ("状态" if zh else "Status") + " | " + ("路径" if zh else "Path") + " |",
        "|---|---|---|---|",
    ]
    for item in inv["items"]:
        title = item["title_zh"] if zh else item["title_en"]
        desc = item["desc_zh"] if zh else item["desc_en"]
        st = _status_label(item["status"], lang)
        path = item.get("path") or ("—" if zh else "—")
        if item["status"] == "missing" and item.get("detail"):
            st += f" ({item['detail']})"
        lines.append(f"| **{title}** | {desc} | {st} | `{path}` |")
    if zh:
        if inv["n_required"] == 0:
            summary = "**汇总**: 当前路线无必需发表级主图（或未执行）。"
        elif inv["complete"]:
            summary = f"**汇总**: {inv['n_present']}/{inv['n_required']} 必需图已就绪。"
        else:
            summary = f"**汇总**: {inv['n_present']}/{inv['n_required']} 已生成；缺失: {', '.join(inv['missing_ids'])}。"
    else:
        if inv["n_required"] == 0:
            summary = "**Summary**: no required publication figures for this route (or not executed)."
        elif inv["complete"]:
            summary = f"**Summary**: {inv['n_present']}/{inv['n_required']} required figures present."
        else:
            summary = f"**Summary**: {inv['n_present']}/{inv['n_required']} present; missing: {', '.join(inv['missing_ids'])}."
    lines.extend(["", summary, ""])
    return "\n".join(lines)
