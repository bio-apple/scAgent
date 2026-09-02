"""Strict code–result dual output: human conclusions and runnable code stay in separate blocks."""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any

PHASES: tuple[tuple[str, str, str], ...] = (
    ("qc", "code_qc", "QC & Preprocessing"),
    ("downstream", "code_downstream", "Clustering & Differential"),
    ("interpret", "code_interpret", "Biological Interpretation"),
)

_FENCE_OPEN = re.compile(r"^```(?:python|py|r|R)?[^\n]*\n")
_FENCE_CLOSE = re.compile(r"\n```\s*$")


def is_r_path(state: dict) -> bool:
    lang = state.get("language") or (state.get("plan") or {}).get("language") or "python"
    return str(lang).lower() == "r" or bool(state.get("r_degraded"))


def report_lang(state: dict) -> str:
    raw = state.get("report_lang") or "zh"
    return "en" if raw == "en" else "zh"


def strip_code_fences(code: str | None) -> str:
    s = (code or "").strip()
    if s.startswith("```"):
        s = _FENCE_OPEN.sub("", s)
        s = _FENCE_CLOSE.sub("", s)
    s = s.strip()
    return s + ("\n" if s else "")


def _executed(state: dict, phase: str) -> bool:
    key = {"qc": "execution_qc", "downstream": "execution_downstream", "interpret": "execution_interpret"}.get(phase)
    if not key:
        return False
    return bool((state.get(key) or {}).get("executed"))


def phase_conclusion(state: dict, phase: str, lang: str | None = None) -> str:
    zh = (lang or report_lang(state)) != "en"
    m = (state.get("artifacts") or {}).get("metrics") or {}
    ran = _executed(state, phase)
    pending = (
        "本阶段代码已生成，尚未执行（未加 `--execute`）。数字仅在执行后有效。"
        if zh
        else "Code was generated but not executed (no `--execute`). Numbers are valid only after a run."
    )
    if phase == "qc":
        qc = state.get("qc_strategy") or {}
        lines = [
            f"{'QC 方法' if zh else 'QC method'}: {qc.get('method') or 'mad'}。"
            + ("不使用默认 mito%<5。" if zh else " No default mito%<5 cutoff."),
        ]
        if qc.get("pct_mt_note"):
            lines.append(str(qc["pct_mt_note"]))
        if qc.get("protocol"):
            lines.append(str(qc["protocol"]).strip())
        if ran and m.get("n_before") is not None:
            lines.append(
                f"{'细胞' if zh else 'cells'}: {m.get('n_before')} → {m.get('n_after')}"
                f"（{'移除' if zh else 'removed'} {m.get('pct_removed')}%）。"
                if zh
                else f"cells: {m.get('n_before')} → {m.get('n_after')} (removed {m.get('pct_removed')}%)."
            )
            if m.get("doublet_rate") is not None:
                lines.append(f"{'双细胞比例' if zh else 'doublet rate'}: {m['doublet_rate']}")
            if m.get("doublet_n_high_conf") is not None:
                lines.append(
                    f"{'双细胞分级' if zh else 'doublet tiers'}: "
                    f"high={m.get('doublet_n_high_conf')} low={m.get('doublet_n_low_conf')} "
                    f"filter={m.get('doublet_filter') or 'high_conf'}"
                )
        elif not ran:
            lines.append(pending)
        return "\n".join(x for x in lines if x).strip() + "\n"
    if phase == "downstream":
        plan = state.get("plan") or {}
        ann = state.get("annotation_plan") or {}
        lines = [
            f"{'整合' if zh else 'integrator'}: {plan.get('integrator') or ('无 / 单样本' if zh else 'none / single sample')}。",
            f"Leiden resolution: {m.get('resolution') or plan.get('resolution') or 'adaptive'}。",
        ]
        if plan.get("integrator_reason"):
            lines.insert(1, f"{'整合理由' if zh else 'integration decision'}: {plan.get('integrator_reason')}。")
        if ran and m.get("n_clusters") is not None:
            lines.append(f"{'簇数' if zh else 'clusters'}: {m.get('n_clusters')}。")
            if m.get("deg_engine"):
                extra = f"{'DEG' if zh else 'DEG'}: {m.get('deg_engine')}（n_sig={m.get('deg_n_sig', '—')}）"
                if m.get("deg_n_overlap") is not None:
                    extra += f"；交叉验证 overlap={m.get('deg_n_overlap')}"
                lines.append(extra + "。")
            if m.get("marker_methods") or m.get("marker_method"):
                lines.append(
                    f"cluster marker: {m.get('marker_methods') or m.get('marker_method')}"
                    + (f"；overlap={m.get('marker_n_overlap')}" if m.get("marker_n_overlap") is not None else "")
                    + "。"
                )
        elif not ran:
            lines.append(pending)
        if ann.get("dual_validation"):
            lines.append(
                "注释：CellTypist + cluster DE∩catalog + marker 双验证，禁止单基因定型。"
                if zh
                else "Annotation: CellTypist + cluster DE∩catalog + marker dual validation; no single-gene labels."
            )
        if m.get("trajectory_verdict"):
            lines.append(
                f"{'轨迹' if zh else 'trajectory'}: {m.get('trajectory_verdict')} methods={m.get('trajectory_methods')}。"
            )
        elif "trajectory" in (plan.get("route") or []):
            lines.append("轨迹：计划评估连续性；未执行则无命运轴。" if zh else "Trajectory planned; no fate axis until executed.")
        return "\n".join(lines).strip() + "\n"
    if phase == "interpret":
        ip = state.get("interpretation_plan") or {}
        lines = [
            f"{'方法' if zh else 'method'}: {ip.get('method') or ('—' if zh else 'n/a')}。",
            f"{'基因集' if zh else 'gene sets'}: {ip.get('gene_sets') or ('—' if zh else 'n/a')}。",
        ]
        if ip.get("instructions"):
            lines.append(str(ip["instructions"]).strip())
        if ran and m.get("n_pathway_terms") is not None:
            lines.append(f"{'通路条目' if zh else 'pathway terms'}: {m.get('n_pathway_terms')}（{m.get('enrichment_engine') or 'n/a'}）。")
        elif not ran:
            lines.append(pending)
        arts = state.get("artifacts") or {}
        chains = arts.get("evidence_chains") or {}
        if chains.get("claims"):
            n_ok = chains.get("n_ok")
            n_all = chains.get("n_claims")
            lines.append(
                f"{'证据链' if zh else 'evidence chains'}: {n_ok}/{n_all} "
                + ("条完整（marker + 通路 p 值 + DOI）。" if zh else "complete (markers + pathway p + DOI).")
            )
            for c in (chains.get("claims") or [])[:5]:
                pw = c.get("pathway") or {}
                dois = ", ".join(
                    str((x or {}).get("doi") or "") for x in (c.get("citations") or []) if (x or {}).get("doi")
                )
                genes = ", ".join(
                    (m.get("gene") if isinstance(m, dict) else str(m)) for m in (c.get("markers") or [])[:4]
                )
                flag = "OK" if c.get("ok") else "INCOMPLETE"
                lines.append(
                    f"- [{flag}] {c.get('assertion')}: {genes}; {pw.get('id')} p={pw.get('pval')}; DOI {dois}"
                )
        return "\n".join(x for x in lines if x).strip() + "\n"
    return pending + "\n"


def build_dual(state: dict, lang: str | None = None) -> dict[str, Any]:
    use = lang or report_lang(state)
    if is_r_path(state):
        blocks = _r_blocks(state, lang=use)
    else:
        blocks = _python_blocks(state, lang=use)
    return {
        "format": "code-result-v1",
        "language": "r" if is_r_path(state) else "python",
        "query": state.get("user_query") or "",
        "blocks": blocks,
    }


def _python_blocks(state: dict, lang: str | None = None) -> list[dict[str, Any]]:
    blocks: list[dict[str, Any]] = []
    for phase, key, title in PHASES:
        code = strip_code_fences(state.get(key) or "")
        if not code:
            continue
        blocks.append(
            {
                "phase": phase,
                "title": title,
                "language": "python",
                "conclusion": phase_conclusion(state, phase, lang),
                "code": code,
                "script": {
                    "qc": "workspace/qc_preprocess.py",
                    "downstream": "workspace/cluster_annotate.py",
                    "cluster": "workspace/cluster_only.py",
                    "annotate": "workspace/annotate_deg.py",
                    "interpret": "workspace/interpret_pathways.py",
                }.get(phase),
            }
        )
    return blocks


def _r_blocks(state: dict, lang: str | None = None) -> list[dict[str, Any]]:
    from scagent.export_nb import seurat_phase_chunks

    zh = (lang or report_lang(state)) != "en"
    note = (
        "scAgent 不执行 R/Seurat kernel。下列代码可在 RStudio 或 IRkernel 中运行；本仓库现有 SOP 为 Scanpy。"
        if zh
        else "scAgent does not execute an R/Seurat kernel. Run the chunks in RStudio or IRkernel. Existing SOPs are Scanpy."
    )
    blocks: list[dict[str, Any]] = []
    for phase, title, code in seurat_phase_chunks(state):
        conclusion = phase_conclusion(state, phase, lang)
        if note not in conclusion:
            conclusion = note + "\n" + conclusion
        blocks.append(
            {
                "phase": phase,
                "title": title,
                "language": "r",
                "conclusion": conclusion,
                "code": strip_code_fences(code),
                "script": "outputs/analysis.Rmd",
            }
        )
    return blocks


def render_dual_markdown(
    state: dict | None = None,
    *,
    dual: dict | None = None,
    lang: str | None = None,
    heading: bool = True,
) -> str:
    use = lang or report_lang(state or {})
    dual = dual or build_dual(state or {}, lang=use)
    zh = use != "en"
    lines: list[str] = []
    if heading:
        lines += [
            "# " + ("代码-结果双重输出" if zh else "Code–result dual output"),
            "",
            (
                "每个阶段先给**结论**，再给无污染、可运行的代码。不要把叙事写进脚本。"
                if zh
                else "Each phase lists **conclusions** first, then unpolluted runnable code. Do not mix narrative into scripts."
            ),
            "",
        ]
    if not dual.get("blocks"):
        lines.append("（无代码。）" if zh else "(No code.)")
        return "\n".join(lines) + "\n"
    for b in dual["blocks"]:
        fence = "r" if b.get("language") == "r" else "python"
        lines += [
            f"## [结论] {b['title']}",
            "",
            (b.get("conclusion") or "").rstrip(),
            "",
            f"## [代码] {b['title']}",
            "",
        ]
        if b.get("script"):
            lines.append(f"`{b['script']}`")
            lines.append("")
        lines += [f"```{fence}", (b.get("code") or "").rstrip(), "```", ""]
    from scagent.publication_figures import render_publication_figure_inventory_markdown

    inv_md = render_publication_figure_inventory_markdown(state or {}, lang=use).rstrip()
    if inv_md:
        lines += [
            f"## {'发表级图表清单' if zh else 'Publication figure checklist'}",
            "",
            inv_md,
            "",
        ]
    return "\n".join(lines).rstrip() + "\n"


def render_dual_console(state: dict, *, preview_lines: int = 24) -> str:
    dual = build_dual(state)
    zh = report_lang(state) != "en"
    parts = ["======== " + ("代码-结果" if zh else "code–result") + " ========"]
    if not dual["blocks"]:
        parts.append("（无代码。）" if zh else "(No code.)")
        return "\n".join(parts)
    for b in dual["blocks"]:
        parts += [f"\n-------- [结论] {b['title']} --------", (b.get("conclusion") or "").rstrip()]
        code_lines = (b.get("code") or "").splitlines()
        head = "\n".join(code_lines[:preview_lines])
        more = f"\n… ({len(code_lines) - preview_lines} more lines)" if len(code_lines) > preview_lines else ""
        parts += [
            f"\n-------- [代码] {b['title']} → {b.get('script') or ''} --------",
            head + more,
        ]
    parts.append("\n完整双重输出: outputs/dual.md")
    parts.append("Notebook: outputs/analysis.ipynb" if dual.get("language") != "r" else "Rmd: outputs/analysis.Rmd")
    return "\n".join(parts) + "\n"


def export_dual(state: dict, out_dir: Path) -> Path:
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    path = out_dir / "dual.md"
    path.write_text(render_dual_markdown(state, lang=state.get("report_lang")), encoding="utf-8")
    return path
