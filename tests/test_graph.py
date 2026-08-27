from pathlib import Path

from langgraph.checkpoint.memory import MemorySaver

from sandbox.executor import write_and_maybe_run
from workflows.scRNA_langgraph import (
    after_hitl_mt,
    after_hitl_res,
    after_planner,
    after_review_down,
    after_review_qc,
    build_graph,
    run_analysis,
)


def test_graph_compiles():
    app = build_graph(checkpointer=MemorySaver())
    assert app is not None


def test_review_retry_routing():
    assert after_review_qc({"review_qc": {"passed": False}, "retry_count_qc": 0}) == "retry"
    assert after_review_qc({"review_qc": {"passed": True}, "retry_count_qc": 0}) == "hitl_res"
    assert after_review_qc({"review_qc": {"passed": False}, "retry_count_qc": 99}) == "review_pub"
    assert after_review_qc({"review_qc": {"passed": True}, "mode": "qc_only"}) == "review_pub"
    assert (
        after_review_qc(
            {
                "review_qc": {"passed": True},
                "execute_code": True,
                "execution_qc": {"executed": True, "ok": False},
                "retry_count_qc": 99,
            }
        )
        == "review_pub"
    )
    assert after_review_qc({"review_qc": {"passed": True}, "interrupt_after_qc": True, "auto_confirm": False}) == "hitl_res"
    assert after_planner({"mode": "annotate_only"}) == "hitl_res"
    assert after_planner({}) == "hitl_mt"
    assert after_hitl_mt({"interrupt_after_qc": True, "auto_confirm": False}) == "review_pub"
    assert after_hitl_mt({"interrupt_after_qc": True, "auto_confirm": False, "qc_choice": "recommended"}) == "qc_expert"
    assert after_hitl_res({"interrupt_after_qc": True, "auto_confirm": False}) == "review_pub"
    assert after_hitl_res({"interrupt_after_qc": True, "auto_confirm": False, "resolution_choice": "coarse"}) == "cluster_deg"
    assert after_review_down({"review_downstream": {"passed": False}, "retry_count_downstream": 0}) == "retry"
    assert after_review_down({"review_downstream": {"passed": True}, "retry_count_downstream": 0}) == "bio_interpret"
    assert after_review_down(
        {
            "review_downstream": {"passed": True},
            "execute_code": True,
            "execution_downstream": {"executed": True, "ok": False},
            "retry_count_downstream": 0,
        }
    ) == "retry"
    assert after_review_down(
        {
            "review_downstream": {"passed": True},
            "execute_code": True,
            "execution_downstream": {"executed": False, "ok": False, "jail": "schema"},
            "retry_count_downstream": 0,
        }
    ) == "retry"


def test_end_to_end_without_llm(tmp_path):
    state = run_analysis(
        "对 PBMC 做标准质控、聚类和注释",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        checkpointer=MemorySaver(),
    )
    assert "scanpy-scrna-seq" in (state.get("skills_used") or [])
    assert (state.get("plan") or {}).get("loop") == "plan-and-solve"
    assert (state.get("plan") or {}).get("collaboration") == "multi-agent"
    assert {a["id"] for a in (state.get("plan") or {}).get("agents") or []} >= {"qc_preprocess", "cluster_deg", "bio_interpret", "code_audit"}
    assert "gsea" in ((state.get("plan") or {}).get("route") or [])
    assert state.get("interpretation_plan", {}).get("role") == "bio_interpret"
    assert Path("workspace/interpret_pathways.py").exists()
    assert (state.get("plan") or {}).get("dag", {}).get("nodes")
    assert state["qc_strategy"]["plots_required"] == ["violin", "scatter", "mad"]
    assert "median_abs_deviation" in (state.get("code_qc") or "")
    assert 'side="high"' in state["code_qc"]
    assert "celltypist" in (state.get("code_downstream") or "").lower()
    assert state["review_qc"]["passed"] is True
    assert state["review_downstream"]["passed"] is True
    card = state.get("review_publication") or {}
    assert card.get("score") is not None
    assert 80 <= int(card["score"]) <= 95
    assert card.get("max_score") == 100
    assert "Overall score" in (state.get("report") or "")
    assert "✅ **QC:** PASS" in (state.get("report") or "")
    assert "证据链" in (state.get("report") or "")
    assert state["annotation_plan"]["forbid_single_gene"] is True
    assert state["annotation_plan"]["dual_validation"] is True
    assert "未执行" in state["report"] or "Not executed" in state["report"]
    assert Path("workspace/qc_preprocess.py").exists()
    assert Path("workspace/run_manifest.json").exists()
    mem = state.get("analysis_memory") or {}
    assert mem.get("qc", {}).get("method") == "mad"
    assert "user_query" not in mem
    assert "```yaml" in (state.get("report") or "")
    assert "analysis.ipynb" in (state.get("notebook") or "") or Path("outputs/analysis.ipynb").exists()
    assert Path("outputs/dual.md").is_file()
    dual_txt = Path("outputs/dual.md").read_text(encoding="utf-8")
    assert "## [结论]" in dual_txt
    assert "## [代码]" in dual_txt
    assert "median_abs_deviation" in dual_txt
    assert "[结论]" in (state.get("report") or "")


def test_r_language_degrades(tmp_path):
    state = run_analysis(
        "seurat 分析",
        data_path="",
        tissue="pbmc",
        language="r",
        execute_code=False,
        checkpointer=MemorySaver(),
    )
    assert state.get("r_degraded") is True
    assert state.get("status") == "r_degraded"
    assert not (state.get("code_qc") or "")
    assert "Seurat" in (state.get("report") or "") or "降级" in (state.get("report") or "")
    assert str(state.get("notebook") or "").endswith(".Rmd") or Path("outputs/analysis.Rmd").exists()


def test_interrupt_stops_before_mt(tmp_path):
    state = run_analysis(
        "QC 后暂停",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        interrupt_after_qc=True,
        auto_confirm=False,
        checkpointer=MemorySaver(),
    )
    assert state.get("status") == "awaiting_mt_confirmation"
    assert not (state.get("code_qc") or "")
    assert not (state.get("code_downstream") or "")
    assert Path("outputs/decisions/mt.html").is_file()
    assert 2 <= len((state.get("hitl_mt") or {}).get("options") or []) <= 3


def test_interrupt_stops_before_resolution(tmp_path):
    state = run_analysis(
        "MT 已确认，等 resolution",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        interrupt_after_qc=True,
        auto_confirm=False,
        qc_choice="recommended",
        checkpointer=MemorySaver(),
    )
    assert state.get("status") == "awaiting_resolution_confirmation"
    assert state.get("code_qc")
    assert not (state.get("code_downstream") or "")
    assert Path("outputs/decisions/resolution.html").is_file()
    assert 2 <= len((state.get("hitl_resolution") or {}).get("options") or []) <= 3


def test_executor_runs_python(tmp_path):
    r = write_and_maybe_run(
        "print('ok-sandbox')",
        workspace=tmp_path,
        execute=True,
        timeout=60,
        filename="analysis.py",
    )
    assert r["ok"] is True
    assert "ok-sandbox" in r["stdout"]
    assert (tmp_path / "analysis.py").exists()
    assert (tmp_path / "run_manifest.json").exists()
