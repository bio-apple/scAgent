from pathlib import Path

from langgraph.checkpoint.memory import MemorySaver

from sandbox.executor import write_and_maybe_run
from workflows.scRNA_langgraph import after_review_down, after_review_qc, build_graph, run_analysis


def test_graph_compiles():
    app = build_graph(checkpointer=MemorySaver())
    assert app is not None


def test_review_retry_routing():
    assert after_review_qc({"review_qc": {"passed": False}, "retry_count_qc": 0}) == "retry"
    assert after_review_qc({"review_qc": {"passed": True}, "retry_count_qc": 0}) == "annotation"
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
    assert after_review_qc({"review_qc": {"passed": True}, "interrupt_after_qc": True, "auto_confirm": False}) == "review_pub"
    assert after_review_down({"review_downstream": {"passed": False}, "retry_count_downstream": 0}) == "retry"
    assert after_review_down({"review_downstream": {"passed": True}, "retry_count_downstream": 0}) == "review_pub"
    assert (
        after_review_down(
            {
                "review_downstream": {"passed": True},
                "execute_code": True,
                "execution_downstream": {"executed": True, "ok": False},
                "retry_count_downstream": 0,
            }
        )
        == "retry"
    )


def test_end_to_end_without_llm(tmp_path):
    state = run_analysis(
        "对 PBMC 做标准质控、聚类和注释",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        checkpointer=MemorySaver(),
    )
    assert "scanpy-scrna-seq" in (state.get("skills_used") or [])
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
    assert state["annotation_plan"]["forbid_single_gene"] is True
    assert state["annotation_plan"]["dual_validation"] is True
    assert "未执行" in state["report"] or "Not executed" in state["report"]
    assert Path("workspace/qc_preprocess.py").exists()
    assert Path("workspace/run_manifest.json").exists()
    mem = state.get("analysis_memory") or {}
    assert mem.get("qc", {}).get("method") == "mad"
    assert "user_query" not in mem
    assert "```yaml" in (state.get("report") or "")


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


def test_interrupt_stops_after_qc(tmp_path):
    state = run_analysis(
        "QC 后暂停",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        interrupt_after_qc=True,
        auto_confirm=False,
        checkpointer=MemorySaver(),
    )
    assert state.get("status") == "awaiting_qc_confirmation"
    assert not (state.get("code_downstream") or "")


def test_executor_runs_python(tmp_path):
    r = write_and_maybe_run(
        "print('ok-sandbox')",
        workspace=tmp_path,
        execute=True,
        timeout=10,
        filename="analysis.py",
    )
    assert r["ok"] is True
    assert "ok-sandbox" in r["stdout"]
    assert (tmp_path / "analysis.py").exists()
    assert (tmp_path / "run_manifest.json").exists()
