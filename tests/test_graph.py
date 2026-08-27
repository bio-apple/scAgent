from sandbox.executor import write_and_maybe_run
from workflows.scRNA_langgraph import after_review, build_graph, run_analysis


def test_graph_compiles():
    app = build_graph()
    assert app is not None


def test_review_retry_routing():
    assert after_review({"review": {"passed": False}, "retry_count": 0}) == "retry"
    assert after_review({"review": {"passed": True}, "retry_count": 0}) == "annotate"
    assert after_review({"review": {"passed": False}, "retry_count": 99}) == "annotate"


def test_end_to_end_without_llm(tmp_path):
    state = run_analysis(
        "对 PBMC 做标准质控、聚类和注释",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
    )
    assert "scanpy-scrna-seq" in (state.get("skills_used") or [])
    assert state["qc_strategy"]["plots_required"] == ["violin", "scatter", "mad"]
    assert "median_abs_deviation" in state["code"]
    assert state["review"]["passed"] is True
    assert state["annotation_plan"]["forbid_single_gene"] is True
    assert "质控" in state["report"] or "QC" in state["report"] or "分析" in state["report"]


def test_executor_runs_python(tmp_path):
    r = write_and_maybe_run(
        "print('ok-sandbox')",
        workspace=tmp_path,
        execute=True,
        timeout=10,
    )
    assert r["ok"] is True
    assert "ok-sandbox" in r["stdout"]
    assert (tmp_path / "analysis.py").exists()
