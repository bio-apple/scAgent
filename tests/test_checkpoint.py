from langgraph.checkpoint.memory import MemorySaver

from workflows.checkpointing import SqlitePickleSaver
from workflows.scRNA_langgraph import _reuse_execution, build_graph, run_analysis


def test_graph_compiles_with_memory_checkpointer():
    app = build_graph(checkpointer=MemorySaver())
    assert app is not None


def test_checkpoint_persists_retry_and_execution_fields(tmp_path):
    saver = MemorySaver()
    tid = "unit-thread-1"
    state = run_analysis(
        "对 PBMC 做标准质控、聚类和注释",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        checkpointer=saver,
        thread_id=tid,
    )
    assert state.get("thread_id") == tid
    snap = build_graph(checkpointer=saver).get_state({"configurable": {"thread_id": tid}})
    vals = snap.values or {}
    assert "retry_count_qc" in vals
    assert vals.get("code_qc") or vals.get("plan")
    assert vals.get("execution_qc") is not None or vals.get("review_qc") is not None


def test_sqlite_pickle_saver_roundtrip(tmp_path):
    path = tmp_path / "ckpt.sqlite"
    s1 = SqlitePickleSaver(path)
    s1.storage["t"]["ns"]["cid"] = (("json", b"{}"), ("json", b"{}"), None)
    s1._flush()
    s2 = SqlitePickleSaver(path)
    assert "t" in s2.storage
    assert "ns" in s2.storage["t"]


def test_sqlite_pickle_graph_invoke(tmp_path):
    saver = SqlitePickleSaver(tmp_path / "ckpt.sqlite")
    tid = "sqlite-thread"
    state = run_analysis(
        "对 PBMC 做标准质控、聚类和注释",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        checkpointer=saver,
        thread_id=tid,
    )
    assert state.get("code_qc")
    snap = build_graph(checkpointer=saver).get_state({"configurable": {"thread_id": tid}})
    assert (snap.values or {}).get("retry_count_qc") is not None
    assert (snap.values or {}).get("execution_qc") is not None


def test_reuse_execution_requires_same_code():
    prev = {"ok": True, "executed": True, "code_fp": "deadbeef"}
    assert _reuse_execution(prev, "print(1)", False) is None
    from workflows.scRNA_langgraph import _code_fp

    prev["code_fp"] = _code_fp("print(1)")
    assert _reuse_execution(prev, "print(1)", True) is prev
    assert _reuse_execution({"ok": True, "executed": False, "code_fp": prev["code_fp"]}, "print(1)", True) is None
