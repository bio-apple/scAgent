from agents.memory import build_memory, dump_memory_yaml, persist_memory


def test_memory_records_steps_not_chat():
    state = {
        "data_path": "/data/PBMC001.h5ad",
        "thread_id": "t-pbmc",
        "metadata": {"tissue": "pbmc", "expression_layer": "counts"},
        "qc_strategy": {"method": "mad", "nmads": 5, "doublets": True, "ambient": "none", "hard": {}},
        "plan": {"integrator": "harmony", "celltypist_model": "Immune_All_Low.pkl", "needs_pseudobulk": False},
        "annotation_plan": {"dual_validation": True},
        "code_qc": "print('qc')",
        "code_downstream": "celltypist.annotate(adata)\npositive = ['MS4A1']\nnegative = ['CD3D']\nref2_label = 'x'\nfuse_annotation(adata)\nsc.tl.rank_genes_groups(adata, 'leiden')\n",
        "execution_qc": {"executed": True, "ok": True, "snapshots": [".cache/steps/qc/adata_qc.h5ad"]},
        "execution_downstream": {"executed": True, "ok": False},
        "artifacts": {
            "metrics": {"n_before": 1000, "n_after": 800, "pct_removed": 20.0},
            "h5ads": {"qc": "workspace/adata_qc.h5ad"},
        },
        "logs": ["this is chat-like log and must not appear as a top-level memory key"],
        "user_query": "帮我分析一下",
    }
    mem = build_memory(state)
    assert mem["sample"] == "PBMC001"
    assert mem["qc"]["method"] == "mad"
    assert mem["qc"]["mt"] is None
    assert mem["qc"]["n_before"] == 1000
    assert mem["normalize"] == "normalize_total+log1p"
    assert mem["integration"] == "harmony"
    assert mem["annotation"] == ["CellTypist", "Marker", "ref2", "fusion"]
    assert mem["deg"] == "wilcox (exploratory)"
    assert mem["steps"]["qc"]["status"] == "ok"
    assert mem["steps"]["downstream"]["status"] == "failed"
    assert mem["resume_from"] == "downstream"
    text = dump_memory_yaml(mem)
    assert "帮我分析一下" not in text
    assert "chat-like" not in text
    assert "qc:" in text
    assert "integration: harmony" in text


def test_persist_memory_yaml(tmp_path):
    from scagent.config import load_config

    cfg = dict(load_config())
    cfg["paths"] = {**cfg["paths"], "cache": str(tmp_path / "cache")}
    state = {
        "data_path": "x.h5ad",
        "thread_id": "t1",
        "metadata": {"tissue": "heart"},
        "qc_strategy": {"method": "mad", "nmads": 6, "hard": {}},
        "plan": {"integrator": None},
        "code_qc": "x",
        "execution_qc": {"executed": False},
    }
    mem = persist_memory(state, extra_dir=tmp_path / "out", cfg=cfg)
    assert (tmp_path / "out" / "memory.yaml").is_file()
    assert (tmp_path / "cache" / "memory.yaml").is_file()
    assert mem["tissue"] == "heart"
    assert mem["qc"]["nmads"] == 6
    assert mem["steps"]["qc"]["status"] == "planned"
