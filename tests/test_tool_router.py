from scagent.tool_router import (
    PREFERENCE_TABLE,
    analysis_language,
    build_tool_route,
    format_route_table,
    resolve_module,
    router_cfg,
)


def test_preference_table_r_first_python_backup():
    expected = {
        "qc": ("seurat", "scanpy"),
        "clustering": ("seurat", "scanpy"),
        "integration": ("harmony_r", "scvi"),
        "annotation": ("azimuth", "celltypist"),
        "cellchat": ("cellchat", "squidpy"),
        "spatial": ("giotto", "squidpy"),
    }
    got = {mod: (r_tool, py_tool) for mod, r_tool, py_tool in PREFERENCE_TABLE}
    for mod, pair in expected.items():
        assert got[mod] == pair


def test_router_defaults_r_first():
    cfg = {"tool_router": {"policy": "r_first"}, "analysis": {"language": "r_first"}}
    r = resolve_module("qc", cfg=cfg)
    assert r["primary_tool"] == "seurat"
    assert r["fallback_tool"] == "scanpy"


def test_integration_python_backup_is_scvi():
    cfg = {"tool_router": {"policy": "r_first"}, "analysis": {"language": "r_first"}}
    r = resolve_module("integration", cfg=cfg, meta={"need_batch_correction": True})
    assert r["primary_tool"] == "harmony_r"
    assert r["fallback_tool"] == "scvi"


def test_annotation_python_backup_is_celltypist():
    cfg = {"tool_router": {"policy": "python_only"}, "analysis": {"language": "python"}}
    r = resolve_module("annotation", cfg=cfg)
    assert r["engine"] == "python"
    assert r["tool"] == "celltypist"
    assert r["primary_tool"] == "azimuth"


def test_llm_does_not_set_language():
    """plan.language / LLM fields must not override the router."""
    cfg = {"tool_router": {"policy": "r_first"}, "analysis": {"language": "r_first"}}
    r = resolve_module(
        "qc",
        cfg=cfg,
        plan={"language": "python", "narrative": "please use scanpy"},
    )
    assert r["primary_tool"] == "seurat"
    assert r["policy"] == "r_first"


def test_python_only_forces_scanpy():
    cfg = {"tool_router": {"policy": "python_only"}, "analysis": {"language": "python"}}
    r = resolve_module("qc", cfg=cfg)
    assert r["engine"] == "python"
    assert r["tool"] == "scanpy"


def test_build_tool_route_table():
    route = build_tool_route({"tissue": "pbmc", "need_batch_correction": True}, {"route": []}, cfg={"tool_router": {"policy": "r_first"}, "analysis": {"language": "r_first"}})
    assert route["policy"] == "r_first"
    assert "qc" in route["routes"]
    assert route["system_prompt"].startswith("Always use R")


def test_format_route_table_zh():
    route = build_tool_route({}, {}, cfg={"tool_router": {"policy": "r_first"}, "analysis": {"language": "r_first"}})
    txt = format_route_table(route, lang="zh")
    assert "Tool Router" in txt
    assert "seurat" in txt.lower() or "scanpy" in txt.lower()


def test_legacy_language_r_is_r_only():
    cfg = {"analysis": {"language": "r"}}
    assert analysis_language(cfg) == "r"
    r = resolve_module("qc", cfg=cfg)
    assert r["policy"] == "r_only"


def test_templates_include_router_bootstrap():
    from agents.templates import cluster_annotate_script, qc_preprocess_script

    qc = qc_preprocess_script({"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"}, {"nmads": 5})
    assert "maybe_run_r_qc" in qc
    dn = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5},
        {"integrator": None},
    )
    assert "maybe_run_r_downstream" in dn
