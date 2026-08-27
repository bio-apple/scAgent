from agents.dependencies import resolve_route
from agents.intent import parse_intent, rule_intent
from scagent.analysis import needs_condition_de


def test_degradation_is_not_deg():
    assert rule_intent("degradation analysis")["condition_comparison"] is False
    assert "deg" not in rule_intent("RNA degradation QC")["intents"]
    assert needs_condition_de("degradation analysis") is False
    assert needs_condition_de("比较对照组 vs 处理组的差异表达") is True
    assert needs_condition_de("对 PBMC 做标准注释") is False


def test_standard_query_intents():
    intent = parse_intent("对 PBMC 做标准质控、聚类和注释")
    assert "qc" in intent["intents"]
    assert "annotation" in intent["intents"]
    assert "enrichment" in intent["intents"]
    assert intent["condition_comparison"] is False


def test_deg_route_requires_annotation_groupby():
    route = resolve_route(["deg"])
    assert "qc" in route
    assert "annotate" in route
    assert "pseudobulk_deg" in route
    assert route.index("qc") < route.index("annotate") < route.index("pseudobulk_deg")
    assert route.index("pca") < route.index("leiden")
    assert route.index("umap") < route.index("cluster_deg")
    assert route.index("leiden") < route.index("cluster_deg")


def test_harmony_before_neighbors():
    route = resolve_route(["clustering", "annotation"], integrator="harmony")
    assert "harmony" in route
    assert route.index("harmony") < route.index("neighbors")
    assert route.index("neighbors") < route.index("leiden")
    assert route.index("umap") < route.index("annotate")


def test_trajectory_requires_umap_and_leiden():
    route = resolve_route(["trajectory"])
    assert route.index("pca") < route.index("neighbors")
    assert route.index("neighbors") < route.index("umap")
    assert route.index("leiden") < route.index("trajectory")
    assert route.index("umap") < route.index("trajectory")
    from agents.dependencies import serialize_dag

    dag = serialize_dag(route)
    assert dag["loop"] == "plan-and-solve"
    assert any(e["to"] == "trajectory" for e in dag["edges"])
