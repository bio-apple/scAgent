from agents.literature import (
    fetch_phase_literature,
    format_literature_report_block,
    literature_recommendations,
)


def test_literature_recommendations_format():
    hits = [
        {
            "stem": "2020-batch-benchmark",
            "section": "methods",
            "text": "Harmony is recommended for batch integration on large datasets. Runtime is shorter than LIGER.",
        }
    ]
    recs = literature_recommendations(hits)
    assert recs
    assert "2020-batch-benchmark" in recs[0]
    assert "methods" in recs[0]
    assert "Harmony" in recs[0]


def test_fetch_phase_literature_qc(monkeypatch):
    fake = [
        {
            "stem": "mito-qc",
            "section": "results",
            "source": "knowledge/papers/.parsed/mito.md",
            "text": "A uniform 5% mtDNA threshold is not valid across tissues.",
            "score": 1.0,
        }
    ]

    def _fake_search(query, **_k):
        assert "mitochondrial" in query.lower() or "quality" in query.lower() or "QC" in query
        return fake

    monkeypatch.setattr("agents.literature.search_paper_knowledge", _fake_search)
    lit = fetch_phase_literature("qc", tissue="heart", user_query="filter dying cells")
    assert lit["paper_recs"]
    assert "mito-qc" in lit["paper_excerpt"]
    assert "5%" in lit["paper_recs"][0]


def test_format_literature_report_block_includes_phases():
    md = format_literature_report_block(
        plan={"paper_recs": ["paper-a [abstract]: Use Harmony first."]},
        qc={"paper_recs": ["paper-b [methods]: Prefer MAD over fixed mito%."]},
        ann={},
        interpret={},
        lang="zh",
    )
    assert "## 文献最佳实践建议" in md
    assert "### Planner" in md
    assert "### QC" in md
    assert "Harmony" in md
    assert "MAD" in md
