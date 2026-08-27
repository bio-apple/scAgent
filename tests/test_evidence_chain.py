from agents.reviewer import publication_review
from agents.templates import interpret_pathways_script
from agents.writer import render_report
from scagent.evidence import assemble_claims, validate_claim, write_evidence_chains
from scagent.enrich import GO_SETS, ora


def test_tex_claim_requires_three_legs():
    incomplete = {
        "assertion": "Cluster 3 处于耗竭 T 细胞状态",
        "markers": [{"gene": "PDCD1", "source": "observed"}],
        "pathway": {},
        "citations": [],
    }
    bad = validate_claim(incomplete)
    assert bad["ok"] is False
    assert any("marker" in x.lower() or "基因" in x for x in bad["issues"])
    assert any("通路" in x or "p-value" in x.lower() for x in bad["issues"])
    assert any("DOI" in x or "PMID" in x for x in bad["issues"])

    complete = {
        "markers": [{"gene": "PDCD1"}, {"gene": "HAVCR2"}],
        "pathway": {"id": "GO:0002429", "pval": 1e-4, "fdr": 0.01},
        "citations": [{"pmid": "21739672", "doi": "10.1038/nri3156", "title": "T cell exhaustion"}],
    }
    ok = validate_claim(complete)
    assert ok["ok"] is True


def test_assemble_tex_chain_from_markers_and_ora():
    rows = [
        {
            "cluster": "3",
            "marker_label": "CD8 Tex",
            "fused": "CD8 Tex",
            "positive": ["CD8A", "PDCD1", "HAVCR2"],
            "negative": ["IL7R"],
            "dual_ok": True,
        }
    ]
    genes = ["PDCD1", "HAVCR2", "LAG3", "CD8A"]
    enrich = {"engine": "ora", "terms": ora(genes, gene_sets=GO_SETS, min_overlap=1)}
    payload = assemble_claims(rows, enrich)
    assert payload["n_claims"] == 1
    claim = payload["claims"][0]
    assert claim["ok"] is True, claim["issues"]
    genes_out = {m["gene"] for m in claim["markers"]}
    assert {"PDCD1", "HAVCR2"} <= genes_out
    assert str(claim["pathway"].get("id")).startswith("GO:0002429")
    assert claim["pathway"].get("pval") is not None
    dois = {c["doi"] for c in claim["citations"]}
    assert "10.1038/nri3156" in dois
    assert "21739672" in {c["pmid"] for c in claim["citations"]}


def test_unknown_state_cannot_invent_doi():
    payload = assemble_claims(
        [{"cluster": "0", "fused": "MysteriousBlob", "positive": ["GENEA", "GENEB"]}],
        {"terms": []},
    )
    assert payload["n_claims"] == 1
    assert payload["claims"][0]["ok"] is False
    assert payload["claims"][0]["citations"] == []


def test_write_evidence_chains_roundtrip(tmp_path):
    (tmp_path / "annotation_evidence.json").write_text(
        '[{"cluster":"3","fused":"CD8 Tex","positive":["PDCD1","HAVCR2","CD8A"]}]',
        encoding="utf-8",
    )
    (tmp_path / "pathway_enrichment.json").write_text(
        '{"terms":[{"term":"GO:0002429","pval":0.001,"fdr":0.02,"overlap":3,"method":"ora"}]}',
        encoding="utf-8",
    )
    out = write_evidence_chains(tmp_path)
    assert (tmp_path / "evidence_chains.json").is_file()
    assert out["all_ok"] is True
    assert out["claims"][0]["pathway"]["id"] == "GO:0002429"


def test_publication_review_fails_unsupported_state_claim():
    card = publication_review(
        {
            "execute_code": True,
            "execution_qc": {"executed": True, "ok": True},
            "execution_downstream": {"executed": True, "ok": True},
            "code_qc": "print(1)",
            "code_downstream": "positive negative fuse_annotation celltypist",
            "review_qc": {"passed": True},
            "review_downstream": {
                "passed": True,
                "has_dual": True,
                "has_fusion": True,
                "has_celltypist": True,
            },
            "metadata": {"n_samples": 1},
            "artifacts": {
                "figures": ["workspace/figures/violin.png", "workspace/figures/umap.png"],
                "evidence_chains": {
                    "n_claims": 1,
                    "n_ok": 0,
                    "claims": [
                        {
                            "ok": False,
                            "assertion": "Cluster 3 处于耗竭 T 细胞状态",
                            "issues": ["需要 PubMed PMID 或 DOI"],
                        }
                    ],
                },
            },
        }
    )
    by = {i["key"]: i for i in card["items"]}
    assert by["evidence"]["status"] == "fail"
    assert card["passed"] is False
    assert card["verdict"] == "FAIL"


def test_report_renders_evidence_chain():
    report = render_report(
        {
            "report_lang": "zh",
            "plan": {},
            "artifacts": {
                "evidence_chains": {
                    "n_claims": 1,
                    "n_ok": 1,
                    "claims": [
                        {
                            "ok": True,
                            "assertion": "Cluster 3 处于 CD8 Tex 状态",
                            "markers": [{"gene": "PDCD1", "source": "observed"}, {"gene": "HAVCR2", "source": "observed"}],
                            "pathway": {"id": "GO:0002429", "name": "immune response-activating cell surface receptor signaling pathway", "pval": 1e-4},
                            "citations": [{"doi": "10.1038/nri3156", "pmid": "21739672", "title": "T cell exhaustion"}],
                            "note": "支持性证据链，不是干预因果。",
                        }
                    ],
                }
            },
        }
    )
    assert "证据链" in report
    assert "PDCD1" in report
    assert "HAVCR2" in report
    assert "GO:0002429" in report
    assert "10.1038/nri3156" in report


def test_interpret_script_writes_evidence_chains():
    code = interpret_pathways_script({}, {}, {"tissue": "pbmc"})
    assert "write_evidence_chains" in code
    assert "evidence_chains" in code
