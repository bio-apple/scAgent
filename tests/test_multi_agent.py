from pathlib import Path

from agents.roles import ROLES, assign_roles
from agents.templates import interpret_pathways_script, qc_preprocess_script
from scagent.enrich import ora, run_enrichment


def test_assign_roles_covers_four_agents():
    ids = {r["id"] for r in ROLES}
    assert ids == {"qc_preprocess", "cluster_deg", "bio_interpret", "code_audit"}
    assigned = assign_roles(["qc", "pca", "leiden", "gsea"])
    assert {a["id"] for a in assigned} == ids


def test_qc_script_includes_pca_not_deg():
    code = qc_preprocess_script({"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"}, {"nmads": 5})
    compile(code, "<qc_pca>", "exec")
    assert "scagent_pca" in code
    assert "rank_genes_groups" not in code


def test_ora_hits_interferon_set():
    genes = ["STAT1", "IRF1", "CXCL9", "CXCL10", "GBP1", "GAPDH"]
    rows = ora(genes, min_overlap=2)
    terms = [r["term"] for r in rows]
    assert "HALLMARK_INTERFERON_GAMMA_RESPONSE" in terms
    hit = next(r for r in rows if r["term"] == "HALLMARK_INTERFERON_GAMMA_RESPONSE")
    assert hit["overlap"] >= 4
    assert 0 <= hit["fdr"] <= 1
    out = run_enrichment(genes)
    assert out["n_terms"] >= 1
    assert out["engine"] in {"ora", "ora_gseapy_offline"}


def test_interpret_script_compiles(tmp_path):
    code = interpret_pathways_script({"tissue": "pbmc"}, {}, {"tissue": "pbmc"})
    compile(code, "<interpret>", "exec")
    assert "enrich_from_workspace" in code
    (tmp_path / "interpret_pathways.py").write_text(code, encoding="utf-8")
    assert (tmp_path / "interpret_pathways.py").is_file()
