from agents.planner import choose_integrator
from agents.reviewer import audit_code, audit_execution, format_review_card, publication_review
from agents.templates import LOCKED_START, cluster_annotate_script, qc_preprocess_script, scanpy_script
from scagent.inspect_data import detect_platform, detect_species_from_genes, gene_composition, inspect_data
from scagent.skills_loader import list_skills, recommend_skills


def test_existing_skills_preserved():
    names = {s.name for s in list_skills()}
    expected = {
        "anndata-data-structure",
        "celltypist-cell-annotation",
        "cellxgene-census",
        "harmony-batch-correction",
        "popv-cell-annotation",
        "scanpy-scrna-seq",
        "scvi-tools-single-cell",
        "single-cell-annotation-guide",
        "single-cell-annotation",
        "cellchat-cell-communication",
    }
    assert expected <= names


def test_recommend_harmony_for_multi_sample():
    skills = recommend_skills({"n_samples": 4, "need_batch_correction": True, "tissue": "pbmc"})
    assert "scanpy-scrna-seq" in skills
    assert "harmony-batch-correction" in skills
    assert "celltypist-cell-annotation" in skills
    assert "single-cell-annotation" in skills


def test_recommend_cellchat_for_communication_query():
    skills = recommend_skills({"tissue": "tumor", "task": "CellChat 分析肿瘤与 T 细胞的配体受体通讯"})
    assert "cellchat-cell-communication" in skills
    assert "scanpy-scrna-seq" in skills


def test_choose_integrator():
    assert choose_integrator({"n_samples": 1, "need_batch_correction": False}) is None
    assert choose_integrator({"n_samples": 3, "n_cells": 8000, "need_batch_correction": True}) == "harmony"
    assert choose_integrator({"n_samples": 10, "n_cells": 1000, "need_batch_correction": True}) == "scvi"
    assert choose_integrator({"n_samples": 2, "n_cells": 200_000, "need_batch_correction": True}) == "scvi"


def test_inspect_parse_platform(tmp_path):
    p = tmp_path / "parse_biosciences_sample.h5ad"
    p.write_text("not-an-h5ad")
    meta = inspect_data(str(p), tissue="PBMC")
    assert meta["platform"] == "parse"
    assert meta["tissue"] == "pbmc"


def test_detect_species_and_platform():
    assert detect_species_from_genes(["MT-ND1", "GAPDH", "MT-CO1"]) == "human"
    assert detect_species_from_genes(["mt-Nd1", "Gapdh"]) == "mouse"
    assert detect_platform(__import__("pathlib").Path("filtered_feature_bc_matrix.h5")) == "10x"


def test_audit_requires_qc_trio():
    bad = "import scanpy as sc\nsc.pp.filter_cells(adata, min_genes=200)\n"
    r = audit_code(bad, {"need_batch_correction": True}, phase="qc")
    assert r["passed"] is False
    assert r["has_violin"] is False
    assert r["has_mad"] is False


def test_qc_template_passes_reviewer():
    meta = {"data_path": "x.h5ad", "species": "human", "need_batch_correction": True, "sample_key": "batch", "tissue": "pbmc"}
    code = qc_preprocess_script(meta, {"nmads": 5, "extra_qc": ["hb"]})
    compile(code, "<qc>", "exec")
    assert LOCKED_START in code
    assert "log1p=True" in code
    assert 'side="high"' in code
    assert "detect_doublets" in code
    assert "scrublet" in code.lower()
    assert "predicted_doublet" in code
    assert "select_hvg" in code
    assert "filter_genes" in code
    r = audit_code(code, meta, phase="qc")
    assert r["has_violin"] and r["has_scatter"] and r["has_mad"] and r["has_locked"]
    assert r["passed"] is True


def test_downstream_template_dual_validation():
    meta = {"data_path": "x.h5ad", "species": "human", "need_batch_correction": True, "sample_key": "batch", "tissue": "pbmc"}
    plan = {"integrator": "harmony"}
    code = cluster_annotate_script(meta, {"nmads": 5}, plan)
    compile(code, "<down>", "exec")
    assert "celltypist" in code.lower()
    assert "ref2_label" in code
    assert "fuse_annotation" in code
    assert "positive" in code and "negative" in code
    assert "harmonypy" in code
    assert "use_raw" in code
    assert "use_highly_variable=True" in code
    r = audit_code(code, meta, phase="downstream")
    assert r["has_celltypist"] and r["has_dual"]
    assert r.get("has_ref2") is True
    assert r.get("has_fusion") is True
    assert r["passed"] is True


def test_combined_script_compiles():
    code = scanpy_script(
        {"data_path": "x.h5ad", "species": "human", "need_batch_correction": True, "sample_key": "batch"},
        {"nmads": 5},
        {"integrator": "harmony"},
    )
    compile(code, "<scanpy_script>", "exec")
    assert "harmony" in code.lower()


def test_overfilter_respects_config_threshold():
    ok = audit_execution(
        {"executed": True, "ok": True},
        {"metrics": {"pct_removed": 40}, "h5ads": {"qc": "/tmp/adata_qc.h5ad"}, "figures": ["/tmp/violin.png", "/tmp/scatter.png"]},
        phase="qc",
        execute_code=True,
        metadata={"overfilter_warn_pct": 80},
    )
    assert ok["passed"] is True
    assert not any("过度过滤" in x for x in ok["issues"])


def test_audit_code_issue_records():
    r = audit_code("import scanpy as sc\n", {}, phase="qc")
    assert r["passed"] is False
    assert r.get("issue_records")
    assert all("id" in rec and "message" in rec for rec in r["issue_records"])
    assert any(rec["id"] == "qc.violin" for rec in r["issue_records"])
    skip = audit_execution({}, {}, phase="qc", execute_code=False)
    assert skip["passed"] is True and skip["skipped"] is True
    fail = audit_execution(
        {"executed": True, "ok": False, "stderr": "Traceback"},
        {"metrics": {}, "h5ads": {}, "figures": []},
        phase="qc",
        execute_code=True,
    )
    assert fail["passed"] is False
    over = audit_execution(
        {"executed": True, "ok": True},
        {"metrics": {"pct_removed": 55}, "h5ads": {"qc": "/tmp/adata_qc.h5ad"}, "figures": ["/tmp/violin.png", "/tmp/scatter.png"]},
        phase="qc",
        execute_code=True,
    )
    assert over["passed"] is False
    assert any("过度过滤" in x for x in over["issues"])


def test_audit_rejects_umap_mix_and_raw_p():
    mix = (
        "import scanpy as sc\nnp.random.seed(0)\n"
        "sc.external.pp.harmony_integrate(adata, 'batch')\n"
        "sc.tl.umap(adata)\n"
        "print('UMAP looks mixed so integration succeeded')\n"
        "celltypist.annotate(adata)\n"
        "ref2_label = adata.obs['leiden']\n"
        "positive = ['MS4A1', 'CD79A']\n"
        "negative = ['CD3D']\n"
        "cell_type_l1 = 'B'\n"
    )
    r = audit_code(mix, {"need_batch_correction": True, "tissue": "pbmc", "species": "human"}, phase="downstream")
    assert r["passed"] is False
    assert any(rec["id"] == "down.umap_mix" for rec in r["issue_records"])
    deg = (
        "import scanpy as sc\nnp.random.seed(0)\n"
        "celltypist.annotate(adata)\n"
        "ref2_label = adata.obs['leiden']\n"
        "positive = ['MS4A1', 'CD79A']\n"
        "negative = ['CD3D']\n"
        "cell_type_l1 = 'B'\n"
        "sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')\n"
        "print(adata.uns['rank_genes_groups']['pvals'])\n"
    )
    r2 = audit_code(deg, {"tissue": "pbmc", "species": "human"}, phase="downstream")
    assert r2["passed"] is False
    assert any(rec["id"] in {"down.padj", "down.deg_note"} for rec in r2["issue_records"])


def test_audit_rejects_wilcoxon_on_scaled_x():
    scaled = (
        "import scanpy as sc\nnp.random.seed(0)\n"
        "sc.pp.scale(adata, max_value=10)\n"
        "celltypist.annotate(adata)\n"
        "ref2_label = adata.obs['leiden']\n"
        "positive = ['MS4A1', 'CD79A']\n"
        "negative = ['CD3D']\n"
        "cell_type_l1 = 'B'\n"
        "sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')\n"
        "print('pvals_adj and exploratory Wilcoxon; use pseudobulk + FDR for condition DE')\n"
    )
    r = audit_code(scaled, {"tissue": "pbmc", "species": "human"}, phase="downstream")
    assert r["passed"] is False
    assert any(rec["id"] == "down.deg_scaled" for rec in r["issue_records"])


def test_gene_composition_flags_hb():
    c = gene_composition(["MT-ND1", "RPS3", "RPL5", "HBA1", "HBB", "HBA2", "GAPDH", "MKI67", "CCNB1"])
    assert c["n_mt_genes"] >= 1
    assert c["n_ribo_genes"] >= 2
    assert c["need_hb_qc"] is True
    assert c["need_cell_cycle"] is True
    none = gene_composition(["GAPDH", "ACTB"])
    assert none["need_hb_qc"] is False
    assert none["need_cell_cycle"] is False


def test_publication_review_card_and_score():
    meta = {
        "data_path": "x.h5ad",
        "species": "human",
        "need_batch_correction": False,
        "n_samples": 1,
        "tissue": "pbmc",
    }
    qc = qc_preprocess_script(meta, {"nmads": 5})
    down = cluster_annotate_script(meta, {"nmads": 5}, {"integrator": None})
    rq = audit_code(qc, meta, phase="qc")
    rd = audit_code(down, meta, phase="downstream")
    card = publication_review(
        {
            "metadata": meta,
            "plan": {"integrator": None},
            "code_qc": qc,
            "code_downstream": down,
            "review_qc": rq,
            "review_downstream": rd,
            "execute_code": False,
            "artifacts": {},
        }
    )
    by_key = {i["key"]: i for i in card["items"]}
    assert by_key["qc"]["status"] == "pass"
    assert by_key["batch_correction"]["status"] == "pass"
    assert by_key["doublet_detection"]["status"] == "pass"
    assert by_key["markers"]["status"] == "pass"
    assert by_key["deg"]["status"] == "pass"
    assert by_key["figures"]["status"] == "missing"
    assert by_key["annotation"]["status"] == "pass"
    assert by_key["evidence"]["status"] == "missing"
    assert 80 <= card["score"] <= 95
    text = format_review_card(card, "zh")
    assert "✅ **QC:** PASS" in text
    assert "⚠️ **Figures:** Missing" in text
    assert "Overall score:" in text
    assert "/ 100" in text
