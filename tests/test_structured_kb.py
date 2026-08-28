from agents.markers import load_marker_catalog
from scagent.enrich import default_gene_sets
from scagent.evidence import CATALOG_PATH, match_cell_state
from scagent.kb import gene_sets_from_kb, load_structured, lookup_structured


def test_structured_dirs_have_records():
    recs = load_structured()
    cols = {r["collection"] for r in recs}
    assert {"cell_ontology", "marker_db", "pathway", "disease_signature", "tissue_reference"} <= cols
    assert any(r.get("id") == "CL:0000084" for r in recs)


def test_lookup_t_cell_returns_ontology_and_markers():
    hits = lookup_structured("T cell CD3D", collections=["cell_ontology", "marker_db"], tissue="pbmc", top_k=8)
    blob = " ".join(h["text"] for h in hits)
    assert "CL:0000084" in blob
    assert "CD3D" in blob
    assert hits[0].get("retrieval") == "structured"


def test_lookup_hallmark_and_tex_signature():
    pw = lookup_structured("HYPOXIA VEGFA", collections=["pathway"], top_k=3)
    assert any("HALLMARK_HYPOXIA" in (h.get("text") or "") for h in pw)
    st = lookup_structured("exhausted T PDCD1", collections=["disease_signature"], top_k=3)
    assert any("PDCD1" in (h.get("text") or "") and "HAVCR2" in (h.get("text") or "") for h in st)


def test_marker_catalog_reads_marker_db():
    cat = load_marker_catalog(tissue="pbmc")
    tex = next(t for t in cat["cell_types"] if t["name"] == "CD8 Tex")
    assert tex["positive"][:2] == ["CD8A", "PDCD1"]
    tcell = next(t for t in cat["cell_types"] if t["name"] == "T cell")
    assert tcell.get("cl_id") == "CL:0000084"


def test_gene_sets_and_evidence_use_kb():
    sets = gene_sets_from_kb() or default_gene_sets()
    assert "HALLMARK_INTERFERON_GAMMA_RESPONSE" in sets
    assert "GO:0002429" in sets
    assert CATALOG_PATH.as_posix().endswith("disease_signature/cell_states.json")
    assert match_cell_state("exhausted T")["name"] == "CD8 Tex"


def test_tissue_reference_heart():
    hits = lookup_structured("heart cardiomyocyte", collections=["tissue_reference"], tissue="heart", top_k=3)
    blob = " ".join(h["text"] for h in hits).lower()
    assert "cardiomyocyte" in blob
    assert "mito" in blob or "tnnt2" in blob or "heart" in blob
