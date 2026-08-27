from rag.ingest import ingest
from rag.retriever import _bm25_bundle, retrieve


def test_rag_papers_harmony():
    ingest()
    _bm25_bundle.cache_clear()
    hits = retrieve("Harmony batch correction PCA embedding", collection="papers", top_k=8)
    assert hits, "expected BM25 hits from knowledge/papers"
    blob = " ".join(h["source"] + " " + h["text"] for h in hits).lower()
    assert "harmony" in blob
    assert "korsunsky" in blob or "pca" in blob


def test_rag_best_practices_qc():
    ingest()
    _bm25_bundle.cache_clear()
    hits = retrieve("MAD-based filtering per sample mitochondrial", collection="best_practices", top_k=5)
    assert hits, "expected BM25 hits from best_practices/reference"
    blob = " ".join(h["source"] + " " + h["text"] for h in hits).lower()
    assert "mad" in blob
    assert "qc.md" in blob or "quality control" in blob


def test_rag_mito_qc_chinese():
    ingest()
    _bm25_bundle.cache_clear()
    hits = retrieve("线粒体 QC MAD 阈值 肿瘤", collection="papers", top_k=5)
    assert hits
    blob = " ".join(h["text"] + h["source"] for h in hits)
    assert "Yates" in blob or "线粒体" in blob
