from rag.ingest import chunk_semantic, ingest
from rag.retriever import clear_retrieve_cache, retrieve
from rag.synonyms import expand_query


def _refresh():
    ingest()
    clear_retrieve_cache()


def test_rag_papers_harmony():
    _refresh()
    hits = retrieve("Harmony batch correction PCA embedding", collection="papers", top_k=8)
    assert hits, "expected BM25 hits from knowledge/papers"
    blob = " ".join(h["source"] + " " + h["text"] for h in hits).lower()
    assert "harmony" in blob
    assert "korsunsky" in blob or "pca" in blob


def test_rag_best_practices_qc():
    _refresh()
    hits = retrieve("MAD-based filtering per sample mitochondrial", collection="best_practices", top_k=5)
    assert hits, "expected BM25 hits from best_practices/reference"
    blob = " ".join(h["source"] + " " + h["text"] for h in hits).lower()
    assert "mad" in blob
    assert "qc.md" in blob or "quality control" in blob


def test_rag_mito_qc_chinese():
    _refresh()
    hits = retrieve("线粒体 QC MAD 阈值 肿瘤", collection="papers", top_k=5)
    assert hits
    blob = " ".join(h["text"] + h["source"] for h in hits)
    assert "Yates" in blob or "线粒体" in blob


def test_hybrid_chinese_batch_hits_harmony():
    assert "harmony" in expand_query("批次效应校正").lower()
    assert "palantir" in expand_query("轨迹分析").lower()
    assert "scvelo" in expand_query("rna velocity").lower()
    _refresh()
    hits = retrieve("批次效应校正", collection="papers", top_k=8)
    assert hits
    blob = " ".join(h["source"] + " " + h["text"] for h in hits).lower()
    assert "harmony" in blob or "batch" in blob or "korsunsky" in blob
    assert hits[0].get("retrieval") == "hybrid"


def test_semantic_chunk_keeps_paragraphs():
    para = "完整段落不被固定窗口切断。" * 3
    text = f"## 第一节\n\n{para}\n\n## 第二节\n\n第二段内容。"
    chunks = chunk_semantic(text, chunk_size=200, overlap=20)
    assert chunks
    assert any(para in c for c in chunks)
    assert any("第一节" in c or "第二节" in c for c in chunks)
