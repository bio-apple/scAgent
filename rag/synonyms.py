"""Bilingual query expansion so 批次效应校正 can hit Harmony / batch correction."""

from __future__ import annotations

# Longer phrases first.
_PHRASES: tuple[tuple[str, str], ...] = (
    ("批次效应校正", "batch effect correction harmony integration scvi"),
    ("批次校正", "batch correction harmony integration"),
    ("批次效应", "batch effect batch correction harmony"),
    ("差异表达", "differential expression deg pseudobulk deseq edger"),
    ("细胞注释", "cell type annotation celltypist marker"),
    ("双细胞", "doublet scrublet scdblfinder"),
    ("线粒体", "mitochondrial pct_mt qc mad"),
    ("质控", "quality control qc mad violin"),
    ("轨迹分析", "trajectory paga pseudotime"),
    ("整合", "integration harmony scvi scanorama batch"),
    ("聚类", "clustering leiden louvain"),
    ("harmony", "batch correction integration 批次"),
    ("scvi", "batch correction integration variational"),
    ("batch correction", "批次效应校正 harmony integration"),
    ("differential expression", "差异表达 pseudobulk deg"),
)


def expand_query(query: str) -> str:
    q = query or ""
    extra: list[str] = []
    low = q.lower()
    for src, dst in _PHRASES:
        if src.lower() in low:
            extra.append(dst)
    if not extra:
        return q
    return q + " " + " ".join(extra)
