"""Named specialist agents. Planner assigns work; no single generalist owns the run."""

from __future__ import annotations

ROLES: tuple[dict[str, str], ...] = (
    {
        "id": "qc_preprocess",
        "name": "QC & Preprocessing Agent",
        "charge": "数据校验、MAD/双细胞/ambient 清洗、HVG 与 PCA 降维",
    },
    {
        "id": "cluster_deg",
        "name": "Clustering & Differential Agent",
        "charge": "邻居图/UMAP/Leiden、注释证据融合、探索性 marker、组间 DEG、轨迹/命运（Palantir/DPT/scVelo）",
    },
    {
        "id": "bio_interpret",
        "name": "Biological Interpretation Agent",
        "charge": "通路富集（GSEA/GSVA 或 ORA）与文献验证",
    },
    {
        "id": "code_audit",
        "name": "Code Audit & Execution Agent",
        "charge": "把各专家指令变成可执行 Python/R、schema 拦截、沙箱运行、失败自修复",
    },
)

# DAG node -> domain agent (code_audit executes every phase)
_OWNER = {
    "qc": "qc_preprocess",
    "normalize": "qc_preprocess",
    "hvg": "qc_preprocess",
    "pca": "qc_preprocess",
    "ambient_soupx": "qc_preprocess",
    "ambient_decontx": "qc_preprocess",
    "impute_magic": "qc_preprocess",
    "impute_alra": "qc_preprocess",
    "harmony": "cluster_deg",
    "scvi": "cluster_deg",
    "cca": "cluster_deg",
    "bbknn": "cluster_deg",
    "neighbors": "cluster_deg",
    "leiden": "cluster_deg",
    "umap": "cluster_deg",
    "annotate": "cluster_deg",
    "cluster_deg": "cluster_deg",
    "pseudobulk_deg": "cluster_deg",
    "trajectory": "cluster_deg",
    "gsea": "bio_interpret",
}


def assign_roles(route: list[str] | None) -> list[dict[str, str]]:
    """Subset of ROLES that this DAG actually needs, plus the shared code_audit worker."""
    needed = {"code_audit"}
    for step in route or []:
        owner = _OWNER.get(step)
        if owner:
            needed.add(owner)
    return [dict(r) for r in ROLES if r["id"] in needed]
