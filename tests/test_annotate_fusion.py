import pandas as pd

from agents.reviewer import audit_code
from agents.templates import cluster_annotate_script
from scagent.annotate import deg_catalog_labels, fuse_annotation, labels_agree


class _A:
    def __init__(self, obs, uns=None):
        self.obs = obs
        self.n_obs = len(obs)
        self.uns = uns or {}


def test_labels_agree_alias():
    assert labels_agree("CD8 T", "CD8 T cells") is True
    assert labels_agree("B cell", "NK") is False
    assert labels_agree("unknown", "B cell") is False


def test_fuse_majority_and_conflict():
    obs = pd.DataFrame(
        {
            "marker_label": ["CD8 T", "CD8 T", "unknown", "B cell"],
            "celltypist_label": ["CD8 T cells", "B cell", "NK", "B cell"],
            "deg_label": ["CD8 T", "B cell", "NK", "unknown"],
        }
    )
    adata = _A(obs)
    fuse_annotation(adata)
    assert adata.obs.loc[0, "cell_type"] == "CD8 T"
    assert adata.obs.loc[0, "annotation_status"] == "fused"
    assert adata.obs.loc[1, "cell_type"] == "mixed"
    assert adata.obs.loc[1, "annotation_status"] == "conflict_mixed"
    assert str(adata.obs.loc[2, "cell_type"]).endswith("|unvalidated")
    assert adata.obs.loc[3, "annotation_status"] == "fused"


def test_deg_catalog_overlap():
    obs = pd.DataFrame({"leiden": ["0", "0", "1", "1"]})
    uns = {
        "rank_genes_groups": {
            "names": pd.DataFrame(
                {
                    "0": ["MS4A1", "CD79A", "CD19", "X"],
                    "1": ["CD3D", "CD8A", "NKG7", "Y"],
                }
            )
        }
    }
    catalog = {
        "cell_types": [
            {"name": "B cell", "positive": ["MS4A1", "CD79A", "CD19"]},
            {"name": "T cell", "positive": ["CD3D", "CD8A", "IL7R"]},
        ]
    }
    adata = _A(obs, uns)
    deg_catalog_labels(adata, catalog)
    assert list(adata.obs["deg_label"]) == ["B cell", "B cell", "T cell", "T cell"]


def test_template_does_not_double_count_deg_as_ref2():
    code = cluster_annotate_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"nmads": 5},
        {"integrator": None},
    )
    assert "deg_catalog_labels" in code
    assert "fuse_annotation" in code
    assert 'adata.obs["ref2_label"] = adata.obs["deg_label"]' not in code


def test_audit_rejects_azimuth_only():
    code = (
        "import scanpy as sc\nnp.random.seed(0)\n"
        "azimuth.RunAzimuth(adata)\n"
        "adata.obs['cell_type'] = adata.obs['predicted.celltype.l2']\n"
    )
    r = audit_code(code, {"tissue": "pbmc", "species": "human"}, phase="downstream")
    assert r["passed"] is False
    ids = {rec["id"] for rec in r["issue_records"]}
    assert "down.fusion" in ids or "down.azimuth_only" in ids
