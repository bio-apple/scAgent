import json
from pathlib import Path

import pytest

from rag.ingest import ingest
from rag.kb import add_doc, update_kb
from scagent.cli import build_parser


def _cfg(tmp_path: Path, collections: list[str], dirs: dict[str, str]) -> dict:
    return {
        "_root": str(tmp_path),
        "paths": {"knowledge": str(tmp_path / "knowledge"), "index": str(tmp_path / "index")},
        "rag": {
            "collections": collections,
            "collection_dirs": dirs,
            "chunk_size": 400,
            "chunk_overlap": 40,
            "chunking": "fixed",
        },
    }


def test_cli_update_kb_and_add_doc_parsers():
    p = build_parser()
    u = p.parse_args(["update-kb"])
    assert u.func.__name__ == "cmd_update_kb"
    assert u.branch == "main"
    a = p.parse_args(["add-doc", "lab.md", "--name", "qc.md"])
    assert a.func.__name__ == "cmd_add_doc"
    assert a.path == "lab.md"
    assert a.name == "qc.md"


def test_add_doc_file_is_indexed(tmp_path):
    src = tmp_path / "lab_qc.md"
    src.write_text("Lab SOP: inspect the mitochondrial histogram; never use a fixed 5 percent cutoff.\n")
    sops = tmp_path / "knowledge" / "sops"
    cfg = _cfg(tmp_path, ["sops"], {"sops": str(sops)})
    info = add_doc(src, dest_dir=sops, cfg=cfg, reindex=True)
    assert info["n_files"] == 1
    assert (sops / "lab_qc.md").is_file()
    chunks = (tmp_path / "index" / "chunks.jsonl").read_text(encoding="utf-8")
    assert "mitochondrial histogram" in chunks
    assert '"collection": "sops"' in chunks


def test_add_doc_directory_and_ipynb(tmp_path):
    folder = tmp_path / "incoming"
    folder.mkdir()
    (folder / "notes.md").write_text("Use MAD n=4 for this tissue.\n")
    nb = {
        "cells": [
            {"cell_type": "markdown", "source": ["# Doublet SOP\n", "Run Scrublet then scDblFinder.\n"]},
            {"cell_type": "code", "source": ["sc.pp.scrublet(adata)\n"]},
        ],
        "metadata": {},
        "nbformat": 4,
        "nbformat_minor": 5,
    }
    (folder / "doublets.ipynb").write_text(json.dumps(nb), encoding="utf-8")
    sops = tmp_path / "sops"
    cfg = _cfg(tmp_path, ["sops"], {"sops": str(sops)})
    info = add_doc(folder, dest_dir=sops, name="our-lab", cfg=cfg, reindex=True)
    assert info["n_files"] == 2
    assert (sops / "our-lab" / "notes.md").is_file()
    blob = (tmp_path / "index" / "chunks.jsonl").read_text(encoding="utf-8")
    assert "Scrublet then scDblFinder" in blob
    assert "MAD n=4" in blob


def test_add_doc_rejects_missing(tmp_path):
    with pytest.raises(FileNotFoundError):
        add_doc(tmp_path / "missing.md", dest_dir=tmp_path / "sops", reindex=False)


def test_update_kb_syncs_local_book_without_network(tmp_path):
    repo = tmp_path / "repo"
    chapter = repo / "jupyter-book" / "preprocessing_visualization"
    chapter.mkdir(parents=True)
    (repo / "jupyter-book" / "_toc.yml").write_text("format: jb-book\n", encoding="utf-8")
    (chapter / "quality_control.md").write_text(
        "# Quality control\nUse MAD-based filtering per sample.\n", encoding="utf-8"
    )
    dest = tmp_path / "knowledge" / "upstream"
    cfg = _cfg(tmp_path, ["upstream"], {"upstream": str(dest)})
    info = update_kb(repo_dir=repo, dest=dest, cfg=cfg, fetch=False, reindex=True)
    assert info["n_files"] >= 1
    assert (dest / "preprocessing_visualization" / "quality_control.md").is_file()
    assert (dest / ".source.json").is_file()
    chunks = (tmp_path / "index" / "chunks.jsonl").read_text(encoding="utf-8")
    assert "MAD-based filtering" in chunks
    assert ingest(cfg, force=False)


def test_retrieve_empty_sops_is_safe():
    from rag.retriever import clear_retrieve_cache, retrieve

    clear_retrieve_cache()
    hits = retrieve("mitochondrial MAD", collection="sops", top_k=3)
    assert hits == []
    merged = retrieve("mitochondrial MAD", collections=["sops", "best_practices"], top_k=3)
    assert merged


def test_config_includes_sops_collection():
    from scagent.config import load_config

    cfg = load_config()
    assert "sops" in (cfg.get("rag") or {}).get("collections")
    assert "marker_db" in (cfg.get("rag") or {}).get("collections")
    assert "cell_ontology" in (cfg.get("rag") or {}).get("collections")
    assert "upstream" in (cfg.get("rag") or {}).get("collections")
    assert (cfg.get("rag") or {}).get("default_collection") == "fused"
