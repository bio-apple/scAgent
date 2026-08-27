from unittest.mock import patch

import pytest

from agents.common import invoke_llm, token_usage
from agents.templates import cluster_annotate_script, qc_preprocess_script
from scagent.config import analysis_params, cfg_get, load_config
from scagent.inspect_data import inspect_data
from scagent.io import read_single_cell
from scagent.logutil import get_logger, timed


def test_analysis_params_from_yaml():
    cfg = load_config()
    p = analysis_params(cfg)
    assert p["n_pcs"] == cfg_get(cfg, "params.n_pcs")
    assert p["n_hvg"] == 2000
    assert p["n_neighbors"] == 15
    assert cfg_get(cfg, "missing.key", "x") == "x"
    assert cfg_get(cfg, "model.api_key_env") == "OPENAI_API_KEY"


def test_templates_inject_config_params():
    p = analysis_params()
    qc = qc_preprocess_script(
        {"data_path": "x.h5ad", "species": "human", "tissue": "pbmc"},
        {"method": "mad", "nmads": 5},
    )
    assert f"n_top_genes={p['n_hvg']}" in qc
    assert f"target_sum={float(p['target_sum'])}" in qc
    down = cluster_annotate_script(
        {"data_path": "x.h5ad", "tissue": "pbmc"},
        {},
        {"integrator": None},
    )
    assert f"n_comps={p['n_pcs']}" in down
    assert f"n_neighbors={p['n_neighbors']}" in down
    rds = qc_preprocess_script({"data_path": "obj.rds", "species": "human", "tissue": "pbmc"}, {"nmads": 5})
    assert "read_single_cell" in rds


def test_io_dispatch_by_suffix(tmp_path):
    missing = tmp_path / "nope.h5ad"
    with pytest.raises(FileNotFoundError):
        read_single_cell(missing)
    bad = tmp_path / "x.csv"
    bad.write_text("a,b\n1,2\n")
    with pytest.raises(ValueError, match="unsupported"):
        read_single_cell(bad)


def test_inspect_rds_notes_without_r(tmp_path):
    p = tmp_path / "obj.rds"
    p.write_bytes(b"not-a-seurat")
    meta = inspect_data(str(p))
    assert meta["platform"] == "seurat"
    assert meta["exists"] is True
    assert any("Seurat" in n for n in meta["notes"])


def test_llm_retry_backoff():
    class Boom:
        def __init__(self):
            self.n = 0

        def invoke(self, _messages):
            self.n += 1
            if self.n < 3:
                raise ConnectionError("429 rate limit")

            class AI:
                content = "ok"
                usage_metadata = {"input_tokens": 10, "output_tokens": 2, "total_tokens": 12}

            return AI()

    model = Boom()
    before = token_usage()["total"]
    with patch("agents.common.time.sleep"):
        ai = invoke_llm(
            model,
            [],
            {
                "model": {"max_retries": 4, "retry_backoff_seconds": 0.01, "rate_limit_rpm": 0},
                "performance": {"cache": False},
            },
        )
    assert ai.content == "ok"
    assert model.n == 3
    assert token_usage()["total"] >= before + 12


def test_logging_timed():
    log = get_logger("test")
    with timed("unit.test", log):
        pass


def test_cli_from_checkpoint_alias():
    from scagent.cli import build_parser

    p = build_parser()
    assert p.parse_args(["run", "q", "--from-checkpoint"]).resume is True
    assert p.parse_args(["run", "q", "--resume"]).resume is True
    assert p.parse_args(["run", "q"]).resume is False
    assert p.parse_args(["memory"]).func.__name__ == "cmd_memory"


def test_writer_html_and_run_log(tmp_path):
    from agents.writer import render_html, write_run_log

    state = {
        "user_query": "demo",
        "report_lang": "zh",
        "thread_id": "t1",
        "metadata": {"species": "human", "tissue": "pbmc"},
        "plan": {"integrator": None, "skills": [], "route": ["qc"]},
        "review_qc": {"passed": True, "issues": [], "issue_records": []},
        "artifacts": {"figure_captions": [{"path": "figures/violin.png", "kind": "violin", "caption": "qc"}]},
    }
    html = render_html(state)
    assert "<!DOCTYPE html>" in html
    assert "figures/violin.png" in html
    write_run_log(state, tmp_path)
    import json

    log = json.loads((tmp_path / "run_log.json").read_text(encoding="utf-8"))
    assert log["thread_id"] == "t1"
    assert "issue_records" in log["review_qc"]
