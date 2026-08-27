from pathlib import Path

from langgraph.checkpoint.memory import MemorySaver

from scagent.cli import build_parser
from scagent.demo import write_tiny_h5ad
from scagent.hitl import (
    build_mt_decision,
    build_resolution_decision,
    pick_option,
    render_decision_html,
    sample_pct_mt,
    write_decision_card,
)
from workflows.scRNA_langgraph import run_analysis


def test_mt_options_are_two_or_three_and_not_fixed_five_percent():
    card = build_mt_decision({"tissue": "heart", "metadata": {"tissue": "heart"}})
    opts = card["options"]
    assert 2 <= len(opts) <= 3
    ids = {o["id"] for o in opts}
    assert "recommended" in ids
    assert all(o.get("nmads") is not None for o in opts)
    html = render_decision_html(card)
    assert "pct_counts_mt" in html
    rec = pick_option(card, "recommended")
    assert rec["nmads"] >= 6
    strict = pick_option(card, "strict")
    assert strict["nmads"] < rec["nmads"]
    assert all("hard" not in o for o in opts)


def test_resolution_options_and_numeric_choice():
    card = build_resolution_decision({"tissue": "pbmc", "metadata": {"n_cells": 8000}})
    opts = card["options"]
    assert 2 <= len(opts) <= 3
    rec = pick_option(card, "recommended")
    assert rec["resolution"] is not None
    custom = pick_option(card, "0.4")
    assert float(custom["resolution"]) == 0.4
    html = render_decision_html(card)
    assert "resolution" in html.lower()


def test_mt_histogram_from_tiny_h5ad(tmp_path):
    path = write_tiny_h5ad(tmp_path / "tiny.h5ad")
    sample = sample_pct_mt(str(path), species="human")
    assert sample["n_sampled"] == 100
    assert sample["n_mt_genes"] >= 1
    assert sample["values"]
    card = build_mt_decision(
        {"tissue": "pbmc", "data_path": str(path), "metadata": {"species": "human", "tissue": "pbmc"}},
        sample=sample,
    )
    assert card["histogram"]["counts"]
    html_path = write_decision_card(card)
    assert html_path.is_file()
    assert (html_path.parent / "mt.json").is_file()


def test_interrupt_both_choices_uses_fixed_resolution(tmp_path):
    state = run_analysis(
        "HITL 两处都确认",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        interrupt_after_qc=True,
        auto_confirm=False,
        qc_choice="strict",
        resolution_choice="coarse",
        checkpointer=MemorySaver(),
    )
    assert state.get("code_qc")
    assert "nmads=4" in (state.get("code_qc") or "")
    assert state.get("code_downstream")
    assert "chosen_resolution = 0.2" in (state.get("code_downstream") or "")
    assert state.get("status") != "awaiting_mt_confirmation"
    assert state.get("status") != "awaiting_resolution_confirmation"


def test_auto_run_still_writes_decision_cards(tmp_path):
    state = run_analysis(
        "自动路径也要有直方图",
        data_path=str(tmp_path / "missing.h5ad"),
        tissue="pbmc",
        execute_code=False,
        checkpointer=MemorySaver(),
    )
    assert state.get("hitl_mt")
    assert state.get("hitl_resolution")
    assert Path("outputs/decisions/mt.html").is_file()
    assert Path("outputs/decisions/resolution.html").is_file()


def test_cli_confirm_parser():
    p = build_parser()
    args = p.parse_args(["confirm", "mt", "recommended"])
    assert args.kind == "mt"
    assert args.choice == "recommended"
    args = p.parse_args(["confirm", "resolution", "0.4"])
    assert args.kind == "resolution"
    assert args.choice == "0.4"
