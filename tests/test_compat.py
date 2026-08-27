"""scAgent version is recorded in run_manifest and checked on --resume."""

from __future__ import annotations

import json
from pathlib import Path

import pytest
from langgraph.checkpoint.memory import MemorySaver

from sandbox.executor import write_and_maybe_run, write_manifest
from scagent import __version__
from scagent.compat import (
    ResumeIncompatibleError,
    assert_resume_compatible,
    check_resume_compat,
    parse_semver,
    scagent_version,
)
from workflows.scRNA_langgraph import run_analysis


def test_parse_semver():
    assert parse_semver("0.1.0") == (0, 1, 0)
    assert parse_semver("v1.2") == (1, 2, 0)
    assert parse_semver("2.0.0+local") == (2, 0, 0)
    assert parse_semver(None) is None
    assert parse_semver("not-a-version") is None


def test_match_and_patch_allow():
    same = check_resume_compat({"scagent_version": "0.1.0"}, current="0.1.0")
    assert same["ok"] is True
    assert same["action"] == "allow"
    assert same["level"] == "match"
    patch = check_resume_compat({"scagent_version": "0.1.0"}, current="0.1.3")
    assert patch["ok"] is True
    assert patch["level"] == "patch"
    assert patch["action"] == "allow"


def test_minor_warns_major_refuses():
    minor = check_resume_compat({"scagent_version": "0.1.0"}, current="0.2.0")
    assert minor["ok"] is True
    assert minor["action"] == "warn"
    assert minor["level"] == "minor"
    assert "次版本" in minor["message"]
    major = check_resume_compat({"scagent_version": "0.1.0"}, current="1.0.0")
    assert major["ok"] is False
    assert major["action"] == "refuse"
    assert major["level"] == "major"
    assert "主版本" in major["message"]
    with pytest.raises(ResumeIncompatibleError, match="主版本"):
        assert_resume_compatible({"scagent_version": "0.1.0"}, current="1.0.0")
    forced = check_resume_compat({"scagent_version": "0.1.0"}, current="1.0.0", force=True)
    assert forced["ok"] is True
    assert forced["level"] == "major"
    assert_resume_compatible({"scagent_version": "0.1.0"}, current="1.0.0", force=True)


def test_legacy_manifest_warns_and_allows(tmp_path: Path):
    missing = check_resume_compat(tmp_path / "nope.json", current="0.1.0")
    assert missing["ok"] is True
    assert missing["level"] == "missing"
    old = tmp_path / "run_manifest.json"
    old.write_text(json.dumps({"python": "3.12", "skills_fingerprint": "abc"}), encoding="utf-8")
    legacy = check_resume_compat(old, current="0.1.0")
    assert legacy["ok"] is True
    assert legacy["level"] == "missing"
    assert "scagent_version" in legacy["message"]


def test_write_manifest_records_scagent_version(tmp_path: Path):
    write_manifest(tmp_path, {"phase": "qc"})
    payload = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert payload["scagent_version"] == scagent_version() == __version__
    r = write_and_maybe_run("print('ok')", workspace=tmp_path, execute=False)
    assert r["ok"] is True
    again = json.loads((tmp_path / "run_manifest.json").read_text(encoding="utf-8"))
    assert again["scagent_version"] == __version__


def test_run_analysis_resume_refuses_major_mismatch(tmp_path: Path, monkeypatch):
    (tmp_path / "run_manifest.json").write_text(
        json.dumps({"scagent_version": "99.0.0", "python": "3.12"}),
        encoding="utf-8",
    )
    monkeypatch.setattr("workflows.scRNA_langgraph.resolve_path", lambda cfg, key: tmp_path)
    with pytest.raises(ResumeIncompatibleError, match="主版本"):
        run_analysis("续跑", resume=True, checkpointer=MemorySaver(), execute_code=False)
