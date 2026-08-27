from pathlib import Path

from sandbox.executor import write_and_maybe_run
from sandbox.jail import resolve_network, sandbox_settings
from sandbox.policy import policy_violations
from scagent.config import load_config


def _cfg(isolation: str = "rlimit"):
    base = load_config()
    analysis = {**(base.get("analysis") or {}), "executor": "subprocess"}
    return {
        **base,
        "analysis": analysis,
        "sandbox": {**sandbox_settings(base), "isolation": isolation, "enabled": True},
    }


def test_policy_blocks_os_system_token():
    hits = policy_violations("import os\nos.system('true')\n")
    assert "os.system" in hits
    assert not policy_violations("import scanpy as sc\nprint(sc.__name__)\n")


def test_executor_refuses_os_system(tmp_path):
    r = write_and_maybe_run(
        "import os\nos.system('true')\n",
        workspace=tmp_path,
        execute=True,
        timeout=5,
        cfg=_cfg(),
    )
    assert r["ok"] is False
    assert r["executed"] is False
    assert r["jail"] == "policy"
    assert "os.system" in r["stderr"]


def test_executor_isolated_home(tmp_path):
    r = write_and_maybe_run(
        "from pathlib import Path\nprint(Path.home())\n(Path.home() / 'marker.txt').write_text('ok')\n",
        workspace=tmp_path,
        execute=True,
        timeout=10,
        cfg=_cfg("rlimit"),
    )
    assert r["ok"] is True, r["stderr"]
    home = Path(r["stdout"].strip().splitlines()[-1])
    assert tmp_path.resolve() in home.resolve().parents or home.resolve() == (tmp_path / "sandbox_home").resolve()
    assert (tmp_path / "sandbox_home" / "marker.txt").read_text() == "ok"


def test_executor_timeout_kills(tmp_path):
    r = write_and_maybe_run(
        "import time\ntime.sleep(30)\n",
        workspace=tmp_path,
        execute=True,
        timeout=1,
        cfg=_cfg("rlimit"),
    )
    assert r["ok"] is False
    assert r["executed"] is True
    assert "timeout" in (r["stderr"] or "").lower() or r["returncode"] != 0


def test_seatbelt_blocks_write_outside_workspace(tmp_path):
    victim = tmp_path / "victim.txt"
    r = write_and_maybe_run(
        f"open({str(victim.resolve())!r}, 'w').write('x')\nprint('wrote')\n",
        workspace=tmp_path / "ws",
        execute=True,
        timeout=10,
        cfg=_cfg("auto"),
    )
    if r.get("jail") != "seatbelt":
        return
    assert r["ok"] is False
    assert not victim.exists()


def test_resolve_network_auto_by_phase():
    s = {"network": "auto"}
    assert resolve_network(s, phase="qc") is False
    assert resolve_network(s, phase="interpret") is False
    assert resolve_network(s, phase="downstream") is True
    assert resolve_network({"network": False}, phase="downstream") is False
