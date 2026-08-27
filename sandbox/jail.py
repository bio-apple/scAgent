"""OS jail + rlimits for generated analysis processes."""

from __future__ import annotations

import os
import re
import shutil
import signal
import subprocess
import sys
from pathlib import Path

from scagent.logutil import get_logger

log = get_logger("sandbox.jail")

_SECRET_RE = re.compile(r"(API_KEY|SECRET|TOKEN|PASSWORD|CREDENTIAL|PRIVATE_KEY)", re.I)


def sandbox_settings(cfg: dict | None) -> dict:
    s = dict((cfg or {}).get("sandbox") or {})
    s.setdefault("enabled", True)
    s.setdefault("isolation", "auto")
    s.setdefault("network", "auto")
    s.setdefault("memory_mb", 16384)
    s.setdefault("max_processes", 64)
    s.setdefault("max_open_files", 4096)
    s.setdefault("static_policy", True)
    s.setdefault("cleanup_tmp", True)
    s.setdefault("docker_image", None)
    return s


def resolve_network(settings: dict, *, phase: str | None = None) -> bool:
    """QC defaults to no network; downstream may download CellTypist models."""
    raw = settings.get("network", "auto")
    if raw in {True, False}:
        return bool(raw)
    val = str(raw or "auto").lower()
    if val in {"1", "true", "yes", "on"}:
        return True
    if val in {"0", "false", "no", "off"}:
        return False
    return (phase or "") in {"downstream", "annotation"}


def isolated_env(workspace: Path, *, seed: int) -> dict[str, str]:
    env = {k: v for k, v in os.environ.items() if not _SECRET_RE.search(k)}
    home = workspace / "sandbox_home"
    tmp = workspace / "sandbox_tmp"
    cache = workspace / "sandbox_cache"
    for p in (home, tmp, cache, cache / "matplotlib"):
        p.mkdir(parents=True, exist_ok=True)
    env["HOME"] = str(home)
    env["TMPDIR"] = str(tmp)
    env["TMP"] = str(tmp)
    env["TEMP"] = str(tmp)
    env["XDG_CACHE_HOME"] = str(cache)
    env["XDG_CONFIG_HOME"] = str(home / ".config")
    env["MPLCONFIGDIR"] = str(cache / "matplotlib")
    env["NUMBA_CACHE_DIR"] = str(cache / "numba")
    env["JOBLIB_TEMP_FOLDER"] = str(tmp)
    env["PYTHONNOUSERSITE"] = "1"
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    env["PYTHONHASHSEED"] = str(seed)
    from scagent.config import REPO_ROOT

    pp = [str(REPO_ROOT)]
    if env.get("PYTHONPATH"):
        pp.append(env["PYTHONPATH"])
    env["PYTHONPATH"] = os.pathsep.join(pp)
    return env


def apply_rlimits(settings: dict) -> None:
    if os.name != "posix":
        return
    import resource

    def _set(name: str, value: int) -> None:
        flag = getattr(resource, name, None)
        if flag is None or value is None:
            return
        try:
            resource.setrlimit(flag, (value, value))
        except (ValueError, OSError) as exc:
            log.debug("rlimit %s=%s skipped: %s", name, value, exc)

    mem = int(settings.get("memory_mb") or 0)
    if mem > 0:
        _set("RLIMIT_AS", mem * 1024 * 1024)
    nproc = int(settings.get("max_processes") or 0)
    if nproc > 0:
        _set("RLIMIT_NPROC", nproc)
    nofile = int(settings.get("max_open_files") or 0)
    if nofile > 0:
        _set("RLIMIT_NOFILE", nofile)
    _set("RLIMIT_CORE", 0)


def _seatbelt_profile(workspace: Path, settings: dict) -> Path:
    ws = str(workspace.resolve())
    lines = [
        "(version 1)",
        "(allow default)",
        "(deny file-write*)",
        f'(allow file-write* (subpath "{ws}"))',
        '(allow file-write* (regex #"^(/private)?/tmp/"))',
        '(allow file-write* (regex #"^/dev/"))',
        '(deny process-exec (literal "/bin/sh"))',
        '(deny process-exec (literal "/bin/bash"))',
        '(deny process-exec (literal "/bin/zsh"))',
        '(deny process-exec (literal "/usr/bin/sudo"))',
        '(deny process-exec (literal "/bin/rm"))',
        '(deny process-exec (literal "/usr/bin/ssh"))',
    ]
    if not settings.get("network"):
        lines.append("(deny network*)")
    path = workspace / ".sandbox.sb"
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def wrap_argv(argv: list[str], workspace: Path, settings: dict) -> tuple[list[str], str]:
    isolation = str(settings.get("isolation") or "auto").lower()
    if isolation in {"off", "none"}:
        return argv, "off"
    if isolation == "rlimit":
        return argv, "rlimit"
    docker = shutil.which("docker") or shutil.which("podman")
    image = settings.get("docker_image") or os.getenv("SCAGENT_DOCKER_IMAGE")
    if isolation in {"docker", "podman"}:
        if docker and image:
            mem = int(settings.get("memory_mb") or 16384)
            pids = int(settings.get("max_processes") or 64)
            cmd = [docker, "run", "--rm", "-w", str(workspace.resolve()), f"--memory={mem}m", f"--pids-limit={pids}"]
            if not settings.get("network"):
                cmd += ["--network", "none"]
            cmd += ["-v", f"{workspace.resolve()}:{workspace.resolve()}", image, *argv]
            return cmd, "docker"
        log.warning("docker isolation requested but docker/image missing; rlimit fallback")
        return argv, "rlimit"
    if sys.platform == "darwin" and shutil.which("sandbox-exec"):
        profile = _seatbelt_profile(workspace, settings)
        return ["sandbox-exec", "-f", str(profile), *argv], "seatbelt"
    bwrap = shutil.which("bwrap")
    if bwrap:
        cmd = [
            bwrap,
            "--die-with-parent",
            "--ro-bind",
            "/",
            "/",
            "--dev",
            "/dev",
            "--proc",
            "/proc",
            "--bind",
            str(workspace.resolve()),
            str(workspace.resolve()),
        ]
        if not settings.get("network"):
            cmd.append("--unshare-net")
        cmd += argv
        return cmd, "bwrap"
    return argv, "rlimit"


def kill_process_group(proc: subprocess.Popen) -> None:
    if proc.poll() is not None:
        return
    if os.name == "posix":
        try:
            os.killpg(proc.pid, signal.SIGKILL)
            return
        except (ProcessLookupError, PermissionError, OSError):
            pass
    proc.kill()


def run_jailed(
    argv: list[str],
    *,
    workspace: Path,
    env: dict[str, str],
    timeout: int,
    settings: dict,
) -> subprocess.CompletedProcess:
    isolation = str(settings.get("isolation") or "auto").lower()
    wrapped, jail = wrap_argv(argv, workspace, settings)
    kwargs: dict = {
        "cwd": str(workspace),
        "env": env,
        "stdout": subprocess.PIPE,
        "stderr": subprocess.PIPE,
        "text": True,
    }
    if os.name == "posix":
        kwargs["start_new_session"] = True
        if isolation != "off":
            kwargs["preexec_fn"] = lambda: apply_rlimits(settings)

    def _launch(cmd: list[str]) -> subprocess.Popen:
        return subprocess.Popen(cmd, **kwargs)

    proc = _launch(wrapped)
    used = jail
    try:
        stdout, stderr = proc.communicate(timeout=timeout)
    except subprocess.TimeoutExpired:
        kill_process_group(proc)
        stdout, stderr = proc.communicate(timeout=5)
        out = subprocess.CompletedProcess(wrapped, -9, stdout or "", (stderr or "") + f"\ntimeout after {timeout}s")
        out.jail = used  # type: ignore[attr-defined]
        return out

    nested_deny = "sandbox_apply: Operation not permitted" in (stderr or "") or proc.returncode == 71
    if used in {"seatbelt", "bwrap"} and nested_deny:
        log.warning("OS jail unavailable (%s); falling back to rlimits", used)
        kwargs.pop("preexec_fn", None)
        if os.name == "posix" and isolation != "off":
            kwargs["preexec_fn"] = lambda: apply_rlimits(settings)
        proc = _launch(argv)
        used = "rlimit"
        try:
            stdout, stderr = proc.communicate(timeout=timeout)
        except subprocess.TimeoutExpired:
            kill_process_group(proc)
            stdout, stderr = proc.communicate(timeout=5)
            out = subprocess.CompletedProcess(argv, -9, stdout or "", (stderr or "") + f"\ntimeout after {timeout}s")
            out.jail = used  # type: ignore[attr-defined]
            return out
        wrapped = argv

    out = subprocess.CompletedProcess(wrapped, proc.returncode or 0, stdout or "", stderr or "")
    out.jail = used  # type: ignore[attr-defined]
    return out
