"""Analysis-version compatibility for run_manifest.json and --resume."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

from scagent import __version__ as _VERSION

_SEMVER = re.compile(r"^\s*v?(\d+)(?:\.(\d+))?(?:\.(\d+))?", re.I)


class ResumeIncompatibleError(RuntimeError):
    """Raised when --resume would continue under a different major scAgent version."""


def scagent_version() -> str:
    return str(_VERSION)


def parse_semver(raw: str | None) -> tuple[int, int, int] | None:
    if raw is None:
        return None
    m = _SEMVER.match(str(raw).strip())
    if not m:
        return None
    return int(m.group(1)), int(m.group(2) or 0), int(m.group(3) or 0)


def _level(saved: tuple[int, int, int], current: tuple[int, int, int]) -> str:
    if saved == current:
        return "match"
    if saved[0] != current[0]:
        return "major"
    if saved[1] != current[1]:
        return "minor"
    return "patch"


def load_manifest(path: Path | dict | None) -> dict[str, Any]:
    if path is None:
        return {}
    if isinstance(path, dict):
        return dict(path)
    p = Path(path)
    if not p.is_file():
        return {}
    try:
        data = json.loads(p.read_text(encoding="utf-8"))
    except Exception:
        return {}
    return data if isinstance(data, dict) else {}


def check_resume_compat(
    manifest: Path | dict | None,
    *,
    current: str | None = None,
    force: bool = False,
) -> dict[str, Any]:
    """Compare saved scagent_version with the running package.

    Major mismatch refuses resume unless force=True. Minor mismatch warns but allows.
    Missing/legacy manifests warn and allow.
    """
    cur = current or scagent_version()
    cur_t = parse_semver(cur)
    data = load_manifest(manifest)
    missing_file = not data and not isinstance(manifest, dict) and (
        manifest is None or not Path(manifest).is_file()
    )
    saved_raw = data.get("scagent_version") if data else None
    warnings: list[str] = []
    if missing_file:
        msg = "无 run_manifest.json，跳过 scAgent 版本兼容检查。"
        warnings.append(msg)
        return {
            "ok": True,
            "action": "warn",
            "level": "missing",
            "saved": None,
            "current": cur,
            "message": msg,
            "warnings": warnings,
        }
    saved_t = parse_semver(saved_raw if saved_raw is not None else None)
    if saved_t is None or cur_t is None:
        msg = (
            f"run_manifest 无有效 scagent_version（saved={saved_raw!r}）；"
            "无法核对分析脚本与当前 scAgent 的兼容性，仍继续续跑。"
        )
        warnings.append(msg)
        return {
            "ok": True,
            "action": "warn",
            "level": "missing",
            "saved": saved_raw,
            "current": cur,
            "message": msg,
            "warnings": warnings,
        }
    level = _level(saved_t, cur_t)
    pair = f"{saved_raw} → {cur}"
    if level == "match":
        return {
            "ok": True,
            "action": "allow",
            "level": level,
            "saved": saved_raw,
            "current": cur,
            "message": "",
            "warnings": [],
        }
    if level == "patch":
        return {
            "ok": True,
            "action": "allow",
            "level": level,
            "saved": saved_raw,
            "current": cur,
            "message": "",
            "warnings": [],
        }
    if level == "minor":
        msg = (
            f"scAgent 次版本不同（{pair}）。生成脚本与 schema 可能已变；"
            "续跑结果需人工核对。"
        )
        warnings.append(msg)
        return {
            "ok": True,
            "action": "warn",
            "level": level,
            "saved": saved_raw,
            "current": cur,
            "message": msg,
            "warnings": warnings,
        }
    msg = (
        f"scAgent 主版本不同（{pair}）。分析脚本与 checkpoint 可能不兼容，拒绝续跑。"
        "若确认仍要继续，请加 --force-resume。"
    )
    if force:
        warnings.append(msg + " （已 --force-resume，继续）")
        return {
            "ok": True,
            "action": "warn",
            "level": level,
            "saved": saved_raw,
            "current": cur,
            "message": warnings[-1],
            "warnings": warnings,
        }
    return {
        "ok": False,
        "action": "refuse",
        "level": level,
        "saved": saved_raw,
        "current": cur,
        "message": msg,
        "warnings": [msg],
    }


def assert_resume_compatible(
    manifest: Path | dict | None,
    *,
    current: str | None = None,
    force: bool = False,
) -> dict[str, Any]:
    result = check_resume_compat(manifest, current=current, force=force)
    if not result.get("ok"):
        raise ResumeIncompatibleError(result.get("message") or "resume incompatible")
    return result
