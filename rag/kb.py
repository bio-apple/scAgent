"""Refresh the local RAG corpus: sc-best-practices book + lab SOPs."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
import urllib.request
import zipfile
from datetime import datetime, timezone
from pathlib import Path

from rag.ingest import collection_dir, ingest
from rag.retriever import clear_retrieve_cache
from scagent.config import REPO_ROOT, load_config, resolve_path

SCBP_GIT = "https://github.com/theislab/single-cell-best-practices.git"
DOC_EXTS = {".md", ".txt", ".pdf", ".ipynb"}
SKIP_DIR_NAMES = {"_build", "figures", "datasets", ".git", "__pycache__", ".ipynb_checkpoints", "node_modules"}


def cache_repo_dir(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    return resolve_path(cfg, "cache") / "sc-best-practices"


def book_dest(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    root = Path(cfg.get("_root") or REPO_ROOT)
    return root / "best_practices" / "upstream"


def sops_dir(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    return collection_dir(cfg, "sops")


def _reindex(cfg: dict | None) -> Path:
    path = ingest(cfg, force=True)
    clear_retrieve_cache()
    return path


def _book_root(repo: Path) -> Path:
    for cand in (repo / "jupyter-book", repo):
        if (cand / "_toc.yml").exists() or (cand / "_config.yml").exists():
            return cand
    return repo


def _zip_url(git_url: str, branch: str) -> str:
    base = git_url.rstrip("/").removesuffix(".git")
    return f"{base}/archive/refs/heads/{branch}.zip"


def _git(*args: str, cwd: Path | None = None) -> subprocess.CompletedProcess:
    git = shutil.which("git")
    if not git:
        raise FileNotFoundError("git")
    return subprocess.run([git, *args], cwd=cwd, capture_output=True, text=True)


def fetch_repo(dest: Path, *, url: str = SCBP_GIT, branch: str = "main") -> Path:
    """Clone or fast-forward the sc-best-practices repo. Zipball fallback if git is missing."""
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    git = shutil.which("git")
    if git:
        if (dest / ".git").is_dir():
            pulled = _git("-C", str(dest), "pull", "--ff-only")
            if pulled.returncode == 0:
                return dest
        else:
            if dest.exists():
                shutil.rmtree(dest)
            cloned = _git("clone", "--depth", "1", "--branch", branch, url, str(dest))
            if cloned.returncode != 0 and branch == "main":
                cloned = _git("clone", "--depth", "1", "--branch", "master", url, str(dest))
            if cloned.returncode == 0:
                return dest
    return download_zipball(dest, zip_url=_zip_url(url, branch))


def download_zipball(dest: Path, *, zip_url: str) -> Path:
    dest = Path(dest)
    dest.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory() as td:
        zpath = Path(td) / "book.zip"
        try:
            urllib.request.urlretrieve(zip_url, zpath)
        except Exception as exc:
            raise RuntimeError(f"下载失败 {zip_url}: {exc}") from exc
        with zipfile.ZipFile(zpath) as zf:
            zf.extractall(td)
        extracted = [p for p in Path(td).iterdir() if p.is_dir()]
        if not extracted:
            raise RuntimeError(f"zip 为空: {zip_url}")
        root = extracted[0]
        if dest.exists():
            shutil.rmtree(dest)
        shutil.copytree(root, dest)
    return dest


def sync_book(src: Path, dest: Path) -> int:
    """Copy md/txt/ipynb chapters out of a local clone into dest (replaces dest)."""
    src = Path(src)
    dest = Path(dest)
    root = _book_root(src)
    if dest.exists():
        shutil.rmtree(dest)
    dest.mkdir(parents=True, exist_ok=True)
    n = 0
    for path in root.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in {".md", ".txt", ".ipynb"}:
            continue
        if any(part in SKIP_DIR_NAMES for part in path.parts):
            continue
        if path.name in {"README.md"}:
            continue
        rel = path.relative_to(root)
        target = dest / rel
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        n += 1
    return n


def _commit_id(repo: Path) -> str | None:
    if not (repo / ".git").is_dir():
        return None
    try:
        r = _git("-C", str(repo), "rev-parse", "--short", "HEAD")
    except FileNotFoundError:
        return None
    if r.returncode != 0:
        return None
    return (r.stdout or "").strip() or None


def update_kb(
    *,
    url: str | None = None,
    branch: str = "main",
    cfg: dict | None = None,
    repo_dir: Path | None = None,
    dest: Path | None = None,
    fetch: bool = True,
    reindex: bool = True,
) -> dict:
    """Pull theislab/single-cell-best-practices and index it under best_practices/upstream."""
    cfg = cfg or load_config()
    url = url or SCBP_GIT
    repo = Path(repo_dir) if repo_dir else cache_repo_dir(cfg)
    out = Path(dest) if dest else book_dest(cfg)
    if fetch:
        fetch_repo(repo, url=url, branch=branch)
    n = sync_book(repo, out)
    commit = _commit_id(repo)
    (out / ".source.json").write_text(
        json.dumps(
            {
                "url": url,
                "branch": branch,
                "commit": commit,
                "n_files": n,
                "pulled_at": datetime.now(timezone.utc).isoformat(),
            },
            ensure_ascii=False,
            indent=2,
        )
        + "\n",
        encoding="utf-8",
    )
    index = str(_reindex(cfg)) if reindex else None
    return {"ok": True, "n_files": n, "repo": str(repo), "dest": str(out), "commit": commit, "index": index}


def add_doc(
    src: str | Path,
    *,
    dest_dir: Path | None = None,
    name: str | None = None,
    cfg: dict | None = None,
    reindex: bool = True,
) -> dict:
    """Copy a lab SOP (file or directory) into knowledge/sops and rebuild the index."""
    cfg = cfg or load_config()
    src = Path(src).expanduser().resolve()
    if not src.exists():
        raise FileNotFoundError(src)
    dest_root = Path(dest_dir) if dest_dir else sops_dir(cfg)
    dest_root.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []

    def _copy_file(path: Path, target: Path) -> None:
        if path.suffix.lower() not in DOC_EXTS:
            return
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        copied.append(str(target))

    if src.is_file():
        if src.suffix.lower() not in DOC_EXTS:
            raise ValueError(f"不支持的文件类型 {src.suffix}（需要 md/txt/pdf/ipynb）")
        target = dest_root / (name or src.name)
        _copy_file(src, target)
    else:
        if name:
            dest_root = dest_root / name
            dest_root.mkdir(parents=True, exist_ok=True)
        for path in src.rglob("*"):
            if not path.is_file() or any(part in SKIP_DIR_NAMES for part in path.parts):
                continue
            _copy_file(path, dest_root / path.relative_to(src))
        if not copied:
            raise ValueError(f"{src} 里没有 md/txt/pdf/ipynb")
    index = str(_reindex(cfg)) if reindex else None
    return {"ok": True, "n_files": len(copied), "files": copied, "dest": str(dest_root), "index": index}
