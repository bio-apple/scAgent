"""Parse academic PDFs into section-structured Markdown for RAG."""

from __future__ import annotations

import json
import re
import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from scagent.config import REPO_ROOT, load_config, resolve_path

# Sections kept for RAG (references excluded).
KEEP_SECTIONS = frozenset(
    {"title", "abstract", "introduction", "methods", "results", "discussion", "conclusion", "other"}
)
DROP_SECTIONS = frozenset({"references", "acknowledgments", "supplementary"})

SECTION_HEADING_RE = re.compile(
    r"^(?:#+\s*)?"
    r"(?:\d+(?:\.\d+)*\.?\s+)?"
    r"(abstract|summary|introduction|background|"
    r"methods?|materials?\s+and\s+methods?|data\s+description|study\s+design|"
    r"statistical\s+analysis|experimental\s+(?:design|procedures?)|"
    r"results?|findings|"
    r"discussion|"
    r"conclusions?|"
    r"references|bibliography|literature\s+cited|"
    r"acknowledgments?|supplementary(?:\s+material)?)"
    r"\s*$",
    re.I | re.M,
)

CANONICAL_SECTION = {
    "abstract": "abstract",
    "summary": "abstract",
    "introduction": "introduction",
    "background": "introduction",
    "method": "methods",
    "methods": "methods",
    "materials and methods": "methods",
    "experimental design": "methods",
    "experimental procedures": "methods",
    "experimental procedure": "methods",
    "result": "results",
    "results": "results",
    "findings": "results",
    "discussion": "discussion",
    "data description": "methods",
    "study design": "methods",
    "statistical analysis": "methods",
    "conclusion": "conclusion",
    "conclusions": "conclusion",
    "references": "references",
    "bibliography": "references",
    "literature cited": "references",
    "acknowledgment": "acknowledgments",
    "acknowledgments": "acknowledgments",
    "acknowledgement": "acknowledgments",
    "acknowledgements": "acknowledgments",
    "supplementary": "supplementary",
    "supplementary material": "supplementary",
}


@dataclass
class ParsedPaper:
    source_pdf: str
    title: str
    sections: dict[str, str] = field(default_factory=dict)
    backend: str = "mineru"
    markdown: str = ""

    def to_markdown(self) -> str:
        if self.markdown:
            return self.markdown
        parts: list[str] = []
        if self.title:
            parts.append(f"# {self.title.strip()}\n")
        order = ("abstract", "introduction", "methods", "results", "discussion", "conclusion", "other")
        seen: set[str] = set()
        for key in order:
            body = (self.sections.get(key) or "").strip()
            if not body:
                continue
            seen.add(key)
            heading = key.replace("_", " ").title()
            if key == "methods":
                heading = "Methods"
            parts.append(f"## {heading}\n\n{body}\n")
        for key, body in self.sections.items():
            if key in seen or key in DROP_SECTIONS or key == "title":
                continue
            body = (body or "").strip()
            if body:
                parts.append(f"## {key.title()}\n\n{body}\n")
        return "\n".join(parts).strip() + "\n"


def _normalize_section(name: str) -> str:
    key = re.sub(r"\s+", " ", name.strip().lower())
    key = re.sub(r"^#+\s*", "", key)
    key = re.sub(r"^\d+(?:\.\d+)*\.?\s*", "", key)
    key = re.sub(r"^(section|appendix)\s+\d+(?:\.\d+)*\.?\s*", "", key)
    return CANONICAL_SECTION.get(key, "other")


def _clean_line(line: str) -> str:
    line = line.replace("\x00", " ")
    line = re.sub(r"\s+", " ", line).strip()
    return line


def _dehyphenate(text: str) -> str:
    return re.sub(r"(\w)-\n(\w)", r"\1\2", text)


# Post-process noisy PDF header/footer lines before section splitting.
_JUNK_LINE_RES = (
    re.compile(r"^(open access|research|review|perspective|brief report|r e s e a r c h)\s*$", re.I),
    re.compile(r"^(nature methods|nature reviews genetics|genome biology|nature communications)\b", re.I),
    re.compile(r"^received:\s*", re.I),
    re.compile(r"^(accepted|published online|article|correspondence):\s*", re.I),
    re.compile(r"^citation:\s*", re.I),
    re.compile(r"^doi:\s*", re.I),
    re.compile(r"^https?://", re.I),
    re.compile(r"^©|^copyright", re.I),
    re.compile(r"^\d+\s*$"),
    re.compile(r"^page \d+", re.I),
)


def _postprocess_text(text: str) -> str:
    text = sanitize_paper_text(text or "")
    text = _dehyphenate(text)
    text = re.sub(r"\r\n?", "\n", text)
    lines: list[str] = []
    for raw in text.split("\n"):
        line = _clean_line(raw)
        if not line:
            if lines and lines[-1] != "":
                lines.append("")
            continue
        if any(p.search(line) for p in _JUNK_LINE_RES):
            continue
        if len(line) < 4 and not line[0].isalpha():
            continue
        lines.append(line)
    return re.sub(r"\n{3,}", "\n\n", "\n".join(lines)).strip()


def sanitize_paper_text(text: str) -> str:
    """Remove PDF extraction garble: control chars, HTML sub/sup, entities, soft hyphens."""
    import html as _html

    text = text or ""
    # Soft hyphen / zero-width / BOM / replacement char
    text = text.replace("\ufeff", "").replace("\u00ad", "").replace("�", "")
    text = text.replace("\u200b", "").replace("\u200c", "").replace("\u200d", "")
    # Drop C0/C1 controls except newline/tab
    text = re.sub(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f-\x9f]", "", text)
    # MinerU / HTML leftovers: keep inner text of sub/sup/span
    text = re.sub(r"</?(?:sub|sup|span|font|b|i|em|strong)\b[^>]*>", "", text, flags=re.I)
    text = re.sub(r"<br\s*/?>", "\n", text, flags=re.I)
    text = _html.unescape(text)
    # Common PDF glyph loss: "10 Genomics" was "10x Genomics"
    text = re.sub(r"\b10\s+Genomics\b", "10x Genomics", text)
    text = re.sub(r"\b10\s*×\s*Genomics\b", "10x Genomics", text, flags=re.I)
    # Collapse runs of weird separators left by glyph stripping
    text = re.sub(r"[ \t]{2,}", " ", text)
    text = re.sub(r" *\n *", "\n", text)
    return text


def _clean_mineru_markdown(text: str) -> str:
    text = re.sub(r"<details>.*?</details>", "", text or "", flags=re.I | re.S)
    text = re.sub(r"!\[[^\]]*\]\([^)]+\)", "", text)
    text = sanitize_paper_text(text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def _extract_pdf_pymupdf(path: Path) -> str:
    import fitz  # pymupdf

    doc = fitz.open(str(path))
    pages: list[str] = []
    for page in doc:
        blocks = page.get_text("blocks") or []
        # blocks: x0, y0, x1, y1, text, block_no, block_type (0=text)
        text_blocks = [b for b in blocks if len(b) >= 7 and b[6] == 0 and str(b[4] or "").strip()]
        text_blocks.sort(key=lambda b: (round(float(b[1]) / 12), float(b[0])))
        page_text = "\n\n".join(str(b[4]).strip() for b in text_blocks)
        if not page_text.strip():
            page_text = page.get_text("text", sort=True) or ""
        pages.append(page_text)
    doc.close()
    return _postprocess_text("\n\n".join(pages))


def _extract_pdf_pdfplumber(path: Path) -> str:
    import pdfplumber

    parts: list[str] = []
    with pdfplumber.open(str(path)) as pdf:
        for page in pdf.pages:
            text = page.extract_text(
                layout=True,
                x_tolerance=2,
                y_tolerance=3,
                keep_blank_chars=False,
            )
            if not (text or "").strip():
                text = page.extract_text(x_tolerance=2, y_tolerance=2) or ""
            parts.append(text or "")
    return _postprocess_text("\n".join(parts))


def _extract_pdf_pypdf(path: Path) -> str:
    from pypdf import PdfReader

    reader = PdfReader(str(path))
    pages = [page.extract_text() or "" for page in reader.pages]
    return _postprocess_text("\n".join(pages))


def _extract_pdf_marker(path: Path) -> str:
    exe = shutil.which("marker_single") or shutil.which("marker")
    if not exe:
        raise RuntimeError("marker CLI not found on PATH")
    out_dir = path.parent / ".marker_tmp"
    out_dir.mkdir(parents=True, exist_ok=True)
    subprocess.run([exe, str(path), "--output_dir", str(out_dir)], check=True, capture_output=True, timeout=600)
    md_files = sorted(out_dir.rglob("*.md"), key=lambda p: p.stat().st_mtime, reverse=True)
    if not md_files:
        raise RuntimeError("marker produced no markdown")
    text = md_files[0].read_text(encoding="utf-8", errors="replace")
    shutil.rmtree(out_dir, ignore_errors=True)
    return text


def _mineru_settings(cfg: dict | None = None) -> dict:
    papers_cfg = ((cfg or load_config()).get("rag") or {}).get("papers") or {}
    base = {
        "cli": "auto",
        "backend": "pipeline",
        "method": "auto",
        "lang": None,
        "api_url": None,
        "formula": True,
        "table": True,
        "timeout_seconds": 1200,
        "keep_tmp": False,
    }
    raw = papers_cfg.get("mineru") or {}
    if isinstance(raw, dict):
        base.update({k: v for k, v in raw.items() if v is not None})
    return base


def _resolve_mineru_cli(cfg: dict | None = None) -> tuple[str, bool]:
    """Return (executable path, is_legacy_magic_pdf)."""
    pref = str(_mineru_settings(cfg).get("cli") or "auto").lower()
    if pref in {"magic-pdf", "magic_pdf"}:
        exe = shutil.which("magic-pdf")
        if not exe:
            raise RuntimeError("magic-pdf not found on PATH (pip install magic-pdf[full])")
        return exe, True
    if pref == "mineru":
        exe = shutil.which("mineru")
        if not exe:
            raise RuntimeError("mineru not found on PATH (pip install mineru)")
        return exe, False
    exe = shutil.which("mineru")
    if exe:
        return exe, False
    exe = shutil.which("magic-pdf")
    if exe:
        return exe, True
    raise RuntimeError(
        "MinerU CLI not found. Install `mineru` (2.x) or legacy `magic-pdf`: "
        "pip install mineru  OR  pip install magic-pdf[full]"
    )


def _build_mineru_cmd(exe: str, path: Path, out_dir: Path, *, cfg: dict | None, legacy: bool) -> list[str]:
    opts = _mineru_settings(cfg)
    cmd = [exe, "-p", str(path), "-o", str(out_dir)]
    if legacy:
        cmd.extend(["-m", str(opts.get("method") or "auto")])
        if opts.get("lang"):
            cmd.extend(["--lang", str(opts["lang"])])
        return cmd
    if opts.get("api_url"):
        cmd.extend(["--api-url", str(opts["api_url"])])
    if opts.get("backend"):
        cmd.extend(["-b", str(opts["backend"])])
    if opts.get("method"):
        cmd.extend(["-m", str(opts["method"])])
    if opts.get("lang"):
        cmd.extend(["-l", str(opts["lang"])])
    if opts.get("formula") is False:
        cmd.extend(["-f", "false"])
    if opts.get("table") is False:
        cmd.extend(["-t", "false"])
    return cmd


def _find_mineru_markdown(out_dir: Path, stem: str) -> Path | None:
    candidates = [p for p in out_dir.rglob("*.md") if p.is_file()]
    if not candidates:
        return None
    exact = [p for p in candidates if p.stem == stem]
    pool = exact or [
        p for p in candidates if not any(p.stem.endswith(s) for s in ("_layout", "_span", "_origin"))
    ]
    pool = pool or candidates
    return max(pool, key=lambda p: p.stat().st_size)


def _clean_mineru_markdown(text: str) -> str:
    text = re.sub(r"<details>.*?</details>", "", text or "", flags=re.I | re.S)
    text = re.sub(r"!\[[^\]]*\]\([^)]+\)", "", text)
    text = re.sub(r"\n{3,}", "\n\n", text)
    return text.strip()


def _extract_pdf_mineru(path: Path, *, cfg: dict | None = None) -> tuple[str, str]:
    cfg = cfg or load_config()
    opts = _mineru_settings(cfg)
    exe, legacy = _resolve_mineru_cli(cfg)
    out_dir = path.parent / ".mineru_tmp" / path.stem
    out_dir.mkdir(parents=True, exist_ok=True)
    cmd = _build_mineru_cmd(exe, path, out_dir, cfg=cfg, legacy=legacy)
    timeout = int(opts.get("timeout_seconds") or 1200)
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout, check=False)
    except subprocess.TimeoutExpired as exc:
        raise RuntimeError(f"MinerU timed out after {timeout}s") from exc
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout or "").strip()
        raise RuntimeError(f"MinerU failed ({proc.returncode}): {err[:500]}")
    md_path = _find_mineru_markdown(out_dir, path.stem)
    if not md_path:
        raise RuntimeError(f"MinerU produced no markdown under {out_dir}")
    text = _clean_mineru_markdown(md_path.read_text(encoding="utf-8", errors="replace"))
    if not opts.get("keep_tmp"):
        shutil.rmtree(out_dir, ignore_errors=True)
    backend_name = "magic-pdf" if legacy else "mineru"
    return text, backend_name


def _auto_backend_order(cfg: dict | None = None) -> list[str]:
    papers_cfg = ((cfg or load_config()).get("rag") or {}).get("papers") or {}
    raw = papers_cfg.get("auto_backends") or ["mineru", "pymupdf", "pdfplumber", "marker"]
    return [str(x).lower() for x in raw]


def extract_pdf_text(path: Path, backend: str = "auto", *, cfg: dict | None = None) -> tuple[str, str]:
    """Return (text, backend_used). Default auto chain excludes pypdf (low quality on two-column PDFs)."""
    cfg = cfg or load_config()
    b = (backend or "auto").lower()
    if b == "auto":
        order = _auto_backend_order(cfg)
    elif b == "high":
        order = ["marker", "mineru", "pymupdf", "pdfplumber"]
    else:
        order = [b]

    errors: list[str] = []
    for name in order:
        try:
            if name == "pymupdf":
                return _extract_pdf_pymupdf(path), "pymupdf"
            if name == "pdfplumber":
                return _extract_pdf_pdfplumber(path), "pdfplumber"
            if name == "pypdf":
                text = _postprocess_text(_extract_pdf_pypdf(path))
                return text, "pypdf"
            if name == "marker":
                return _postprocess_text(_extract_pdf_marker(path)), "marker"
            if name in {"mineru", "magic-pdf", "magic_pdf"}:
                text, used = _extract_pdf_mineru(path, cfg=cfg)
                return text, used
        except Exception as exc:
            errors.append(f"{name}: {exc}")
    hint = (
        "Install MinerU: pip install mineru (or legacy magic-pdf[full]), run mineru-models-download, "
        "then scagent parse-papers --backend mineru --force."
    )
    raise RuntimeError(f"PDF extraction failed ({path.name}): " + "; ".join(errors) + f". {hint}")


def _guess_title(lines: list[str]) -> str:
    for line in lines[:40]:
        line = _clean_line(line)
        if not line or len(line) < 8:
            continue
        if re.match(r"^(doi|http|www\.|arxiv|manuscript|preprint)", line, re.I):
            continue
        if line.lower() in CANONICAL_SECTION or SECTION_HEADING_RE.match(line):
            continue
        if sum(c.isalpha() for c in line) < len(line) * 0.5:
            continue
        return line[:240]
    return "Untitled"


def split_into_sections(raw: str) -> dict[str, str]:
    """Split plain or markdown text into canonical sections; drop references."""
    text = _dehyphenate(raw or "")
    text = re.sub(r"\r\n?", "\n", text)
    lines = text.split("\n")

    # If markdown ## headers exist, prefer them.
    if re.search(r"(?m)^#{1,3}\s+\w", text):
        return _split_markdown_sections(text)

    sections: dict[str, list[str]] = {}
    current = "other"
    title_lines: list[str] = []
    started = False

    for line in lines:
        clean = _clean_line(line)
        if not clean:
            if started:
                sections.setdefault(current, []).append("")
            continue
        m = SECTION_HEADING_RE.match(clean)
        if m:
            name = _normalize_section(m.group(1))
            if name in DROP_SECTIONS:
                current = name
                break
            if name in KEEP_SECTIONS:
                current = name
                started = True
                continue
        if not started and current == "other" and not sections:
            title_lines.append(clean)
            if len(title_lines) >= 3:
                started = True
                current = "other"
            continue
        started = True
        sections.setdefault(current, []).append(clean)

    out: dict[str, str] = {}
    title = _guess_title(title_lines) if title_lines else ""
    if title:
        out["title"] = title
    for key, buf in sections.items():
        if key in DROP_SECTIONS:
            continue
        body = "\n".join(buf).strip()
        if body:
            out[key] = body
    return out


def _split_markdown_sections(text: str) -> dict[str, str]:
    sections: dict[str, str] = {}
    title = ""
    current = "other"
    buf: list[str] = []

    def flush() -> None:
        nonlocal buf
        body = "\n".join(buf).strip()
        if body and current not in DROP_SECTIONS:
            sections[current] = body
        buf = []

    for line in text.splitlines():
        if line.startswith("# ") and not title:
            title = line.lstrip("# ").strip()
            continue
        m = re.match(r"^#{1,3}\s+(.+?)\s*$", line)
        if m:
            flush()
            current = _normalize_section(m.group(1))
            if current in DROP_SECTIONS:
                break
            continue
        buf.append(line)
    flush()
    if title:
        sections["title"] = title
    return sections


def parse_paper(path: str | Path, *, backend: str = "auto", cfg: dict | None = None) -> ParsedPaper:
    pdf = Path(path)
    cfg = cfg or load_config()
    papers_cfg = (cfg.get("rag") or {}).get("papers") or {}
    use_backend = backend or str(papers_cfg.get("parser_backend") or "mineru")
    raw, used = extract_pdf_text(pdf, backend=use_backend, cfg=cfg)
    raw = sanitize_paper_text(raw)
    sections = split_into_sections(raw)
    title = sanitize_paper_text(sections.pop("title", "") or _guess_title(raw.splitlines()))
    cleaned_sections = {k: sanitize_paper_text(v) for k, v in sections.items()}
    paper = ParsedPaper(source_pdf=str(pdf), title=title, sections=cleaned_sections, backend=used)
    paper.markdown = sanitize_paper_text(paper.to_markdown())
    return paper


def parsed_paths(cfg: dict | None = None) -> tuple[Path, Path]:
    cfg = cfg or load_config()
    papers_cfg = (cfg.get("rag") or {}).get("papers") or {}
    papers_dir = resolve_path(cfg, "knowledge") / "papers"
    rel = papers_cfg.get("parsed_dir") or "knowledge/papers/.parsed"
    parsed_root = Path(rel)
    if not parsed_root.is_absolute():
        parsed_root = Path(cfg.get("_root") or REPO_ROOT) / rel
    return papers_dir, parsed_root


def write_parsed_artifacts(paper: ParsedPaper, out_dir: Path) -> tuple[Path, Path]:
    out_dir.mkdir(parents=True, exist_ok=True)
    stem = Path(paper.source_pdf).stem
    md_path = out_dir / f"{stem}.md"
    meta_path = out_dir / f"{stem}.meta.json"
    md_path.write_text(sanitize_paper_text(paper.to_markdown()), encoding="utf-8")
    meta = {
        "source_pdf": paper.source_pdf,
        "title": sanitize_paper_text(paper.title),
        "backend": paper.backend,
        "sections": list(paper.sections.keys()),
    }
    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
    return md_path, meta_path


def sanitize_parsed_dir(cfg: dict | None = None) -> list[dict]:
    """Re-clean existing .parsed/*.md without re-running PDF backends."""
    _, parsed_dir = parsed_paths(cfg)
    rows: list[dict] = []
    if not parsed_dir.is_dir():
        return rows
    for md_path in sorted(parsed_dir.glob("*.md")):
        before = md_path.read_text(encoding="utf-8", errors="replace")
        after = sanitize_paper_text(before)
        after = re.sub(r"\n{3,}", "\n\n", after).strip() + "\n"
        changed = after != before
        if changed:
            md_path.write_text(after, encoding="utf-8")
        meta_path = md_path.with_suffix("").with_name(md_path.stem + ".meta.json")
        # stem.meta.json sits beside stem.md
        meta_path = parsed_dir / f"{md_path.stem}.meta.json"
        if meta_path.exists():
            try:
                meta = json.loads(meta_path.read_text(encoding="utf-8"))
                if isinstance(meta.get("title"), str):
                    meta["title"] = sanitize_paper_text(meta["title"])
                    meta_path.write_text(json.dumps(meta, ensure_ascii=False, indent=2), encoding="utf-8")
            except Exception:
                pass
        rows.append({"md": str(md_path), "changed": changed, "chars": len(after)})
    return rows


def parse_papers_dir(
    papers_dir: str | Path | None = None,
    *,
    backend: str = "auto",
    force: bool = False,
    cfg: dict | None = None,
) -> list[dict]:
    """Parse all PDFs under knowledge/papers/ into .parsed/*.md (+ meta)."""
    cfg = cfg or load_config()
    root, parsed_dir = parsed_paths(cfg)
    src = Path(papers_dir) if papers_dir else root
    papers_cfg = (cfg.get("rag") or {}).get("papers") or {}
    backend = backend or str(papers_cfg.get("parser_backend") or "mineru")
    results: list[dict] = []

    for pdf in sorted(src.glob("*.pdf")):
        md_out = parsed_dir / f"{pdf.stem}.md"
        meta_out = parsed_dir / f"{pdf.stem}.meta.json"
        if (
            not force
            and md_out.exists()
            and md_out.stat().st_mtime >= pdf.stat().st_mtime
        ):
            results.append({"pdf": str(pdf), "md": str(md_out), "skipped": True})
            continue
        paper = parse_paper(pdf, backend=backend, cfg=cfg)
        md_path, meta_path = write_parsed_artifacts(paper, parsed_dir)
        results.append(
            {
                "pdf": str(pdf),
                "md": str(md_path),
                "meta": str(meta_path),
                "title": paper.title,
                "backend": paper.backend,
                "sections": list(paper.sections.keys()),
                "skipped": False,
            }
        )
    return results


def chunk_paper_markdown(
    text: str,
    *,
    source: str,
    stem: str,
    chunk_size: int = 900,
    overlap: int = 150,
) -> list[dict]:
    """Section-first chunking with recursive fallback inside long sections."""
    from rag.ingest import chunk_document

    sections = _split_markdown_sections(text) if "##" in text else {"body": text}
    if not sections:
        sections = {"body": text}
    chunks: list[dict] = []
    idx = 0
    for section, body in sections.items():
        if section in DROP_SECTIONS or section == "title":
            continue
        body = (body or "").strip()
        if not body:
            continue
        for part in chunk_document(body, chunk_size, overlap, mode="semantic"):
            if not part.strip():
                continue
            chunks.append(
                {
                    "chunk_index": idx,
                    "section": section,
                    "text": part.strip(),
                }
            )
            idx += 1
    if not chunks and text.strip():
        for part in chunk_document(text, chunk_size, overlap, mode="semantic"):
            chunks.append({"chunk_index": idx, "section": "other", "text": part.strip()})
            idx += 1
    for c in chunks:
        c["source"] = source
        c["stem"] = stem
    return chunks


def ensure_papers_parsed(cfg: dict | None = None, *, force: bool = False) -> int:
    cfg = cfg or load_config()
    papers_cfg = (cfg.get("rag") or {}).get("papers") or {}
    if not papers_cfg.get("auto_parse", False):
        return 0
    rows = parse_papers_dir(force=force, cfg=cfg)
    return sum(1 for r in rows if not r.get("skipped"))
