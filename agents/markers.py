from __future__ import annotations

import csv
import json
from pathlib import Path

from scagent.config import REPO_ROOT

DEFAULT_CATALOG = REPO_ROOT / "knowledge" / "marker_db" / "catalog.json"
IMMUNE_TISSUES = {"pbmc", "blood", "immune"}

# Official CellTypist filenames (celltypist.org / models_description). Organ atlas names for heart/liver.
CELLTYPIST_MODELS = {
    "pbmc": "Immune_All_Low.pkl",
    "blood": "Immune_All_Low.pkl",
    "immune": "Immune_All_Low.pkl",
    "lung": "Human_Lung_Atlas.pkl",
    "airway": "Cells_Lung_Airway.pkl",
    "brain": "Developing_Human_Brain.pkl",
    "gut": "Cells_Intestinal_Tract.pkl",
    "intestine": "Cells_Intestinal_Tract.pkl",
    "colon": "Cells_Intestinal_Tract.pkl",
    "tumor": "Human_Colorectal_Cancer.pkl",
    "cancer": "Human_Colorectal_Cancer.pkl",
    "crc": "Human_Colorectal_Cancer.pkl",
    "embryo": "Pan_Fetal_Human.pkl",
    "fetal": "Pan_Fetal_Human.pkl",
    "heart": "Adult_Human_Heart.pkl",
    "liver": "Adult_Human_Liver.pkl",
    "kidney": "Adult_Human_Kidney.pkl",
}


def choose_celltypist_model(tissue: str | None, species: str | None = None) -> str | None:
    """Map tissue → CellTypist model. None = do not use Immune_All on a mismatched organ."""
    t = str(tissue or "default").lower()
    if t in {"default", "unknown", ""}:
        return None
    if str(species or "").lower() == "mouse" and t in {"brain"}:
        return "Developing_Mouse_Brain.pkl"
    if str(species or "").lower() == "mouse" and t in {"gut", "intestine"}:
        return "Adult_Mouse_Gut.pkl"
    return CELLTYPIST_MODELS.get(t)


def _normalize_entry(row: dict) -> dict:
    name = row.get("name") or row.get("cell_type") or "unknown"
    lin = row.get("lineage") or []
    if isinstance(lin, str):
        lin = [x.strip() for x in lin.replace("|", ";").split(";") if x.strip()]
    if not lin:
        lin = [name]
    pos = row.get("positive") or []
    neg = row.get("negative") or []
    if isinstance(pos, str):
        pos = [g.strip() for g in pos.replace("|", ";").split(";") if g.strip()]
    if isinstance(neg, str):
        neg = [g.strip() for g in neg.replace("|", ";").split(";") if g.strip()]
    return {
        "name": name,
        "lineage": lin,
        "level": len(lin),
        "positive": pos,
        "negative": neg,
        **{k: row[k] for k in ("cl_id", "source") if row.get(k)},
    }


def load_marker_catalog(path: str | Path | None = None, tissue: str | None = None) -> dict:
    """Load tissue marker catalog. Never silently map non-immune organs onto PBMC markers."""
    catalog_path = Path(path) if path else DEFAULT_CATALOG
    if catalog_path.suffix.lower() in {".csv", ".tsv"}:
        types = [_normalize_entry(r) for r in _from_table(catalog_path)]
        return {"tissue": tissue or "custom", "cell_types": types, "warning": None}
    if catalog_path.exists():
        data = json.loads(catalog_path.read_text(encoding="utf-8"))
        tissues = data.get("tissues") or {}
        aliases = {str(k).lower(): k for k in tissues}
        # Common aliases → catalog keys
        alias_map = {
            "blood": "pbmc",
            "immune": "pbmc",
            "airway": "lung",
            "cancer": "tumor",
            "crc": "tumor",
            "tme": "tumor",
        }
        key = None
        warning = None
        t_raw = (tissue or "").strip()
        t_low = t_raw.lower()
        if t_low:
            key = aliases.get(t_low) or aliases.get(alias_map.get(t_low, ""))
        if key is None:
            if t_low in IMMUNE_TISSUES or t_low in {"", "default"}:
                # Immune / unspecified demo default only — never for lung/tumor/brain/…
                key = aliases.get("pbmc") or (next(iter(tissues), None) if not t_low else None)
                if t_low in IMMUNE_TISSUES and key and key != "pbmc":
                    warning = f"immune tissue {t_raw!r} using catalog {key!r}"
            else:
                # Refuse wrong-organ fallback (was: silent PBMC for lung/tumor/embryo)
                return {
                    "tissue": t_raw or "unknown",
                    "cell_types": [],
                    "warning": (
                        f"no marker catalog for tissue={t_raw!r}; "
                        "refusing PBMC fallback — use reference mapping + literature markers"
                    ),
                }
        types = [_normalize_entry(r) for r in (tissues.get(key) or [])] if key else []
        return {"tissue": key or t_raw or "unknown", "cell_types": types, "warning": warning}
    return {"tissue": tissue or "unknown", "cell_types": [], "warning": "marker catalog file missing"}


def _from_table(path: Path) -> list[dict]:
    delim = "\t" if path.suffix.lower() == ".tsv" else ","
    rows = []
    with path.open(encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter=delim)
        for row in reader:
            rows.append(
                {
                    "name": row.get("cell_type") or row.get("name") or "unknown",
                    "lineage": row.get("lineage") or "",
                    "positive": row.get("positive") or row.get("positives") or "",
                    "negative": row.get("negative") or row.get("negatives") or "",
                }
            )
    return rows


def catalog_as_python(catalog: dict) -> str:
    return json.dumps(catalog.get("cell_types") or [], ensure_ascii=False, indent=2)


def max_lineage_depth(catalog: dict) -> int:
    types = catalog.get("cell_types") or []
    return max((len(t.get("lineage") or []) for t in types), default=1)


def is_immune_tissue(tissue: str | None) -> bool:
    return str(tissue or "").lower() in IMMUNE_TISSUES
