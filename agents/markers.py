from __future__ import annotations

import csv
import json
from pathlib import Path

from scagent.config import REPO_ROOT

DEFAULT_CATALOG = REPO_ROOT / "knowledge" / "markers" / "catalog.json"
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
    }


def load_marker_catalog(path: str | Path | None = None, tissue: str | None = None) -> dict:
    catalog_path = Path(path) if path else DEFAULT_CATALOG
    if catalog_path.suffix.lower() in {".csv", ".tsv"}:
        types = [_normalize_entry(r) for r in _from_table(catalog_path)]
        return {"tissue": tissue or "custom", "cell_types": types}
    if catalog_path.exists():
        data = json.loads(catalog_path.read_text(encoding="utf-8"))
        tissues = data.get("tissues") or {}
        key = None
        if tissue:
            key = next((k for k in tissues if k.lower() == tissue.lower()), None)
        if key is None:
            key = "pbmc" if "pbmc" in tissues else (next(iter(tissues), None))
        types = [_normalize_entry(r) for r in (tissues.get(key) or [])]
        return {"tissue": key or tissue or "unknown", "cell_types": types}
    return {"tissue": tissue or "unknown", "cell_types": []}


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
