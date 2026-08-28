# Immune / PBMC canonical markers (≥2 required)

Never assign a cell type from a single gene. Authoritative table: `catalog.json` in this directory (CellMarker 2.0 / PanglaoDB subset).

| Cell type | Positive markers (≥2) | Negative |
|-----------|----------------------|----------|
| T cell | CD3D, CD3E, CD3G | MS4A1, CD14 |
| CD4 T | CD4, IL7R | CD8A |
| CD8 T | CD8A, CD8B | CD4 |
| NK | NKG7, GNLY, KLRD1 | CD3D |
| B cell | MS4A1, CD79A, CD19 | CD3D |
| Classical mono | CD14, LYZ, S100A8 | low FCGR3A |
| Non-classical mono | FCGR3A, MS4A7 | low CD14 |
| cDC | FCER1A, CST3 | CD14 |
| pDC | IL3RA, CLEC4C, LILRA4 | |
| Platelet | PPBP, PF4 | |
| Epithelial | EPCAM, KRT18 | PTPRC |

Human mito prefix `MT-`, mouse `mt-`. Hemoglobin `HB*` flags blood contamination.
