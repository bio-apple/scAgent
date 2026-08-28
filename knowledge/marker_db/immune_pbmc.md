# 免疫 / PBMC canonical markers（至少双标记）

禁止只用一个基因决定细胞类型。权威表见同目录 `catalog.json`（CellMarker 2.0 / PanglaoDB 子集）。

| 细胞类型 | 阳性 marker（≥2） | 阴性 |
|----------|-------------------|------|
| T 细胞 | CD3D, CD3E, CD3G | MS4A1, CD14 |
| CD4 T | CD4, IL7R | CD8A |
| CD8 T | CD8A, CD8B | CD4 |
| NK | NKG7, GNLY, KLRD1 | CD3D |
| B 细胞 | MS4A1, CD79A, CD19 | CD3D |
| 经典单核 | CD14, LYZ, S100A8 | FCGR3A 低 |
| 非经典单核 | FCGR3A, MS4A7 | CD14 低 |
| cDC | FCER1A, CST3 | CD14 |
| pDC | IL3RA, CLEC4C, LILRA4 | |
| 血小板 | PPBP, PF4 | |
| 上皮 | EPCAM, KRT18 | PTPRC |

人线粒体前缀 `MT-`，小鼠 `mt-`。血红蛋白 `HB*` 用于血源污染。
