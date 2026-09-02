# scAgent Hard Rules

These rules are **non-negotiable**. They override ad-hoc prompts, skill excerpts, and model improvisation.
Violations must **fail review** or **block downstream phases** — never silently proceed.

## 1. QC must pass before clustering

- Do **not** run neighbors / UMAP / Leiden / annotation until QC review passes.
- QC must inspect **`nFeature_RNA` / `n_genes_by_counts`**, **`nCount_RNA` / `total_counts`**, and **`percent.mt` / `pct_counts_mt`** (or tissue-appropriate equivalents).
- QC must produce **Violin / VlnPlot** (plus scatter + MAD diagnostics per SOP).
- If QC fails or is missing, route to **review only** — no clustering script generation or execution.

## 2. Differential expression must use pseudobulk for multi-sample designs

- **Multi-sample / replicate-aware group comparisons** → **sample × cell-type pseudobulk** + **DESeq2 / edgeR + FDR**.
- **Forbidden** for between-group conclusions: direct **`FindMarkers`**, **`sc.tl.rank_genes_groups(groupby=condition)`**, or any **cell-level Wilcoxon / MAST** on the condition column when `n_replicates ≥ 2`.
- Exploratory **cluster markers** (Wilcoxon on `leiden`) remain allowed; label them exploratory only.

## 3. All figures must be publication-ready

- Export **≥ 300 dpi** (PNG and/or PDF).
- Use a **consistent font family and size** across figures in one run.
- Prefer **colorblind-friendly** palettes (e.g. Okabe–Ito); avoid red–green-only encoding for critical comparisons.
- Every saved figure needs a **caption** describing what is shown and what it is *not* evidence for.

## 4. Every step must record provenance

- Persist **software versions**, **parameters**, and **input/output paths** for each phase.
- Write `outputs/memory.yaml` (step timeline) and `run_manifest.json` (execution bundle) when code runs.
- Do not discard intermediate `.h5ad` checkpoints without documenting why.

---

*Prefer these four rules over long per-turn prompts. Load via `agents.common.system_rules()`; enforced in Reviewer + graph routing.*
