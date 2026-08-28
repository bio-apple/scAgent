# QC decisions: data-driven, not fixed cutoffs

Required outputs:

1. Violin: `n_genes_by_counts`, `total_counts`, `pct_counts_mt` (plus ribo/hb when relevant)
2. Scatter: `total_counts` vs `n_genes_by_counts`; `total_counts` vs `pct_counts_mt`
3. MAD: `median ± n * MAD`, n typically 3–5, tissue-adjusted

Tissue priors: `config.yaml` → `qc_profiles`. Relax mito filtering for heart, tumor, metabolic studies (Yates 2025).
Log cells removed at each step. Empty droplets via barcode rank, not min_genes alone.
