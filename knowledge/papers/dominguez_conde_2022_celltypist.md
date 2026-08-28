# Domínguez Conde et al. 2022 — Cross-tissue immune cell analysis (CellTypist)

Science. Logistic-regression reference annotation; multi-tissue immune models.

## Usage

- Input: log1p-normalized AnnData (~10^4 counts/cell).
- Output: per-cell labels + optional majority vote (Leiden-smoothed).
- Confidence < 0.5: manual review—doublet, transition state, or missing reference type.

## Layered evidence (no single-gene labels)

1. Unbiased clustering (Leiden; resolution from stability + marker interpretability).
2. Reference mapping (CellTypist / Azimuth / scANVI / popV).
3. ≥2 independent canonical positive markers + negative markers.
4. Naming consistent with tissue/disease context (HCA community names).

## Failure modes

- Tumor microenvironment, developmental mismatch, cross-species adult immune models applied blindly.
- popV (Luecken 2024) when ensemble uncertainty is needed.
