# Luecken et al. 2022 — scIB: Benchmarking atlas-level data integration

Nature Methods. Benchmark of 68 integration methods at atlas scale.

## Key conclusions

- No universal integrator. Evaluate **batch removal** and **bio-conservation** together.
- HVG selection often helps integration; scaling can erase biology.
- Harmony, scVI, Scanorama win on different data shapes—choose by task.
- UMAP mixing is not success evidence. Use kBET, iLISI, cLISI, NMI/ARI, or marker separability to detect overcorrection.

## scAgent rules

- Multi-sample: is batch technical noise or biology (treatment, disease, donor)?
- Default Python: Harmony first; scVI for nonlinear batch or uncertainty needs.
- After integration: check whether types merge while condition differences disappear.

## Typical queries

“Harmony or scVI?”, “Should I integrate?”, “How to judge integration?”
