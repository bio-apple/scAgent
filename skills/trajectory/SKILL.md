---
name: trajectory
description: "Pseudotime / lineage inference (PAGA, DPT, Palantir, Monocle3, scVelo when applicable)."
---

# Scientific task: Trajectory / fate

## Goal
Infer continuous structure or differentiation axes when biology is not purely discrete clusters.

## Methods
| Tool | Role |
|------|------|
| PAGA | Graph abstraction / continuity check |
| DPT / Palantir | Pseudotime on continuous manifolds |
| Monocle3 | R trajectory (when R path) |
| scVelo | RNA velocity only if spliced/unspliced present |

## Outputs
- Pseudotime values
- Trajectory embedding / PAGA graph
- Gene-vs-pseudotime trends for key genes
- Confidence / verdict (discrete vs continuous)

## Gates
- Do **not** force a fate axis on clearly discrete populations.
- Velocity without layers is invalid — skip and say so.

## Related
- `knowledge/best_practices/trajectory.md`

