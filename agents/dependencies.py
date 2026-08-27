"""Explicit DAG of analysis steps. DEG requires a groupby/annotation column, not a hard-coded list."""

from __future__ import annotations

from collections import defaultdict, deque

# node -> prerequisites
DEPS: dict[str, tuple[str, ...]] = {
    "qc": (),
    "normalize": ("qc",),
    "hvg": ("normalize",),
    "pca": ("hvg",),
    "ambient_soupx": ("qc",),
    "ambient_decontx": ("qc",),
    "impute_magic": ("normalize",),
    "impute_alra": ("normalize",),
    "harmony": ("pca",),
    "scvi": ("pca",),
    "cca": ("pca",),
    "neighbors": ("pca",),
    "leiden": ("neighbors",),
    "umap": ("neighbors",),
    "annotate": ("leiden",),
    "pseudobulk_deg": ("annotate",),  # needs cell_type / leiden groupby
    "trajectory": ("neighbors",),
}

INTENT_GOALS: dict[str, tuple[str, ...]] = {
    "qc": ("qc",),
    "clustering": ("leiden", "umap"),
    "annotation": ("annotate",),
    "deg": ("pseudobulk_deg",),
    "trajectory": ("trajectory",),
}


def _neighbors_prereq(integrator: str | None) -> str:
    if integrator in {"harmony", "scvi", "cca"}:
        return integrator
    return "pca"


def expand_goals(goals: set[str], *, integrator: str | None = None) -> list[str]:
    """Topological expansion. neighbors wait on the integrator when one is selected."""
    deps = {k: list(v) for k, v in DEPS.items()}
    deps["neighbors"] = [_neighbors_prereq(integrator)]
    if integrator:
        goals = set(goals) | {integrator}
    needed: set[str] = set()

    def visit(node: str) -> None:
        if node in needed:
            return
        for p in deps.get(node, []):
            visit(p)
        needed.add(node)

    for g in goals:
        visit(g)
    indeg: dict[str, int] = {n: 0 for n in needed}
    children: dict[str, list[str]] = defaultdict(list)
    for n in needed:
        for p in deps.get(n, []):
            if p in needed:
                indeg[n] += 1
                children[p].append(n)
    q = deque(sorted(n for n, d in indeg.items() if d == 0))
    order: list[str] = []
    while q:
        n = q.popleft()
        order.append(n)
        for c in sorted(children[n]):
            indeg[c] -= 1
            if indeg[c] == 0:
                q.append(c)
    if len(order) != len(needed):
        return sorted(needed)
    return order


def resolve_route(
    intents: list[str] | None,
    *,
    integrator: str | None = None,
    imputation: str | None = None,
    ambient: str | None = None,
    r_degraded: bool = False,
) -> list[str]:
    if r_degraded:
        return ["plan_only"]
    goals: set[str] = set()
    for intent in intents or []:
        goals.update(INTENT_GOALS.get(intent) or ())
    if not goals:
        goals = {"qc", "leiden", "umap", "annotate"}
    if imputation and imputation not in {"none", "off"}:
        goals.add(f"impute_{imputation}")
    if ambient and ambient not in {"none", "off"}:
        goals.add(f"ambient_{ambient}")
    return expand_goals(goals, integrator=integrator)
