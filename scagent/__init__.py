"""scAgent — single-cell RNA-seq analysis agent."""

__version__ = "0.1.0"

__all__ = ["__version__", "load_config", "analysis_params", "read_single_cell"]


def __getattr__(name: str):
    if name in {"load_config", "analysis_params"}:
        from scagent.config import analysis_params, load_config

        return {"load_config": load_config, "analysis_params": analysis_params}[name]
    if name == "read_single_cell":
        from scagent.io import read_single_cell

        return read_single_cell
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
