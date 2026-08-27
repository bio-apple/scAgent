from pathlib import Path

from agents.writer import render_report
from sandbox.executor import analysis_executor, write_and_maybe_run
from scagent.config import load_config
from scagent.dual import build_dual, export_dual, render_dual_markdown, strip_code_fences
from scagent.export_nb import build_notebook, execute_via_jupyter


def test_strip_fences_keeps_plain_python():
    assert strip_code_fences("print(1)\n") == "print(1)\n"
    fenced = "```python\nimport scanpy as sc\nprint(sc.__name__)\n```"
    assert "```" not in strip_code_fences(fenced)
    assert "import scanpy as sc" in strip_code_fences(fenced)


def test_dual_markdown_separates_conclusion_and_code(tmp_path):
    state = {
        "user_query": "QC",
        "report_lang": "zh",
        "qc_strategy": {"method": "mad", "pct_mt_note": "先看分布"},
        "code_qc": "import scanpy as sc\nprint('qc-ok')\n",
        "code_downstream": "print('down-ok')\n",
        "artifacts": {"metrics": {}},
    }
    md = render_dual_markdown(state)
    assert "## [结论] QC & Preprocessing" in md
    assert "## [代码] QC & Preprocessing" in md
    assert "## [结论] Clustering & Differential" in md
    assert "import scanpy as sc" in md
    assert md.index("[结论] QC") < md.index("[代码] QC")
    dual = build_dual(state)
    assert dual["format"] == "code-result-v1"
    assert dual["blocks"][0]["language"] == "python"
    path = export_dual(state, tmp_path)
    assert path.name == "dual.md"
    assert "[结论]" in path.read_text(encoding="utf-8")


def test_report_embeds_dual_blocks():
    report = render_report(
        {
            "user_query": "x",
            "report_lang": "zh",
            "code_qc": "print(1)\n",
            "plan": {"route": ["qc"]},
            "artifacts": {},
        }
    )
    assert "[结论]" in report
    assert "[代码]" in report
    assert "```python" in report


def test_notebook_cells_are_conclusion_then_code():
    nb = build_notebook(
        {
            "code_qc": "print('qc')\n",
            "code_downstream": "print('down')\n",
            "plan": {"route": ["qc"]},
            "artifacts": {"metrics": {"seed": 0}},
        }
    )
    types = [c["cell_type"] for c in nb["cells"]]
    text = [c["source"] for c in nb["cells"]]
    md_i = next(i for i, t in enumerate(text) if "[结论] QC" in t)
    code_i = next(i for i, (kind, src) in enumerate(zip(types, text)) if kind == "code" and "print('qc')" in src)
    assert types[md_i] == "markdown"
    assert "[代码]" in text[md_i + 1]
    assert code_i == md_i + 2


def test_default_executor_is_jupyter():
    assert analysis_executor(load_config()) == "jupyter"


def test_jupyter_executor_runs_print_and_writes_ipynb(tmp_path):
    r = write_and_maybe_run(
        "print('ok-jupyter')\n",
        workspace=tmp_path,
        execute=True,
        timeout=30,
        filename="analysis.py",
        cfg={**load_config(), "analysis": {**load_config().get("analysis", {}), "executor": "jupyter"}},
    )
    assert r["ok"] is True, r["stderr"]
    assert "ok-jupyter" in r["stdout"]
    assert r["jail"] in {"jupyter", "jupyter-subprocess"}
    assert (tmp_path / "analysis.py").is_file()
    assert (tmp_path / "analysis.ipynb").is_file()


def test_jupyter_still_blocks_os_system(tmp_path):
    r = write_and_maybe_run(
        "import os\nos.system('true')\n",
        workspace=tmp_path,
        execute=True,
        timeout=5,
        cfg={**load_config(), "analysis": {**load_config().get("analysis", {}), "executor": "jupyter"}},
    )
    assert r["ok"] is False
    assert r["executed"] is False
    assert r["jail"] == "policy"


def test_execute_via_jupyter_fallback(tmp_path):
    script = tmp_path / "cell.py"
    script.write_text("print('cell-ok')\n", encoding="utf-8")
    ran = execute_via_jupyter(
        "print('cell-ok')\n",
        workspace=tmp_path,
        timeout=20,
        filename="cell.py",
        script=script,
        env={"PYTHONHASHSEED": "0", **{k: v for k, v in __import__("os").environ.items()}},
    )
    assert ran["ok"] is True, ran.get("stderr")
    assert "cell-ok" in ran["stdout"]
    assert Path(ran["notebook"]).is_file()
