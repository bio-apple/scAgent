# knowledge/sops

实验室标准作业程序（SOP）。用下面命令加入 RAG，Planner / QC / 注释会检索这个 collection：

```bash
python -m scagent add-doc /path/to/lab_qc.md
python -m scagent add-doc ./sops/ --name our-lab
```

支持 `.md` / `.txt` / `.pdf` / `.ipynb`。不要用 `README.md` 当 SOP 文件名（ingest 会跳过）。

社区书用 `scagent update-kb`（写入 `knowledge/upstream/`），不要把 sc-best-practices 整本拷到这里。
