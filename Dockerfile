FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential \
    && rm -rf /var/lib/apt/lists/*

COPY pyproject.toml README.md requirements.txt requirements-analysis.txt ./
COPY scagent scagent
COPY agents agents
COPY workflows workflows
COPY rag rag
COPY sandbox sandbox
COPY prompts prompts
COPY config.yaml ./
COPY knowledge knowledge
COPY report_templates report_templates
COPY tests tests

RUN pip install --no-cache-dir -r requirements.txt \
    && pip install --no-cache-dir -e ".[dev]"

ENV PYTHONUNBUFFERED=1 \
    SCAGENT_LOG_LEVEL=INFO

CMD ["python", "-m", "scagent", "--help"]
