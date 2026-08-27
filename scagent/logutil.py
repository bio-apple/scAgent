from __future__ import annotations

import logging
import os
import time
from contextlib import contextmanager
from pathlib import Path

_LOGGER_NAME = "scagent"
_configured = False


def setup_logging(level: str | None = None, log_file: str | Path | None = None) -> logging.Logger:
    """Idempotent logging setup. Level/file come from config or env SCAGENT_LOG_LEVEL."""
    global _configured
    logger = logging.getLogger(_LOGGER_NAME)
    if _configured:
        if level:
            logger.setLevel(getattr(logging, str(level).upper(), logging.INFO))
        return logger
    raw = (level or os.getenv("SCAGENT_LOG_LEVEL") or "INFO").upper()
    logger.setLevel(getattr(logging, raw, logging.INFO))
    fmt = logging.Formatter("%(asctime)s %(levelname)s [%(name)s] %(message)s", datefmt="%H:%M:%S")
    sh = logging.StreamHandler()
    sh.setFormatter(fmt)
    logger.addHandler(sh)
    if log_file:
        path = Path(log_file)
        path.parent.mkdir(parents=True, exist_ok=True)
        fh = logging.FileHandler(path, encoding="utf-8")
        fh.setFormatter(fmt)
        logger.addHandler(fh)
    logger.propagate = False
    _configured = True
    return logger


def get_logger(name: str | None = None) -> logging.Logger:
    base = logging.getLogger(_LOGGER_NAME)
    if not base.handlers:
        setup_logging()
    return logging.getLogger(f"{_LOGGER_NAME}.{name}") if name else base


@contextmanager
def timed(label: str, logger: logging.Logger | None = None):
    log = logger or get_logger()
    t0 = time.perf_counter()
    log.info("%s: start", label)
    try:
        yield
    except Exception:
        log.exception("%s: failed after %.2fs", label, time.perf_counter() - t0)
        raise
    else:
        log.info("%s: done in %.2fs", label, time.perf_counter() - t0)
