"""LangGraph checkpointer: official SqliteSaver/RedisSaver when installed, else sqlite+pickle."""

from __future__ import annotations

import pickle
import sqlite3
import threading
from collections import defaultdict
from pathlib import Path

from langgraph.checkpoint.memory import InMemorySaver, MemorySaver

from scagent.config import load_config, resolve_path
from scagent.logutil import get_logger

log = get_logger("checkpoint")

_SAVERS: dict[str, object] = {}
_SAVERS_LOCK = threading.Lock()


class SqlitePickleSaver(InMemorySaver):
    """Durable InMemorySaver snapshot in SQLite. Works without langgraph-checkpoint-sqlite."""

    def __init__(self, path: str | Path):
        super().__init__()
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._lock = threading.RLock()
        self._conn = sqlite3.connect(str(self.path), check_same_thread=False, timeout=60.0)
        try:
            self._conn.execute("PRAGMA journal_mode=WAL")
        except sqlite3.Error:
            pass
        self._conn.execute("CREATE TABLE IF NOT EXISTS lg_ckpt (k TEXT PRIMARY KEY, v BLOB)")
        self._conn.commit()
        self._load()

    def _load(self) -> None:
        with self._lock:
            row = self._conn.execute("SELECT v FROM lg_ckpt WHERE k='full'").fetchone()
            if not row:
                return
            data = pickle.loads(row[0])
            storage = defaultdict(lambda: defaultdict(dict))
            for tid, nsmap in (data.get("storage") or {}).items():
                for ns, ckpts in nsmap.items():
                    storage[tid][ns].update(ckpts)
            writes = defaultdict(dict)
            writes.update(data.get("writes") or {})
            self.storage = storage
            self.writes = writes
            self.blobs = dict(data.get("blobs") or {})

    def _flush(self) -> None:
        with self._lock:
            payload = pickle.dumps(
                {
                    "storage": {
                        tid: {ns: dict(ckpts) for ns, ckpts in nsmap.items()} for tid, nsmap in self.storage.items()
                    },
                    "writes": dict(self.writes),
                    "blobs": dict(self.blobs),
                }
            )
            self._conn.execute("INSERT OR REPLACE INTO lg_ckpt (k, v) VALUES ('full', ?)", (payload,))
            self._conn.commit()

    def put(self, config, checkpoint, metadata, new_versions):
        with self._lock:
            out = super().put(config, checkpoint, metadata, new_versions)
            self._flush()
            return out

    def put_writes(self, *args, **kwargs):
        with self._lock:
            out = super().put_writes(*args, **kwargs)
            self._flush()
            return out

    def delete_thread(self, thread_id: str) -> None:
        with self._lock:
            super().delete_thread(thread_id)
            self._flush()


def last_thread_path(cfg: dict | None = None) -> Path:
    cfg = cfg or load_config()
    cache = resolve_path(cfg, "cache")
    cache.mkdir(parents=True, exist_ok=True)
    return cache / "last_thread_id"


def remember_thread(thread_id: str, cfg: dict | None = None) -> None:
    last_thread_path(cfg).write_text(thread_id, encoding="utf-8")


def load_last_thread(cfg: dict | None = None) -> str | None:
    p = last_thread_path(cfg)
    if not p.exists():
        return None
    text = p.read_text(encoding="utf-8").strip()
    return text or None


def get_checkpointer(cfg: dict | None = None, *, backend: str | None = None):
    cfg = cfg or load_config()
    ck = cfg.get("checkpoint") or {}
    backend = (backend or ck.get("backend") or "sqlite").lower()
    if backend == "memory":
        return MemorySaver()
    if backend == "redis":
        url = ck.get("redis_url")
        cache_key = f"redis:{url}"
        with _SAVERS_LOCK:
            if cache_key in _SAVERS:
                return _SAVERS[cache_key]
        try:
            from langgraph.checkpoint.redis import RedisSaver

            log.info("checkpointer=RedisSaver")
            saver = RedisSaver.from_conn_string(url)
            with _SAVERS_LOCK:
                _SAVERS[cache_key] = saver
            return saver
        except Exception as exc:
            log.warning("RedisSaver unavailable (%s); sqlite fallback", exc)
    path = Path(ck.get("sqlite_path") or ".cache/checkpoints.sqlite")
    if not path.is_absolute():
        path = Path(cfg.get("_root") or ".") / path
    cache_key = f"sqlite:{path.resolve()}"
    with _SAVERS_LOCK:
        if cache_key in _SAVERS:
            return _SAVERS[cache_key]
    try:
        from langgraph.checkpoint.sqlite import SqliteSaver

        log.info("checkpointer=SqliteSaver path=%s", path)
        if hasattr(SqliteSaver, "from_conn_string"):
            saver = SqliteSaver.from_conn_string(str(path))
            if hasattr(saver, "setup"):
                try:
                    saver.setup()
                except Exception:
                    pass
        else:
            conn = sqlite3.connect(str(path), check_same_thread=False, timeout=60.0)
            saver = SqliteSaver(conn)
            if hasattr(saver, "setup"):
                saver.setup()
        with _SAVERS_LOCK:
            _SAVERS[cache_key] = saver
        return saver
    except Exception as exc:
        log.info("official SqliteSaver missing (%s); using sqlite pickle saver %s", exc, path)
        saver = SqlitePickleSaver(path)
        with _SAVERS_LOCK:
            _SAVERS[cache_key] = saver
        return saver
