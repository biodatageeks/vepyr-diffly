from __future__ import annotations

import os
from pathlib import Path

from dotenv import load_dotenv


ROOT = Path(__file__).resolve().parents[2]
DOTENV_PATH = ROOT / ".env"


def load_repo_env() -> None:
    load_dotenv(DOTENV_PATH, override=False)


def env_str(name: str) -> str | None:
    value = os.getenv(name)
    if value is None:
        return None
    stripped = value.strip()
    return stripped or None


def env_path(name: str) -> Path | None:
    value = env_str(name)
    return None if value is None else Path(value)


def env_int(name: str) -> int | None:
    value = env_str(name)
    return None if value is None else int(value)


def env_int_or(name: str, default: int) -> int:
    value = env_int(name)
    return default if value is None else value


def env_bool(name: str) -> bool | None:
    value = env_str(name)
    if value is None:
        return None
    normalized = value.lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"invalid boolean value for {name}: {value}")
