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

