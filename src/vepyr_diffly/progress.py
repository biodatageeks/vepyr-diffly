from __future__ import annotations

from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from threading import Event, Lock, Thread
from typing import Iterable

try:
    from rich.console import Console
except ModuleNotFoundError:  # pragma: no cover

    class Console:  # type: ignore[override]
        def print(self, value: object = "") -> None:
            print(value)


def _format_size(num_bytes: int) -> str:
    value = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if value < 1024.0 or unit == "TB":
            return f"{value:.1f}{unit}"
        value /= 1024.0
    return f"{num_bytes}B"


@dataclass
class ProgressReporter:
    log_path: Path
    console: Console | None = None
    heartbeat_seconds: float = 15.0
    _current_stage: str = field(init=False, default="idle")
    _tracked_paths: tuple[Path, ...] = field(init=False, default_factory=tuple)
    _stop_event: Event = field(init=False, default_factory=Event)
    _lock: Lock = field(init=False, default_factory=Lock)
    _heartbeat_thread: Thread | None = field(init=False, default=None)

    def start(self) -> None:
        if self._heartbeat_thread is not None:
            return
        self._heartbeat_thread = Thread(target=self._heartbeat_loop, daemon=True)
        self._heartbeat_thread.start()

    def stop(self) -> None:
        self._stop_event.set()
        if self._heartbeat_thread is not None:
            self._heartbeat_thread.join(timeout=1.0)
            self._heartbeat_thread = None

    def log(self, message: str) -> None:
        self._write_line(message)

    def stage(self, message: str, *, tracked_paths: Iterable[Path] = ()) -> None:
        with self._lock:
            self._current_stage = message
            self._tracked_paths = tuple(tracked_paths)
        self._write_line(message)

    def _heartbeat_loop(self) -> None:
        while not self._stop_event.wait(self.heartbeat_seconds):
            with self._lock:
                stage = self._current_stage
                tracked_paths = self._tracked_paths
            suffix = self._tracked_paths_status(tracked_paths)
            self._write_line(f"heartbeat: {stage}{suffix}")

    def _tracked_paths_status(self, tracked_paths: tuple[Path, ...]) -> str:
        if not tracked_paths:
            return ""
        parts: list[str] = []
        for path in tracked_paths:
            if path.exists():
                if path.is_dir():
                    file_count, total_size = self._directory_status(path)
                    parts.append(f"{path.name}={_format_size(total_size)}/{file_count} files")
                else:
                    parts.append(f"{path.name}={_format_size(path.stat().st_size)}")
            else:
                parts.append(f"{path.name}=missing")
        return f" [{' '.join(parts)}]"

    def _directory_status(self, path: Path) -> tuple[int, int]:
        file_count = 0
        total_size = 0
        for entry in path.rglob("*"):
            if not entry.is_file():
                continue
            file_count += 1
            total_size += entry.stat().st_size
        return file_count, total_size

    def _write_line(self, message: str) -> None:
        timestamp = datetime.now().astimezone().isoformat(timespec="seconds")
        line = f"[{timestamp}] {message}"
        self.log_path.parent.mkdir(parents=True, exist_ok=True)
        with self.log_path.open("a", encoding="utf-8") as handle:
            handle.write(line + "\n")
        if self.console is not None:
            self.console.print(line)
