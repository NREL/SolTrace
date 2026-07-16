"""Shared process, path, and GitHub Actions helpers for release commands."""

from __future__ import annotations

import os
import shlex
import subprocess
from pathlib import Path
from typing import Iterable, Mapping, Sequence


class ReleaseError(RuntimeError):
    """A release artifact could not be prepared or validated."""


def run(
    command: Sequence[str | os.PathLike[str]],
    *,
    cwd: Path | None = None,
    env: Mapping[str, str] | None = None,
    capture: bool = False,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a command, log it, and convert process failures into ``ReleaseError``."""
    args = [os.fspath(part) for part in command]
    print(f"+ {shlex.join(args)}", flush=True)
    try:
        return subprocess.run(
            args,
            cwd=cwd,
            env=env,
            check=check,
            text=True,
            stdout=subprocess.PIPE if capture else None,
            stderr=subprocess.PIPE if capture else None,
        )
    except FileNotFoundError as error:
        raise ReleaseError(f"Required command was not found: {args[0]}") from error
    except subprocess.CalledProcessError as error:
        details = (error.stderr or error.stdout or "").strip()
        suffix = f"\n{details}" if details else ""
        raise ReleaseError(
            f"Command failed with exit code {error.returncode}: {shlex.join(args)}{suffix}"
        ) from error


def require_path(path: Path, description: str) -> Path:
    """Return an existing path or raise an error that identifies the missing input."""
    if not path.exists():
        raise ReleaseError(f"{description} was not found: {path}")
    return path


def append_github_file(variable: str, values: Iterable[str | os.PathLike[str]]) -> None:
    """Append lines to a GitHub Actions environment or PATH command file."""
    output = os.environ.get(variable)
    if not output:
        raise ReleaseError(f"{variable} is not set; this command expects GitHub Actions")
    with Path(output).open("a", encoding="utf-8") as stream:
        for value in values:
            stream.write(f"{os.fspath(value)}\n")


def enabled(value: str | bool) -> bool:
    """Interpret common CMake and environment representations of a true value."""
    return value is True or str(value).upper() in {"1", "ON", "TRUE", "YES"}
