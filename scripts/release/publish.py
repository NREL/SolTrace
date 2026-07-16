"""Publish one matrix artifact to the tag's shared draft GitHub release."""

from __future__ import annotations

from pathlib import Path

from .common import ReleaseError, require_path, run


def publish(*, tag: str, asset: Path) -> None:
    """Create the draft release if needed, then upload or replace one asset.

    Matrix jobs may attempt release creation concurrently. A failed create is
    accepted only when a follow-up lookup confirms that another job won the race.
    """
    require_path(asset, "release asset")
    existing = run(["gh", "release", "view", tag], check=False, capture=True)
    if existing.returncode != 0:
        created = run(
            [
                "gh",
                "release",
                "create",
                tag,
                "--draft",
                "--generate-notes",
                "--title",
                tag,
            ],
            check=False,
            capture=True,
        )
        if created.returncode != 0:
            # Matrix jobs can race while creating the shared draft release.
            raced = run(["gh", "release", "view", tag], check=False, capture=True)
            if raced.returncode != 0:
                raise ReleaseError(
                    f"Could not create or find release {tag}: {created.stderr.strip()}"
                )
    run(["gh", "release", "upload", tag, asset, "--clobber"])
