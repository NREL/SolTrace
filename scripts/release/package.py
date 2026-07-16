"""Create final macOS and Windows release archives from staged install trees."""

from __future__ import annotations

import zipfile
from pathlib import Path

from .common import ReleaseError, require_path, run


def package_macos(*, install_dir: Path, asset: Path, app_name: str) -> None:
    """Archive a macOS bundle with ``ditto`` so Apple metadata is preserved."""
    app = require_path(install_dir / f"{app_name}.app", "macOS application bundle")
    asset.parent.mkdir(parents=True, exist_ok=True)
    asset.unlink(missing_ok=True)
    run(["ditto", "-c", "-k", "--sequesterRsrc", "--keepParent", app, asset])


def package_windows(*, install_dir: Path, asset: Path, app_name: str) -> None:
    """Create a ZIP whose root matches the Windows CMake installation tree."""
    require_path(install_dir / "bin" / f"{app_name}.exe", "Windows executable")
    asset.parent.mkdir(parents=True, exist_ok=True)
    asset.unlink(missing_ok=True)
    with zipfile.ZipFile(asset, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for source in sorted(install_dir.rglob("*")):
            if source.is_file():
                archive.write(source, source.relative_to(install_dir))
    if not asset.is_file():
        raise ReleaseError(f"Windows archive was not created: {asset}")
