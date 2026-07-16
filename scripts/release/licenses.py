"""Discover and bundle project and third-party license material."""

from __future__ import annotations

import json
import os
import shutil
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Mapping

from .common import ReleaseError, enabled, require_path


LICENSE_PREFIXES = ("license", "copying", "notice", "eula", "third-party")


def is_license_file(path: Path) -> bool:
    """Return whether a file name uses one of the recognized license prefixes."""
    return path.is_file() and path.name.lower().startswith(LICENSE_PREFIXES)


def safe_license_name(root: Path, path: Path, prefix: str) -> str:
    """Flatten a source-relative license path into a collision-resistant file name."""
    relative = path.relative_to(root)
    flattened = "__".join(relative.parts)
    return f"{prefix}__{flattened}"


def copy_discovered(root: Path | None, destination: Path, prefix: str) -> int:
    """Copy recognized license files below ``root`` and return the copied count."""
    if root is None or not root.is_dir():
        return 0
    count = 0
    for source in sorted(root.rglob("*")):
        if not is_license_file(source):
            continue
        shutil.copy2(source, destination / safe_license_name(root, source, prefix))
        count += 1
    return count


def find_qt_root(environment: Mapping[str, str] = os.environ) -> Path:
    """Resolve the Qt installation root from variables set by install-qt-action.

    aqt binary packages do not consistently contain a ``LICENSES`` directory,
    so locating the installation and locating license texts are separate tasks.
    """
    configured = environment.get("QT_ROOT_DIR")
    if configured and Path(configured).is_dir():
        return Path(configured)

    qt6_dir = environment.get("Qt6_DIR")
    if qt6_dir:
        path = Path(qt6_dir).resolve()
        if len(path.parents) >= 3:
            root = path.parents[2]
            if root.is_dir():
                return root
    raise ReleaseError("Could not locate the Qt installation")


def find_qt_license_dir(qt_root: Path) -> Path | None:
    """Find license texts shipped in or immediately above an aqt installation."""
    candidates = [qt_root / "LICENSES"]
    candidates.extend(parent / "LICENSES" for parent in list(qt_root.parents)[:3])
    for candidate in candidates:
        if candidate.is_dir():
            return candidate
    for candidate in qt_root.rglob("LICENSES"):
        if candidate.is_dir():
            return candidate
    return None


def download_qt_licenses(
    destination: Path,
    version: str,
    *,
    opener=urllib.request.urlopen,
) -> int:
    """Download canonical license texts from Qt's official GitHub mirror.

    Qt's aqt binary archives can omit their source ``LICENSES`` directory. The
    GitHub contents API is used to enumerate the exact files at the pinned Qt
    tag rather than maintaining a potentially stale license-name list here.
    """
    tag = f"v{version}"
    query = urllib.parse.urlencode({"ref": tag})
    api_url = f"https://api.github.com/repos/qt/qtbase/contents/LICENSES?{query}"
    headers = {
        "Accept": "application/vnd.github+json",
        "User-Agent": "SolTrace-release-tooling",
        "X-GitHub-Api-Version": "2022-11-28",
    }
    token = os.environ.get("GH_TOKEN") or os.environ.get("GITHUB_TOKEN")
    if token:
        headers["Authorization"] = f"Bearer {token}"

    try:
        with opener(urllib.request.Request(api_url, headers=headers)) as response:
            entries = json.load(response)
        files = [entry for entry in entries if entry.get("type") == "file"]
        if not files:
            raise ReleaseError(f"Qt {version} returned no canonical license files")
        destination.mkdir(parents=True, exist_ok=True)
        for entry in files:
            request = urllib.request.Request(entry["download_url"], headers=headers)
            with opener(request) as response, (destination / entry["name"]).open("wb") as output:
                shutil.copyfileobj(response, output)
    except (urllib.error.URLError, urllib.error.HTTPError, KeyError, ValueError) as error:
        raise ReleaseError(f"Could not download Qt {version} license texts: {error}") from error
    return len(files)


def destination_for(platform: str, install_dir: Path, app_name: str, appdir: Path | None) -> Path:
    """Choose the platform-native license directory inside a staged artifact."""
    if platform == "macos":
        return install_dir / f"{app_name}.app" / "Contents" / "Resources" / "licenses"
    if platform == "windows":
        return install_dir / "licenses"
    if platform == "linux" and appdir is not None:
        return appdir / "usr" / "share" / "licenses"
    raise ReleaseError(f"Unsupported license destination for platform {platform!r}")


def bundle(
    *,
    platform: str,
    workspace: Path,
    build_dir: Path,
    install_dir: Path,
    app_name: str,
    optix_enabled: str | bool,
    appdir: Path | None = None,
) -> Path:
    """Collect SolTrace, Qt, Embree, build dependency, CUDA, and OptiX licenses.

    Linux callers provide ``appdir`` because licenses are inserted before
    linuxdeploy creates the AppImage. macOS and Windows destinations are derived
    directly from the CMake installation tree.
    """
    destination = destination_for(platform, install_dir, app_name, appdir)
    soltrace_dir = destination / "SolTrace"
    soltrace_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(require_path(workspace / "LICENSE.md", "SolTrace license"), soltrace_dir)

    qt_root = find_qt_root()
    qt_destination = destination / "Qt"
    qt_license_dir = find_qt_license_dir(qt_root)
    if qt_license_dir:
        shutil.copytree(qt_license_dir, qt_destination, dirs_exist_ok=True)
        print(f"Copied Qt license files from {qt_license_dir}")
    else:
        qt_version = os.environ.get("QT_VERSION")
        if not qt_version:
            raise ReleaseError("QT_VERSION is required to retrieve Qt license texts")
        count = download_qt_licenses(qt_destination, qt_version)
        print(f"Downloaded {count} canonical Qt {qt_version} license file(s)")

    sources: list[tuple[Path | None, str]] = [
        (_environment_path("EMBREE_INSTALL_DIR"), "Embree"),
        (build_dir / "_deps", "Dependencies"),
    ]
    if enabled(optix_enabled):
        sources.extend(
            [
                (workspace / "optix-dev", "OptiX"),
                (_environment_path("CUDA_PATH"), "CUDA"),
            ]
        )

    for root, prefix in sources:
        count = copy_discovered(root, soltrace_dir, prefix)
        print(f"Copied {count} {prefix} license file(s)")
    return destination


def _environment_path(name: str) -> Path | None:
    """Convert an optional environment variable into a path."""
    value = os.environ.get(name)
    return Path(value) if value else None
