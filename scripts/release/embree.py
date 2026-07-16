"""Download and expose a pinned Embree binary release to subsequent CI steps."""

from __future__ import annotations

import json
import os
import re
import shutil
import stat
import tarfile
import zipfile
from pathlib import Path

from .common import ReleaseError, append_github_file, run


def select_asset(asset_names: list[str], pattern: str) -> str:
    """Return the first release asset matching a case-insensitive regular expression."""
    matcher = re.compile(pattern, re.IGNORECASE)
    for name in asset_names:
        if matcher.search(name):
            return name
    available = "\n  ".join(asset_names)
    raise ReleaseError(
        f"No Embree asset matched regex {pattern!r}. Available assets:\n  {available}"
    )


def _extract(archive: Path, destination: Path) -> None:
    """Extract a supported Embree ZIP or compressed tar archive."""
    name = archive.name.lower()
    if name.endswith(".zip"):
        _extract_zip(archive, destination)
    elif name.endswith((".tar.gz", ".tgz")):
        with tarfile.open(archive, "r:gz") as bundle:
            bundle.extractall(destination, filter="data")
    else:
        raise ReleaseError(f"Unsupported Embree archive format: {archive}")


def _extract_zip(archive: Path, destination: Path) -> None:
    """Extract a trusted release ZIP while retaining Unix dylib symlinks."""
    destination = destination.resolve()
    with zipfile.ZipFile(archive) as bundle:
        for member in bundle.infolist():
            target = (destination / member.filename).resolve()
            if target != destination and destination not in target.parents:
                raise ReleaseError(f"Archive member escapes destination: {member.filename}")
            if member.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue

            target.parent.mkdir(parents=True, exist_ok=True)
            mode = member.external_attr >> 16
            if stat.S_ISLNK(mode):
                link_target = bundle.read(member).decode("utf-8")
                target.symlink_to(link_target)
            else:
                with bundle.open(member) as source, target.open("wb") as output:
                    shutil.copyfileobj(source, output)
                if mode:
                    target.chmod(mode & 0o777)


def _find_cmake_dir(install_dir: Path, version: str) -> Path:
    """Locate the versioned CMake package directory required by ``find_package``."""
    expected = f"embree-{version}"
    matches = sorted(
        path
        for path in install_dir.rglob(expected)
        if path.is_dir()
        and path.parent.name == "cmake"
        and path.parent.parent.name in {"lib", "lib64"}
    )
    if not matches:
        raise ReleaseError(
            f"Could not find lib[/64]/cmake/{expected} below {install_dir}"
        )
    return matches[0]


def _find_runtime_dir(install_dir: Path) -> Path:
    """Locate the directory containing Embree's platform runtime library."""
    windows = os.name == "nt"
    for path in sorted(install_dir.rglob("*")):
        if not path.is_file():
            continue
        name = path.name.lower()
        if windows and name.startswith("embree") and name.endswith(".dll"):
            return path.parent
        if not windows and name.startswith("libembree") and (
            ".so" in name or name.endswith(".dylib")
        ):
            return path.parent
    raise ReleaseError(f"Could not find an Embree runtime library below {install_dir}")


def install(*, version: str, asset_regex: str, runner_temp: Path, repository: str) -> None:
    """Download Embree and publish its CMake and runtime paths to GitHub Actions.

    The selected archive is unpacked below ``runner_temp``. The function writes
    ``embree_DIR``, ``EMBREE_RUNTIME_DIR``, and ``EMBREE_INSTALL_DIR`` to
    ``GITHUB_ENV`` and adds the runtime directory to ``GITHUB_PATH``.
    """
    release = f"v{version}"
    archive_dir = runner_temp / "embree-archive"
    install_dir = runner_temp / "embree"
    shutil.rmtree(archive_dir, ignore_errors=True)
    shutil.rmtree(install_dir, ignore_errors=True)
    archive_dir.mkdir(parents=True)
    install_dir.mkdir(parents=True)

    metadata = run(
        ["gh", "release", "view", release, "--repo", repository, "--json", "assets"],
        capture=True,
    )
    assets = [asset["name"] for asset in json.loads(metadata.stdout)["assets"]]
    print("Available Embree release assets:")
    for name in assets:
        print(f"  {name}")

    asset_name = select_asset(assets, asset_regex)
    run(
        [
            "gh",
            "release",
            "download",
            release,
            "--repo",
            repository,
            "--pattern",
            asset_name,
            "--dir",
            archive_dir,
        ]
    )
    _extract(archive_dir / asset_name, install_dir)

    cmake_dir = _find_cmake_dir(install_dir, version)
    runtime_dir = _find_runtime_dir(install_dir)
    append_github_file(
        "GITHUB_ENV",
        [
            f"embree_DIR={cmake_dir}",
            f"EMBREE_RUNTIME_DIR={runtime_dir}",
            f"EMBREE_INSTALL_DIR={install_dir}",
        ],
    )
    append_github_file("GITHUB_PATH", [runtime_dir])
    print(f"Embree CMake package: {cmake_dir}")
    print(f"Embree runtime directory: {runtime_dir}")
