"""Construct and statically validate Linux AppImage release artifacts."""

from __future__ import annotations

import os
import re
import shutil
import stat
import urllib.request
from pathlib import Path

from .common import ReleaseError, enabled, require_path, run
from .licenses import bundle as bundle_licenses


LINUXDEPLOY_URL = (
    "https://github.com/linuxdeploy/linuxdeploy/releases/download/continuous/"
    "linuxdeploy-x86_64.AppImage"
)
QT_PLUGIN_URL = (
    "https://github.com/linuxdeploy/linuxdeploy-plugin-qt/releases/download/continuous/"
    "linuxdeploy-plugin-qt-x86_64.AppImage"
)
GLIBC_PATTERN = re.compile(rb"GLIBC_[0-9]+(?:\.[0-9]+)*")


def glibc_versions(data: bytes) -> set[str]:
    """Extract referenced GLIBC symbol versions from bytes belonging to an ELF file."""
    return {match.decode("ascii") for match in GLIBC_PATTERN.findall(data)}


def _download(url: str, destination: Path) -> None:
    """Download an AppImage tool and mark it executable."""
    print(f"Downloading {url}")
    with urllib.request.urlopen(url) as response, destination.open("wb") as output:
        shutil.copyfileobj(response, output)
    destination.chmod(destination.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)


def _runtime_environment() -> dict[str, str]:
    """Build the environment linuxdeploy needs to discover Embree and CUDA."""
    environment = os.environ.copy()
    library_paths: list[str] = []
    embree_runtime = environment.get("EMBREE_RUNTIME_DIR")
    if embree_runtime:
        library_paths.append(embree_runtime)
    cuda_path = environment.get("CUDA_PATH")
    if cuda_path and (Path(cuda_path) / "lib64").is_dir():
        library_paths.append(os.fspath(Path(cuda_path) / "lib64"))
    existing = environment.get("LD_LIBRARY_PATH")
    if existing:
        library_paths.append(existing)
    environment["LD_LIBRARY_PATH"] = os.pathsep.join(library_paths)
    environment["APPIMAGE_EXTRACT_AND_RUN"] = "1"
    return environment


def _prepare_appdir(
    *,
    workspace: Path,
    build_dir: Path,
    install_dir: Path,
    appdir: Path,
    app_name: str,
    optix_enabled: str | bool,
) -> tuple[Path, Path, Path]:
    """Create the AppDir skeleton and return its executable, desktop file, and icon."""
    shutil.rmtree(appdir, ignore_errors=True)
    executable = appdir / "usr" / "bin" / app_name
    desktop = appdir / "usr" / "share" / "applications" / f"{app_name}.desktop"
    icon = appdir / "usr" / "share" / "icons" / "hicolor" / "256x256" / "apps" / f"{app_name}.png"
    executable.parent.mkdir(parents=True)
    desktop.parent.mkdir(parents=True)
    icon.parent.mkdir(parents=True)

    shutil.copy2(
        require_path(install_dir / "bin" / app_name, "installed Linux executable"),
        executable,
    )
    shutil.copy2(
        require_path(workspace / "packaging" / "linux" / "SolTrace.desktop", "desktop file"),
        desktop,
    )
    run(
        [
            "convert",
            workspace / "gui" / "appicon" / "soltrace_icon.png",
            "-resize",
            "256x256",
            icon,
        ]
    )

    if enabled(optix_enabled):
        ptx_source = require_path(install_dir / "bin" / "ptx", "installed PTX directory")
        shutil.copytree(ptx_source, executable.parent / "ptx")

    bundle_licenses(
        platform="linux",
        workspace=workspace,
        build_dir=build_dir,
        install_dir=install_dir,
        app_name=app_name,
        optix_enabled=optix_enabled,
        appdir=appdir,
    )
    return executable, desktop, icon


def _validate_extracted(
    *,
    workspace: Path,
    extracted: Path,
    app_name: str,
    optix_enabled: str | bool,
    environment: dict[str, str],
) -> None:
    """Verify extracted files, dynamic dependencies, PTX assets, and glibc baseline."""
    executable = require_path(extracted / "usr" / "bin" / app_name, "extracted AppImage executable")
    if enabled(optix_enabled):
        for shader in ("intersection.ptx", "materials.ptx", "sun.ptx"):
            require_path(executable.parent / "ptx" / shader, f"PTX shader {shader}")

    validation_environment = environment.copy()
    # Do not let build-time Embree/CUDA paths conceal a missing bundled library.
    validation_environment.pop("LD_LIBRARY_PATH", None)
    dependencies = run(["ldd", executable], capture=True, env=validation_environment)
    (workspace / "appimage-ldd.txt").write_text(dependencies.stdout, encoding="utf-8")
    print(dependencies.stdout)
    if "not found" in dependencies.stdout:
        raise ReleaseError("AppImage contains unresolved runtime libraries")

    versions: set[str] = set()
    for candidate in extracted.rglob("*"):
        if not candidate.is_file():
            continue
        with candidate.open("rb") as stream:
            if stream.read(4) != b"\x7fELF":
                continue
            versions.update(glibc_versions(stream.read()))
    ordered = sorted(
        versions,
        key=lambda version: tuple(int(part) for part in version[6:].split(".")),
    )
    (workspace / "appimage-glibc-versions.txt").write_text(
        "\n".join(ordered) + "\n", encoding="utf-8"
    )
    print(f"Highest required glibc symbol version: {ordered[-1] if ordered else 'none found'}")


def package_appimage(
    *,
    workspace: Path,
    build_dir: Path,
    install_dir: Path,
    asset: Path,
    app_name: str,
    optix_enabled: str | bool,
) -> None:
    """Build an AppImage from the CMake install tree and validate its contents.

    linuxdeploy and its Qt plugin are downloaded from their continuous releases.
    The completed AppImage is extracted again so validation operates on the
    shipped filesystem rather than on build-tree files.
    """
    appdir = workspace / "AppDir"
    tools = workspace / "appimage-tools"
    extracted = workspace / "squashfs-root"
    shutil.rmtree(tools, ignore_errors=True)
    shutil.rmtree(extracted, ignore_errors=True)
    tools.mkdir(parents=True)
    asset.unlink(missing_ok=True)

    executable, desktop, icon = _prepare_appdir(
        workspace=workspace,
        build_dir=build_dir,
        install_dir=install_dir,
        appdir=appdir,
        app_name=app_name,
        optix_enabled=optix_enabled,
    )
    linuxdeploy = tools / "linuxdeploy-x86_64.AppImage"
    qt_plugin = tools / "linuxdeploy-plugin-qt-x86_64.AppImage"
    _download(LINUXDEPLOY_URL, linuxdeploy)
    _download(QT_PLUGIN_URL, qt_plugin)

    environment = _runtime_environment()
    environment["QML_SOURCES_PATHS"] = os.fspath(workspace / "gui" / "ui")
    before = set(workspace.glob("*.AppImage"))
    run(
        [
            linuxdeploy,
            "--appdir",
            appdir,
            "--executable",
            executable,
            "--desktop-file",
            desktop,
            "--icon-file",
            icon,
            "--plugin",
            "qt",
            "--output",
            "appimage",
        ],
        cwd=workspace,
        env=environment,
    )
    generated = set(workspace.glob("*.AppImage")) - before
    if len(generated) != 1:
        raise ReleaseError(f"Expected one generated AppImage, found: {sorted(generated)}")
    shutil.move(next(iter(generated)), asset)

    run([asset, "--appimage-extract"], cwd=workspace, env=environment, capture=True)
    _validate_extracted(
        workspace=workspace,
        extracted=extracted,
        app_name=app_name,
        optix_enabled=optix_enabled,
        environment=environment,
    )
