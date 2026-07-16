"""Stage and statically validate Windows runtime dependencies."""

from __future__ import annotations

import os
import re
import shutil
from pathlib import Path

from .common import ReleaseError, enabled, require_path, run


INSTALLER_EXCLUDED_ROOTS = frozenset({"include", "lib"})
MSI_VERSION_PATTERN = re.compile(r"^[0-9]+\.[0-9]+\.[0-9]+$")
WIX_GUID_PATTERN = re.compile(
    r"^[0-9A-Fa-f]{8}-[0-9A-Fa-f]{4}-[0-9A-Fa-f]{4}-"
    r"[0-9A-Fa-f]{4}-[0-9A-Fa-f]{12}$"
)


def _copy_matches(source: Path, pattern: str, destination: Path, description: str) -> None:
    """Copy all matching files, failing when a required runtime family is absent."""
    matches = sorted(path for path in source.glob(pattern) if path.is_file())
    if not matches:
        raise ReleaseError(f"No {description} matching {pattern!r} found in {source}")
    destination.mkdir(parents=True, exist_ok=True)
    for path in matches:
        shutil.copy2(path, destination)


def stage_runtime(*, install_dir: Path, optix_enabled: str | bool) -> None:
    """Copy Embree/TBB and, for OptiX builds, CUDA runtime DLLs beside the app."""
    runtime_value = os.environ.get("EMBREE_RUNTIME_DIR")
    if not runtime_value:
        raise ReleaseError("EMBREE_RUNTIME_DIR is not set")
    bin_dir = install_dir / "bin"
    _copy_matches(Path(runtime_value), "*.dll", bin_dir, "Embree runtime DLLs")

    if enabled(optix_enabled):
        cuda_value = os.environ.get("CUDA_PATH")
        if not cuda_value:
            raise ReleaseError("CUDA_PATH is not set for an OptiX build")
        _copy_matches(Path(cuda_value) / "bin", "cudart64_*.dll", bin_dir, "CUDA runtime DLLs")


def _require_match(directory: Path, pattern: str, description: str) -> None:
    """Require at least one staged file matching a dependency pattern."""
    if not any(path.is_file() for path in directory.glob(pattern)):
        raise ReleaseError(f"{description} is missing from {directory} (pattern {pattern})")


def _find_dumpbin() -> Path:
    """Locate the x64 ``dumpbin`` belonging to the latest installed MSVC toolset."""
    program_files = os.environ.get("ProgramFiles(x86)")
    if not program_files:
        raise ReleaseError("ProgramFiles(x86) is not set")
    vswhere = require_path(
        Path(program_files) / "Microsoft Visual Studio" / "Installer" / "vswhere.exe",
        "vswhere.exe",
    )
    result = run(
        [
            vswhere,
            "-latest",
            "-products",
            "*",
            "-requires",
            "Microsoft.VisualStudio.Component.VC.Tools.x86.x64",
            "-property",
            "installationPath",
        ],
        capture=True,
    )
    visual_studio = Path(result.stdout.strip())
    if not visual_studio:
        raise ReleaseError("Visual Studio with x64 C++ tools was not found")
    version_file = require_path(
        visual_studio / "VC" / "Auxiliary" / "Build" / "Microsoft.VCToolsVersion.default.txt",
        "MSVC tools version file",
    )
    tools_version = version_file.read_text(encoding="utf-8").strip()
    return require_path(
        visual_studio
        / "VC"
        / "Tools"
        / "MSVC"
        / tools_version
        / "bin"
        / "Hostx64"
        / "x64"
        / "dumpbin.exe",
        "dumpbin.exe",
    )


def validate(*, install_dir: Path, app_name: str, optix_enabled: str | bool) -> None:
    """Require key DLL families and print the executable's PE dependency table."""
    bin_dir = require_path(install_dir / "bin", "Windows binary directory")
    executable = require_path(bin_dir / f"{app_name}.exe", "SolTrace executable")
    _require_match(bin_dir, "Qt6Core.dll", "Qt6Core.dll")
    _require_match(bin_dir, "embree*.dll", "Embree runtime")
    _require_match(bin_dir, "tbb*.dll", "TBB runtime")
    if enabled(optix_enabled):
        _require_match(bin_dir, "cudart64_*.dll", "CUDA runtime")
        ptx_dir = require_path(bin_dir / "ptx", "OptiX PTX directory")
        for shader in ("intersection.ptx", "materials.ptx", "sun.ptx"):
            require_path(ptx_dir / shader, f"OptiX PTX shader {shader}")
    run([_find_dumpbin(), "/DEPENDENTS", executable])


def prepare_installer_root(
    *, install_dir: Path, installer_root: Path, app_name: str
) -> None:
    """Create a consumer runtime tree without changing CMake install rules.

    The full install tree remains available for the portable ZIP and developer
    use. Top-level development directories are omitted only from the MSI input.
    """
    require_path(install_dir / "bin" / f"{app_name}.exe", "Windows executable")
    require_path(install_dir / "licenses", "bundled license directory")
    shutil.rmtree(installer_root, ignore_errors=True)
    installer_root.mkdir(parents=True)

    for source in sorted(install_dir.iterdir()):
        if source.name.lower() in INSTALLER_EXCLUDED_ROOTS:
            continue
        destination = installer_root / source.name
        if source.is_dir():
            shutil.copytree(source, destination, symlinks=True)
        else:
            shutil.copy2(source, destination)

    require_path(
        installer_root / "bin" / f"{app_name}.exe",
        "installer runtime executable",
    )


def _cmake_value(value: str | Path) -> str:
    """Escape a value for a double-quoted CMake string."""
    text = value.as_posix() if isinstance(value, Path) else value
    return text.replace("\\", "/").replace('"', '\\"')


def _render_cpack_config(template: Path, output: Path, values: dict[str, str | Path]) -> None:
    """Render the standalone WiX CPack template and reject unresolved tokens."""
    rendered = require_path(template, "WiX CPack template").read_text(encoding="utf-8")
    for key, value in values.items():
        rendered = rendered.replace(f"@{key}@", _cmake_value(value))
    unresolved = sorted(set(re.findall(r"@[A-Z0-9_]+@", rendered)))
    if unresolved:
        raise ReleaseError(f"Unresolved WiX CPack template values: {unresolved}")
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(rendered, encoding="utf-8")


def package_msi(
    *,
    workspace: Path,
    build_dir: Path,
    install_dir: Path,
    installer_root: Path,
    asset: Path,
    app_name: str,
    package_name: str,
    install_directory: str,
    package_version: str,
    upgrade_guid: str,
) -> None:
    """Package a filtered, pre-staged runtime tree with standalone CPack/WiX.

    CPack consumes ``installer_root`` through ``CPACK_INSTALLED_DIRECTORIES``;
    it does not execute this project's install rules or dependency discovery.
    """
    if not MSI_VERSION_PATTERN.fullmatch(package_version):
        raise ReleaseError(
            f"MSI package version must contain three numeric fields: {package_version!r}"
        )
    if not WIX_GUID_PATTERN.fullmatch(upgrade_guid):
        raise ReleaseError(f"MSI UpgradeCode is not a GUID: {upgrade_guid!r}")
    prepare_installer_root(
        install_dir=install_dir,
        installer_root=installer_root,
        app_name=app_name,
    )

    asset.parent.mkdir(parents=True, exist_ok=True)
    asset.unlink(missing_ok=True)
    config_dir = build_dir / "windows-msi"
    config = config_dir / f"CPackConfig-{asset.stem}.cmake"
    license_file = config_dir / "LICENSE.txt"
    config_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(require_path(workspace / "LICENSE.md", "SolTrace license"), license_file)
    template = workspace / "packaging" / "windows" / "CPackWixConfig.cmake.in"
    _render_cpack_config(
        template,
        config,
        {
            "APP_NAME": app_name,
            "INSTALL_DIRECTORY": install_directory,
            "INSTALLER_ROOT": installer_root,
            "LICENSE_FILE": license_file,
            "PACKAGE_DIRECTORY": asset.parent,
            "PACKAGE_FILE_NAME": asset.stem,
            "PACKAGE_NAME": package_name,
            "PACKAGE_VERSION": package_version,
            "UPGRADE_GUID": upgrade_guid.upper(),
        },
    )
    run(["cpack", "--config", config, "-G", "WIX", "-B", asset.parent])
    require_path(asset, "generated Windows MSI")
