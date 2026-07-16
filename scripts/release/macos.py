"""Stage and validate relocatable dependencies in the macOS application bundle."""

from __future__ import annotations

import os
import shutil
from pathlib import Path

from .common import ReleaseError, require_path, run


EXPECTED_RPATH = "@executable_path/../Frameworks"


def parse_otool_dependencies(output: str) -> list[str]:
    """Parse unique dependency install names from ``otool -L`` output.

    Dependency rows are indented. Thin- and universal-binary section headers
    are not, so filtering on indentation avoids treating an architecture header's
    absolute input path as a runtime dependency.
    """
    dependencies: list[str] = []
    for line in output.splitlines():
        if not line or not line[0].isspace():
            continue
        stripped = line.strip()
        dependency = stripped.split(" (", 1)[0]
        if dependency and dependency not in dependencies:
            dependencies.append(dependency)
    return dependencies


def dependency_is_relocatable(dependency: str) -> bool:
    """Allow bundle-relative references and Apple-provided system libraries."""
    return dependency.startswith(
        ("@rpath/", "@loader_path/", "@executable_path/", "/System/Library/", "/usr/lib/")
    )


def bundled_dependency_path(*, dependency: str, binary: Path, app: Path) -> Path | None:
    """Map a bundle-relative Mach-O install name to its expected on-disk path.

    The release layout keeps all ``@rpath`` libraries in ``Contents/Frameworks``.
    Loader- and executable-relative dependencies are resolved from the Mach-O
    file and the app's main executable directory respectively. Apple system
    dependencies return ``None`` because they are not bundled.
    """
    # Qt deploys QML metadata under Resources/qml and places plugin binaries in
    # Contents/PlugIns, leaving symlinks in the QML tree. dyld resolves
    # @loader_path from the real Mach-O location, not the QML symlink location.
    real_binary = binary.resolve()
    prefixes = {
        "@rpath/": app / "Contents" / "Frameworks",
        "@loader_path/": real_binary.parent,
        "@executable_path/": app / "Contents" / "MacOS",
    }
    for prefix, base in prefixes.items():
        if dependency.startswith(prefix):
            relative = dependency.removeprefix(prefix)
            return Path(os.path.normpath(base / relative))
    return None


def _copy_entry(source: Path, destination: Path) -> None:
    """Copy a regular file or recreate a relative dylib symlink without dereferencing it."""
    target = destination / source.name
    if target.is_symlink() or target.exists():
        target.unlink()
    if source.is_symlink():
        target.symlink_to(os.readlink(source))
    else:
        shutil.copy2(source, target)


def stage_runtime(*, install_dir: Path, app_name: str) -> None:
    """Copy Embree/TBB dylib chains into the bundle and ensure its runtime RPATH."""
    runtime_value = os.environ.get("EMBREE_RUNTIME_DIR")
    if not runtime_value:
        raise ReleaseError("EMBREE_RUNTIME_DIR is not set")
    runtime_dir = require_path(Path(runtime_value), "Embree runtime directory")
    app = require_path(install_dir / f"{app_name}.app", "macOS application bundle")
    executable = require_path(app / "Contents" / "MacOS" / app_name, "application executable")
    frameworks = app / "Contents" / "Frameworks"
    frameworks.mkdir(parents=True, exist_ok=True)

    libraries = sorted(
        path
        for path in runtime_dir.iterdir()
        if path.name.startswith(("libembree", "libtbb")) and path.name.endswith(".dylib")
    )
    if not libraries:
        raise ReleaseError(f"No Embree/TBB dylibs were found in {runtime_dir}")
    for library in libraries:
        _copy_entry(library, frameworks)

    load_commands = run(["otool", "-l", executable], capture=True).stdout
    if EXPECTED_RPATH not in load_commands:
        run(["install_name_tool", "-add_rpath", EXPECTED_RPATH, executable])


def ad_hoc_sign(*, install_dir: Path, app_name: str) -> None:
    """Ad-hoc sign embedded Mach-O code inside-out, then verify the app bundle.

    This makes every staged executable carry a consistent signature after the
    release pipeline has modified the bundle. It does not provide a trusted
    developer identity or satisfy Gatekeeper's notarization policy.
    """
    app = require_path(install_dir / f"{app_name}.app", "macOS application bundle")

    macho_files: list[Path] = []
    for candidate in sorted(path for path in app.rglob("*") if path.is_file()):
        if candidate.is_symlink():
            continue
        if "Mach-O" in run(["file", candidate], capture=True).stdout:
            macho_files.append(candidate)

    if not macho_files:
        raise ReleaseError(f"No Mach-O files found to sign in {app}")

    # Sign code files first, followed by embedded bundles from deepest to
    # shallowest. Signing the outer app last seals the final nested signatures.
    embedded_bundles = sorted(
        (
            path
            for path in app.rglob("*")
            if path.is_dir()
            and path.suffix.lower()
            in {".app", ".appex", ".bundle", ".framework", ".xpc"}
        ),
        key=lambda path: (-len(path.parts), str(path)),
    )
    for target in [*macho_files, *embedded_bundles, app]:
        run(["codesign", "--force", "--sign", "-", "--timestamp=none", target])

    run(["codesign", "--verify", "--deep", "--strict", "--verbose=4", app])


def validate(*, install_dir: Path, app_name: str, architecture: str) -> None:
    """Audit architecture, RPATHs, install names, and required Embree/TBB files."""
    app = require_path(install_dir / f"{app_name}.app", "macOS application bundle")
    executable = require_path(app / "Contents" / "MacOS" / app_name, "application executable")
    frameworks = require_path(app / "Contents" / "Frameworks", "Frameworks directory")

    file_output = run(["file", executable], capture=True).stdout
    if architecture not in file_output:
        raise ReleaseError(f"Expected {architecture} executable, got: {file_output.strip()}")
    if EXPECTED_RPATH not in run(["otool", "-l", executable], capture=True).stdout:
        raise ReleaseError(f"Application does not contain required RPATH {EXPECTED_RPATH}")

    macho_count = 0
    for candidate in sorted(path for path in app.rglob("*") if path.is_file()):
        if "Mach-O" not in run(["file", candidate], capture=True).stdout:
            continue
        macho_count += 1
        output = run(["otool", "-L", candidate], capture=True).stdout
        print(f"Dependencies for {candidate}:\n{output}")
        for dependency in parse_otool_dependencies(output):
            if not dependency_is_relocatable(dependency):
                raise ReleaseError(
                    f"Non-relocatable dependency in {candidate}: {dependency}"
                )
            bundled_path = bundled_dependency_path(
                dependency=dependency,
                binary=candidate,
                app=app,
            )
            if bundled_path is not None and not bundled_path.exists():
                raise ReleaseError(
                    f"Missing bundled dependency in {candidate}: {dependency} "
                    f"(expected {bundled_path})"
                )
            name = Path(dependency).name
            if name.startswith(("libembree", "libtbb")) and not (frameworks / name).exists():
                raise ReleaseError(f"Bundled dependency is missing: {name}")

    if macho_count == 0:
        raise ReleaseError(f"No Mach-O files found in {app}")
    run(["vtool", "-show-build", executable])
