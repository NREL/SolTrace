"""Generate the GUI-embedded third-party license manifest."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path


LICENSE_PREFIXES = ("license", "copying", "notice", "eula", "third-party")


def is_license_file(path: Path) -> bool:
    return path.is_file() and path.name.lower().startswith(LICENSE_PREFIXES)


def read_text(path: Path) -> str:
    return path.read_text(encoding="utf-8", errors="replace")


def add_file(entries: list[dict[str, str]], component: str, root: Path, path: Path) -> None:
    entries.append(
        {
            "component": component,
            "path": path.relative_to(root).as_posix(),
            "name": path.name,
            "text": read_text(path),
        }
    )


def add_discovered(entries: list[dict[str, str]], component: str, root: Path | None) -> int:
    if root is None or not root.is_dir():
        return 0

    root = root.resolve()
    count = 0
    for path in sorted(root.rglob("*")):
        if is_license_file(path):
            add_file(entries, component, root, path.resolve())
            count += 1
    return count


def resolve_optional_path(value: str | None) -> Path | None:
    if not value:
        return None
    path = Path(value)
    return path.resolve() if path.exists() else None


def qt_license_dir(qt_dir: Path | None) -> Path | None:
    if qt_dir is None:
        return None

    candidates: list[Path] = []
    for base in [qt_dir, *qt_dir.parents]:
        candidates.append(base / "LICENSES")
        candidates.append(base / "Licenses")
    for candidate in candidates:
        if candidate.is_dir():
            return candidate.resolve()
    return None


def inferred_embree_roots(embree_dir: Path | None) -> list[Path]:
    roots: list[Path] = []
    if embree_dir is not None:
        roots.append(embree_dir)
        for parent in embree_dir.parents:
            roots.append(parent)

    if embree_dir is not None:
        for parent in embree_dir.parents:
            if parent == Path("/opt/homebrew"):
                roots.append(parent / "opt" / "embree")

    for candidate in (
        os.environ.get("EMBREE_INSTALL_DIR"),
        os.environ.get("embree_DIR"),
    ):
        path = resolve_optional_path(candidate)
        if path is not None:
            roots.append(path)

    seen: set[Path] = set()
    unique_roots: list[Path] = []
    for root in roots:
        resolved = root.resolve()
        if resolved not in seen:
            seen.add(resolved)
            unique_roots.append(resolved)
    return unique_roots


def add_first_available(
    entries: list[dict[str, str]], component: str, roots: list[Path]
) -> int:
    for root in roots:
        count = add_discovered(entries, component, root)
        if count > 0:
            return count
    return 0


def parser() -> argparse.ArgumentParser:
    result = argparse.ArgumentParser(description=__doc__)
    result.add_argument("--workspace", required=True, type=Path)
    result.add_argument("--build-dir", required=True, type=Path)
    result.add_argument("--qt-dir")
    result.add_argument("--embree-dir")
    result.add_argument("--optix-dir")
    result.add_argument("--cuda-dir")
    result.add_argument("--output", required=True, type=Path)
    return result


def main() -> int:
    args = parser().parse_args()
    workspace = args.workspace.resolve()
    build_dir = args.build_dir.resolve()
    output = args.output.resolve()

    entries: list[dict[str, str]] = []
    warnings: list[str] = []

    soltrace_license = workspace / "LICENSE.md"
    if soltrace_license.is_file():
        add_file(entries, "SolTrace", workspace, soltrace_license)
    else:
        warnings.append("SolTrace license file was not found.")

    qt_root = resolve_optional_path(os.environ.get("QT_ROOT_DIR")) or resolve_optional_path(args.qt_dir)
    qt_licenses = qt_license_dir(qt_root)
    if add_discovered(entries, "Qt", qt_licenses) == 0:
        warnings.append("Qt license files were not found in the local Qt installation.")

    if add_discovered(entries, "Build dependencies", build_dir / "_deps") == 0:
        warnings.append("Build dependency license files were not found.")

    if args.embree_dir is not None:
        embree_count = add_first_available(
            entries, "Embree", inferred_embree_roots(resolve_optional_path(args.embree_dir))
        )
        if embree_count == 0:
            warnings.append("Embree was configured, but license files were not found.")

    if args.optix_dir is not None:
        optix_root = resolve_optional_path(os.environ.get("OptiX_INSTALL_DIR")) or resolve_optional_path(args.optix_dir)
        if optix_root is not None and add_discovered(entries, "OptiX", optix_root) == 0:
            warnings.append("OptiX was configured, but license files were not found.")

    if args.cuda_dir is not None:
        cuda_root = resolve_optional_path(os.environ.get("CUDA_PATH")) or resolve_optional_path(args.cuda_dir)
        if cuda_root is not None and add_discovered(entries, "CUDA", cuda_root) == 0:
            warnings.append("CUDA was configured, but license files were not found.")

    entries.sort(key=lambda item: (item["component"].lower(), item["path"].lower()))
    manifest = {
        "schema": 1,
        "licenses": entries,
        "warnings": warnings,
    }

    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(manifest, indent=2, sort_keys=True), encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
