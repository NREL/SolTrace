"""Command-line entry point for SolTrace release packaging and validation.

The GitHub Actions workflow invokes this module with ``python -m``. Keeping the
command definitions here makes the workflow declarative while the platform
modules remain directly importable and unit-testable.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from . import embree, licenses, linux, macos, package, publish, windows
from .common import ReleaseError


def _path(value: str) -> Path:
    """Normalize a command-line path relative to the caller's working directory."""
    return Path(value).resolve()


def _parser() -> argparse.ArgumentParser:
    """Define the stable command-line interface used by the release workflow."""
    parser = argparse.ArgumentParser(description="Build and validate SolTrace release artifacts")
    commands = parser.add_subparsers(dest="command", required=True)

    install_embree = commands.add_parser(
        "install-embree", help="Download and expose an Embree release"
    )
    install_embree.add_argument("--version", required=True)
    install_embree.add_argument("--asset-regex", required=True)
    install_embree.add_argument("--runner-temp", required=True, type=_path)
    install_embree.add_argument("--repository", default="RenderKit/embree")

    stage = commands.add_parser(
        "stage-runtime", help="Copy non-Qt runtime libraries into an artifact"
    )
    stage.add_argument("--platform", required=True, choices=("macos", "windows"))
    stage.add_argument("--install-dir", required=True, type=_path)
    stage.add_argument("--app-name", required=True)
    stage.add_argument("--optix-enabled", default="OFF")

    sign_macos = commands.add_parser(
        "sign-macos", help="Ad-hoc sign and verify a staged macOS application"
    )
    sign_macos.add_argument("--install-dir", required=True, type=_path)
    sign_macos.add_argument("--app-name", required=True)

    license_parser = commands.add_parser(
        "bundle-licenses", help="Collect redistributable license files"
    )
    _add_artifact_arguments(license_parser, include_asset=False)

    validate = commands.add_parser("validate", help="Statically validate a staged artifact")
    validate.add_argument("--platform", required=True, choices=("macos", "windows"))
    validate.add_argument("--install-dir", required=True, type=_path)
    validate.add_argument("--app-name", required=True)
    validate.add_argument("--optix-enabled", default="OFF")
    validate.add_argument("--architecture", default="arm64")

    package_parser = commands.add_parser("package", help="Create the platform release archive")
    _add_artifact_arguments(package_parser, include_asset=True)

    msi_parser = commands.add_parser(
        "package-windows-msi",
        help="Create a WiX MSI from the staged Windows runtime tree",
    )
    msi_parser.add_argument("--workspace", required=True, type=_path)
    msi_parser.add_argument("--build-dir", required=True, type=_path)
    msi_parser.add_argument("--install-dir", required=True, type=_path)
    msi_parser.add_argument("--installer-root", required=True, type=_path)
    msi_parser.add_argument("--asset", required=True, type=_path)
    msi_parser.add_argument("--app-name", required=True)
    msi_parser.add_argument("--package-name", required=True)
    msi_parser.add_argument("--install-directory", required=True)
    msi_parser.add_argument("--package-version", required=True)
    msi_parser.add_argument("--upgrade-guid", required=True)

    publish_parser = commands.add_parser(
        "publish", help="Upload an asset to a draft GitHub release"
    )
    publish_parser.add_argument("--tag", required=True)
    publish_parser.add_argument("--asset", required=True, type=_path)
    return parser


def _add_artifact_arguments(parser: argparse.ArgumentParser, *, include_asset: bool) -> None:
    """Add paths and build switches shared by artifact-oriented commands."""
    parser.add_argument("--platform", required=True, choices=("macos", "windows", "linux"))
    parser.add_argument("--workspace", required=True, type=_path)
    parser.add_argument("--build-dir", required=True, type=_path)
    parser.add_argument("--install-dir", required=True, type=_path)
    parser.add_argument("--app-name", required=True)
    parser.add_argument("--optix-enabled", default="OFF")
    if include_asset:
        parser.add_argument("--asset", required=True, type=_path)


def _dispatch(args: argparse.Namespace) -> None:
    """Route parsed arguments to the appropriate platform implementation."""
    if args.command == "install-embree":
        embree.install(
            version=args.version,
            asset_regex=args.asset_regex,
            runner_temp=args.runner_temp,
            repository=args.repository,
        )
    elif args.command == "stage-runtime":
        if args.platform == "macos":
            macos.stage_runtime(install_dir=args.install_dir, app_name=args.app_name)
        else:
            windows.stage_runtime(install_dir=args.install_dir, optix_enabled=args.optix_enabled)
    elif args.command == "sign-macos":
        macos.ad_hoc_sign(install_dir=args.install_dir, app_name=args.app_name)
    elif args.command == "bundle-licenses":
        licenses.bundle(
            platform=args.platform,
            workspace=args.workspace,
            build_dir=args.build_dir,
            install_dir=args.install_dir,
            app_name=args.app_name,
            optix_enabled=args.optix_enabled,
        )
    elif args.command == "validate":
        if args.platform == "macos":
            macos.validate(
                install_dir=args.install_dir,
                app_name=args.app_name,
                architecture=args.architecture,
            )
        else:
            windows.validate(
                install_dir=args.install_dir,
                app_name=args.app_name,
                optix_enabled=args.optix_enabled,
            )
    elif args.command == "package":
        if args.platform == "macos":
            package.package_macos(
                install_dir=args.install_dir, asset=args.asset, app_name=args.app_name
            )
        elif args.platform == "windows":
            package.package_windows(
                install_dir=args.install_dir, asset=args.asset, app_name=args.app_name
            )
        else:
            linux.package_appimage(
                workspace=args.workspace,
                build_dir=args.build_dir,
                install_dir=args.install_dir,
                asset=args.asset,
                app_name=args.app_name,
                optix_enabled=args.optix_enabled,
            )
    elif args.command == "package-windows-msi":
        windows.package_msi(
            workspace=args.workspace,
            build_dir=args.build_dir,
            install_dir=args.install_dir,
            installer_root=args.installer_root,
            asset=args.asset,
            app_name=args.app_name,
            package_name=args.package_name,
            install_directory=args.install_directory,
            package_version=args.package_version,
            upgrade_guid=args.upgrade_guid,
        )
    elif args.command == "publish":
        publish.publish(tag=args.tag, asset=args.asset)
    else:
        raise ReleaseError(f"Unsupported command: {args.command}")


def main() -> int:
    """Run the selected command and present expected failures without a traceback."""
    try:
        _dispatch(_parser().parse_args())
    except ReleaseError as error:
        print(f"release error: {error}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
