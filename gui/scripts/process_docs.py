#!/usr/bin/env python3

import argparse
import json
import re
import shutil
import subprocess
import sys
from pathlib import Path

MATH_TEXT_SIZE = "18pt"
TYPST_TEMPLATE = (
    "#set page(width: auto, height: auto, margin: 0pt, fill: none)\n"
    "#set text(size: {math_text_size})\n"
    "$ {math_source} $\n"
)


def default_gui_dir():
    return Path(__file__).resolve().parents[1]


def parse_frontmatter(text):
    if not text.lstrip().startswith("---"):
        return {}, text.strip()

    start = text.find("---")
    end = text.find("---", start + 3)
    if end == -1:
        raise RuntimeError("Unclosed frontmatter block.")

    metadata = {}
    frontmatter = text[start + 3 : end].strip()
    body = text[end + 3 :].strip()

    for line in frontmatter.splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        metadata[key.strip().lower()] = value.strip().strip("\"'")

    return metadata, body


def run_typst(compiler, typst_source):
    result = subprocess.run(
        [compiler, "compile", "-", "-", "--format", "svg"],
        input=typst_source,
        capture_output=True,
        text=True,
        check=False,
    )

    if result.returncode != 0:
        details = [
            "Typst failed to compile math markup.",
            f"Exit code: {result.returncode}",
        ]
        if result.stderr:
            details.append(f"stderr:\n{result.stderr.rstrip()}")
        if result.stdout:
            details.append(f"stdout:\n{result.stdout.rstrip()}")
        raise RuntimeError("\n".join(details))

    return result.stdout


def make_math_block(compiler, content, source, line_number):
    math_source = content.strip()
    typst_source = TYPST_TEMPLATE.format(
        math_text_size=MATH_TEXT_SIZE,
        math_source=math_source,
    )

    try:
        svg = run_typst(compiler, typst_source)
    except RuntimeError as error:
        raise RuntimeError(
            f"Failed to process math block in {source}:{line_number}\n{error}"
        ) from error

    return {
        "type": "math",
        "source": math_source,
        "svg": svg,
    }


def split_blocks(compiler, body, source):
    blocks = []
    cursor = 0

    for match in re.finditer(r"\$\$(.*?)\$\$", body, flags=re.DOTALL):
        text = body[cursor : match.start()]
        if text:
            blocks.append({"type": "text", "content": text})

        line_number = body.count("\n", 0, match.start()) + 1
        blocks.append(make_math_block(compiler, match.group(1), source, line_number))
        cursor = match.end()

    text = body[cursor:]
    if text:
        blocks.append({"type": "text", "content": text})

    return blocks


def processed_path(source_root, dest_root, source):
    relative = source.relative_to(source_root)
    if relative.parts and relative.parts[0] == "script":
        return dest_root / relative
    return (dest_root / relative).with_suffix(".json")


def process_markdown(compiler, source_root, dest_root, source):
    destination = processed_path(source_root, dest_root, source)
    destination.parent.mkdir(parents=True, exist_ok=True)

    if source.relative_to(source_root).parts[0:1] == ("script",):
        shutil.copy2(source, destination)
        return

    raw = source.read_text(encoding="utf-8")
    metadata, body = parse_frontmatter(raw)

    processed = {
        "metadata": metadata,
        "blocks": split_blocks(compiler, body, source),
    }
    destination.write_text(
        json.dumps(processed, ensure_ascii=False) + "\n", encoding="utf-8"
    )


def process_docs(source_root, dest_root):
    compiler = shutil.which("typst")
    if compiler is None:
        raise RuntimeError("Unable to find typst. Install typst to process docs.")

    if dest_root.exists():
        shutil.rmtree(dest_root)
    dest_root.mkdir(parents=True, exist_ok=True)

    for source in sorted(source_root.rglob("*.md")):
        process_markdown(compiler, source_root, dest_root, source)


def main():
    gui_dir = default_gui_dir()
    parser = argparse.ArgumentParser(
        description="Process SolTrace GUI documentation sources."
    )
    parser.add_argument(
        "--source",
        type=Path,
        default=gui_dir / "docs_source",
        help="Editable documentation source directory.",
    )
    parser.add_argument(
        "--dest",
        type=Path,
        default=gui_dir / "docs",
        help="Generated documentation resource directory.",
    )
    args = parser.parse_args()

    process_docs(args.source.resolve(), args.dest.resolve())


if __name__ == "__main__":
    try:
        main()
    except RuntimeError as error:
        print(error, file=sys.stderr)
        raise SystemExit(1)
