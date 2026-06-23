#!/usr/bin/env bash
set -euo pipefail

# Builds a macOS .icns file from soltrace_icon.png.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

INPUT_PNG="$SCRIPT_DIR/soltrace_icon.png"
ICONSET_DIR="$SCRIPT_DIR/SolTrace.iconset"
OUTPUT_ICNS="$SCRIPT_DIR/SolTrace.icns"

if [[ ! -f "$INPUT_PNG" ]]; then
  echo "Error: expected input image not found:"
  echo "  $INPUT_PNG"
  exit 1
fi

if ! command -v sips >/dev/null 2>&1; then
  echo "Error: sips not found. This script must be run on macOS."
  exit 1
fi

if ! command -v iconutil >/dev/null 2>&1; then
  echo "Error: iconutil not found. This script must be run on macOS."
  exit 1
fi

echo "Input:  $INPUT_PNG"
echo "Output: $OUTPUT_ICNS"

rm -rf "$ICONSET_DIR"
mkdir -p "$ICONSET_DIR"

sips -z 16 16       "$INPUT_PNG" --out "$ICONSET_DIR/icon_16x16.png" >/dev/null
sips -z 32 32       "$INPUT_PNG" --out "$ICONSET_DIR/icon_16x16@2x.png" >/dev/null
sips -z 32 32       "$INPUT_PNG" --out "$ICONSET_DIR/icon_32x32.png" >/dev/null
sips -z 64 64       "$INPUT_PNG" --out "$ICONSET_DIR/icon_32x32@2x.png" >/dev/null
sips -z 128 128     "$INPUT_PNG" --out "$ICONSET_DIR/icon_128x128.png" >/dev/null
sips -z 256 256     "$INPUT_PNG" --out "$ICONSET_DIR/icon_128x128@2x.png" >/dev/null
sips -z 256 256     "$INPUT_PNG" --out "$ICONSET_DIR/icon_256x256.png" >/dev/null
sips -z 512 512     "$INPUT_PNG" --out "$ICONSET_DIR/icon_256x256@2x.png" >/dev/null
sips -z 512 512     "$INPUT_PNG" --out "$ICONSET_DIR/icon_512x512.png" >/dev/null
sips -z 1024 1024   "$INPUT_PNG" --out "$ICONSET_DIR/icon_512x512@2x.png" >/dev/null

rm -f "$OUTPUT_ICNS"
iconutil -c icns "$ICONSET_DIR" -o "$OUTPUT_ICNS"

echo "Created:"
echo "  $ICONSET_DIR"
echo "  $OUTPUT_ICNS"

rm -rf "$ICONSET_DIR"