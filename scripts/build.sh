#!/bin/bash
# OpenRTX Shader Build Script
# Mirrors the CI/CD workflow for local development
#
# Usage:
#   ./scripts/build.sh              # Build with default settings
#   ./scripts/build.sh --clean      # Clean build directory first
#   ./scripts/build.sh --verbose    # Show verbose output

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Script directory
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"

# Default paths (can be overridden by environment)
DXC_PATH="${DXC_PATH:-/opt/dxc/bin/x64/dxc}"
SHADERC_PATH="${SHADERC_PATH:-/usr/local/bin/shaderc}"
BUILD_DIR="${BUILD_DIR:-$PROJECT_ROOT/build}"
PROJECT_DIR="${PROJECT_DIR:-$PROJECT_ROOT/project}"
VANILLA_DIR="${VANILLA_DIR:-$PROJECT_ROOT/vanilla}"

# Parse arguments
CLEAN=false
VERBOSE=false
for arg in "$@"; do
    case $arg in
        --clean)
            CLEAN=true
            shift
            ;;
        --verbose)
            VERBOSE=true
            shift
            ;;
        --help|-h)
            echo "OpenRTX Shader Build Script"
            echo ""
            echo "Usage: $0 [options]"
            echo ""
            echo "Options:"
            echo "  --clean     Clean build directory before building"
            echo "  --verbose   Show verbose output"
            echo "  --help      Show this help message"
            echo ""
            echo "Environment variables:"
            echo "  DXC_PATH     Path to DXC compiler (default: /opt/dxc/bin/x64/dxc)"
            echo "  SHADERC_PATH Path to shaderc (default: /usr/local/bin/shaderc)"
            echo "  BUILD_DIR    Output directory (default: ./build)"
            exit 0
            ;;
    esac
done

echo -e "${GREEN}=== OpenRTX Shader Build ===${NC}"
echo ""

# Check for required tools
echo "Checking build tools..."

if ! command -v lazurite &> /dev/null; then
    echo -e "${RED}Error: lazurite not found. Install with: pip install lazurite${NC}"
    exit 1
fi
echo -e "  ${GREEN}✓${NC} Lazurite: $(pip show lazurite 2>/dev/null | grep Version | cut -d' ' -f2)"

if [ -f "$DXC_PATH" ] || command -v dxc &> /dev/null; then
    echo -e "  ${GREEN}✓${NC} DXC: $DXC_PATH"
else
    echo -e "${YELLOW}  ⚠ DXC not found at $DXC_PATH (may still work if in PATH)${NC}"
fi

if [ -f "$SHADERC_PATH" ] || command -v shaderc &> /dev/null; then
    echo -e "  ${GREEN}✓${NC} Shaderc: $SHADERC_PATH"
else
    echo -e "${YELLOW}  ⚠ Shaderc not found at $SHADERC_PATH (may still work if in PATH)${NC}"
fi

# Check for vanilla materials
if [ ! -d "$VANILLA_DIR" ] || [ -z "$(ls -A "$VANILLA_DIR" 2>/dev/null)" ]; then
    echo -e "${RED}Error: Vanilla materials not found in $VANILLA_DIR${NC}"
    echo "Please download vanilla .material.bin files or set VANILLA_DIR"
    exit 1
fi
echo -e "  ${GREEN}✓${NC} Vanilla materials: $(ls "$VANILLA_DIR"/*.material.bin 2>/dev/null | wc -l) files"

echo ""

# Run template generation if script exists
if [ -f "$PROJECT_ROOT/src/script.py" ]; then
    echo "Running template generation..."
    cd "$PROJECT_ROOT/src"
    python script.py
    cd "$PROJECT_ROOT"
    echo -e "${GREEN}✓${NC} Template generation complete"
    echo ""
fi

# Clean if requested
if [ "$CLEAN" = true ]; then
    echo "Cleaning build directory..."
    rm -rf "$BUILD_DIR"
    echo -e "${GREEN}✓${NC} Cleaned"
    echo ""
fi

# Create build directory
mkdir -p "$BUILD_DIR"

# Build shaders
echo "Building shaders..."
echo "  Project: $PROJECT_DIR"
echo "  Output:  $BUILD_DIR"
echo ""

BUILD_CMD="lazurite build \"$PROJECT_DIR\" -o \"$BUILD_DIR\""

# Add compiler paths if they exist
if [ -f "$SHADERC_PATH" ]; then
    BUILD_CMD="$BUILD_CMD --shaderc \"$SHADERC_PATH\""
fi
if [ -f "$DXC_PATH" ]; then
    BUILD_CMD="$BUILD_CMD --dxc \"$DXC_PATH\""
fi

if [ "$VERBOSE" = true ]; then
    echo "Command: $BUILD_CMD"
    echo ""
fi

# Execute build
eval $BUILD_CMD

echo ""
echo -e "${GREEN}=== Build Complete ===${NC}"
echo ""
echo "Output files:"
ls -la "$BUILD_DIR"/*.material.bin 2>/dev/null || echo "  (no .material.bin files found)"
