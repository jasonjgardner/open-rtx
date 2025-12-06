#!/bin/bash
# OpenRTX Shader Build Script
# Builds the shader using Docker with DXC compiler

set -e

VANILLA_BASE_URL="https://static.bedrock.graphics.t3.storage.dev/vanilla/v23"
BUILD_DIR="./build"
VANILLA_DIR="./vanilla"

echo "==================================="
echo "OpenRTX Shader Build Script"
echo "==================================="

# Check for Docker
if ! command -v docker &> /dev/null; then
    echo "Error: Docker is required but not installed."
    echo "Please install Docker and try again."
    exit 1
fi

# Create directories
mkdir -p "$BUILD_DIR"
mkdir -p "$VANILLA_DIR"

# Download vanilla materials if not present
echo ""
echo "Checking for vanilla material.bin files..."

download_if_missing() {
    local file=$1
    local url="${VANILLA_BASE_URL}/${file}"
    local path="${VANILLA_DIR}/${file}"

    if [ ! -f "$path" ]; then
        echo "Downloading ${file}..."
        curl -sL "$url" -o "$path"
    else
        echo "Found ${file}"
    fi
}

download_if_missing "RTXStub.material.bin"
download_if_missing "RTXPostFX.material.bin"
download_if_missing "RTXPostFX.Tonemapping.material.bin"
download_if_missing "RTXPostFX.Bloom.material.bin"

echo ""
echo "Building shader with Docker..."
echo ""

# Build using Docker
docker run --rm \
    -v "$(pwd):/shader" \
    -w /shader \
    jasongardner/dxc:latest \
    sh -c "
        apt-get update -qq && apt-get install -y -qq python3 python3-pip > /dev/null 2>&1 &&
        pip3 install -q lazurite-mc > /dev/null 2>&1 &&
        echo 'Compiling shaders...' &&
        lazurite build project/ -o build/ &&
        echo '' &&
        echo 'Build complete!'
    "

echo ""
echo "==================================="
echo "Build Output:"
echo "==================================="
ls -la "$BUILD_DIR"/*.material.bin 2>/dev/null || echo "No .material.bin files found in build/"

echo ""
echo "Done! Shader files are in: $BUILD_DIR/"
