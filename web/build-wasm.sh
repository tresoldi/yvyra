#!/bin/bash
#
# Build yvyra WASM for browser use.
# Requires Emscripten SDK (emsdk) to be available.
#
# Usage:
#   source ~/emsdk/emsdk_env.sh
#   ./web/build-wasm.sh
#

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
BUILD_DIR="$REPO_DIR/build-wasm-web"

echo "Building yvyra WASM for browser..."

# Configure with Emscripten
emcmake cmake -B "$BUILD_DIR" -S "$REPO_DIR" \
    -DENABLE_SSE=OFF -DENABLE_AVX=OFF -DENABLE_FMA=OFF \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_EXE_LINKER_FLAGS="\
        -s MODULARIZE=1 \
        -s EXPORT_NAME=createYvyra \
        -s ALLOW_MEMORY_GROWTH=1 \
        -s EXIT_RUNTIME=0 \
        -s INVOKE_RUN=0 \
        -s EXPORTED_FUNCTIONS=['_main'] \
        -s EXPORTED_RUNTIME_METHODS=['callMain','FS'] \
        -s FORCE_FILESYSTEM=1 \
        -s ENVIRONMENT=web,worker,node \
        -s STACK_SIZE=262144 \
    "

# Build
cmake --build "$BUILD_DIR" -j"$(nproc)"

# Copy artifacts to web/
cp "$BUILD_DIR/src/yvyra.js" "$SCRIPT_DIR/"
cp "$BUILD_DIR/src/yvyra.wasm" "$SCRIPT_DIR/"

echo ""
echo "WASM build complete:"
ls -lh "$SCRIPT_DIR/yvyra.js" "$SCRIPT_DIR/yvyra.wasm"
echo ""
echo "To test: cd web && python3 -m http.server 8080"
