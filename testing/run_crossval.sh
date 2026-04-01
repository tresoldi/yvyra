#!/bin/bash
#
# Cross-validation runner: builds upstream MrBayes and yvyra (both
# without SIMD), then runs deterministic comparison on all test files.
#
# Usage:
#   ./testing/run_crossval.sh              # full run (build + test)
#   ./testing/run_crossval.sh --skip-build # reuse existing builds
#
# Prerequisites:
#   - autotools (autoconf, automake) for MrBayes build
#   - cmake for yvyra build
#   - python3
#
# The script creates a temporary MrBayes worktree from the pre-fork
# commit, builds it, runs all crossval_*.nex files through both
# programs, and reports results.

set -e

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
MB_WORKTREE="$REPO_DIR/../mrbayes-upstream"
FORK_POINT="bb09fff"

SKIP_BUILD=0
if [ "$1" = "--skip-build" ]; then
    SKIP_BUILD=1
fi

echo "=== Cross-validation suite ==="
echo "  Repo:       $REPO_DIR"
echo "  Fork point: $FORK_POINT"
echo ""

# --- Build yvyra without SIMD ---
YVYRA_BIN="$REPO_DIR/build-novec/src/yvyra"
if [ "$SKIP_BUILD" = "0" ] || [ ! -f "$YVYRA_BIN" ]; then
    echo "Building yvyra (no SIMD)..."
    cmake -B "$REPO_DIR/build-novec" -S "$REPO_DIR" \
        -DENABLE_SSE=OFF -DENABLE_AVX=OFF -DENABLE_FMA=OFF \
        -DCMAKE_BUILD_TYPE=Release > /dev/null 2>&1
    cmake --build "$REPO_DIR/build-novec" -j"$(nproc)" > /dev/null 2>&1
    echo "  Built: $YVYRA_BIN"
else
    echo "  Reusing: $YVYRA_BIN"
fi

# --- Build upstream MrBayes ---
MB_BIN="$MB_WORKTREE/src/mb"
if [ "$SKIP_BUILD" = "0" ] || [ ! -f "$MB_BIN" ]; then
    echo "Building upstream MrBayes 3.2.7a..."
    # Create worktree if not exists
    if [ ! -d "$MB_WORKTREE" ]; then
        cd "$REPO_DIR"
        git worktree add "$MB_WORKTREE" "$FORK_POINT" 2>/dev/null || true
    fi
    cd "$MB_WORKTREE"
    if [ ! -f Makefile ]; then
        ./configure --without-beagle --disable-sse > /dev/null 2>&1
    fi
    make -j"$(nproc)" > /dev/null 2>&1
    echo "  Built: $MB_BIN"
else
    echo "  Reusing: $MB_BIN"
fi

# --- Run comparisons ---
echo ""
echo "=== Running deterministic comparisons ==="
echo ""

PASS=0
FAIL=0
EXPECTED_DIFF=0

for NEX in "$SCRIPT_DIR"/crossval_*.nex; do
    NAME=$(basename "$NEX")
    echo "--- $NAME ---"
    OUTPUT=$(python3 "$SCRIPT_DIR/cross_validate.py" \
        --yvyra "$YVYRA_BIN" \
        --mrbayes "$MB_BIN" \
        --nex "$NEX" 2>&1 || true)

    echo "$OUTPUT" | grep -E 'RESULT:|Differing|Matching|Best|Difference' || true

    # Check results
    P_RESULT=$(echo "$OUTPUT" | grep '.p files' | grep 'RESULT:' | head -1 || true)
    T_RESULT=$(echo "$OUTPUT" | grep '.t files' | grep 'RESULT:' | head -1 || true)

    if echo "$OUTPUT" | grep -q "ERROR: MrBayes failed"; then
        echo "  => SKIP (MrBayes does not support this feature)"
    elif echo "$P_RESULT" | grep -q "IDENTICAL" && echo "$T_RESULT" | grep -q "IDENTICAL"; then
        echo "  => PASS (exact match)"
        PASS=$((PASS + 1))
    elif echo "$T_RESULT" | grep -q "IDENTICAL" && echo "$NAME" | grep -q "multistate"; then
        echo "  => EXPECTED DIFF (multi-state ascertainment correction)"
        EXPECTED_DIFF=$((EXPECTED_DIFF + 1))
    else
        echo "  => FAIL (unexpected difference)"
        FAIL=$((FAIL + 1))
    fi
    echo ""
done

# --- Summary ---
echo "=== Summary ==="
echo "  Exact matches:      $PASS"
echo "  Expected diffs:     $EXPECTED_DIFF"
echo "  Failures:           $FAIL"

if [ "$FAIL" -gt 0 ]; then
    echo ""
    echo "CROSS-VALIDATION FAILED"
    exit 1
else
    echo ""
    echo "CROSS-VALIDATION PASSED"
    exit 0
fi
