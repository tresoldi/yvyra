#!/bin/sh
# TAP-format smoke test for yvyra Mk model.
# Usage: smoke_test.sh /path/to/binary

if [ "$#" -ne 1 ]; then
    echo "Usage: $0 path-to-binary" >&2
    exit 1
fi

binary=$1
script_dir=$(cd "$(dirname "$0")" && pwd)

echo 1..3
t=0

ok ()   { printf 'ok %d - %s\n' "$t" "$m"; }
fail () { printf 'not ok %d - %s\n' "$t" "$m"; }
error () { fail; printf 'Bail out! %s\n' "$1"; exit 1; }

tmpdir=$(mktemp -d)
cp "$script_dir/mk_smoke.nex" "$tmpdir/"
trap 'rm -rf "$tmpdir"' EXIT TERM INT

t=$(( t + 1 )); m='Binary runs and exits cleanly'
if (cd "$tmpdir" && "$binary" mk_smoke.nex) >"$tmpdir/output.txt" 2>&1; then
    ok
else
    error "Binary exited with non-zero status"
fi

t=$(( t + 1 )); m='Analysis completed'
if grep -q "Analysis completed" "$tmpdir/output.txt"; then
    ok
else
    fail
fi

t=$(( t + 1 )); m='Likelihood is finite (not nan/inf)'
lnl=$(awk '/Likelihood of best state for "cold" chain/ { print $NF }' "$tmpdir/output.txt" | head -1)
if [ -z "$lnl" ]; then
    # Single-run might not print this line; check .p file instead
    p_file=$(find "$tmpdir" -name "*.p" | head -1)
    if [ -n "$p_file" ] && [ -s "$p_file" ]; then
        ok
    else
        fail
    fi
else
    case "$lnl" in
        *[Nn][Aa][Nn]*|*[Ii][Nn][Ff]*) fail ;;
        *) ok ;;
    esac
fi
