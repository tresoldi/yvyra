#!/usr/bin/env python3
"""
Cross-validation: compare yvyra-c output against upstream MrBayes 3.2.7a.

Runs both programs on the same nexus file with the same seed and compares
the .p (parameters) and .t (trees) files generation-by-generation.

Usage:
    python3 cross_validate.py --yvyra ./build-novec/src/yvyra \
                              --mrbayes ../mrbayes-upstream/src/mb \
                              --nex testing/crossval_simple.nex
"""

import argparse
import os
import re
import shutil
import subprocess
import sys
import tempfile


def run_program(binary, nex_file, workdir):
    """Run a program on a nexus file, return (stdout+stderr, returncode)."""
    result = subprocess.run(
        [binary, os.path.basename(nex_file)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=300,
    )
    return result.stdout + result.stderr, result.returncode


def read_p_file(directory, suffix=".p"):
    """Read .p file, return (header, rows). Rows are lists of strings."""
    p_files = sorted([f for f in os.listdir(directory) if f.endswith(suffix)])
    if not p_files:
        return None, None
    path = os.path.join(directory, p_files[0])
    with open(path) as f:
        lines = f.readlines()
    data_lines = [l for l in lines if l.strip() and not l.strip().startswith("[")]
    if len(data_lines) < 2:
        return None, None
    header = data_lines[0].strip().split("\t")
    rows = [line.strip().split("\t") for line in data_lines[1:]]
    return header, rows


def read_t_file(directory):
    """Read .t file, return list of (gen, newick) tuples."""
    t_files = sorted([f for f in os.listdir(directory) if f.endswith(".t")])
    if not t_files:
        return None
    path = os.path.join(directory, t_files[0])
    trees = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith("tree "):
                # Format: "tree gen.N = [&U] (newick);"
                m = re.match(r"tree gen\.(\d+)\s*=\s*(\[.*?\]\s*)?(.*)", line)
                if m:
                    gen = int(m.group(1))
                    newick = m.group(3).strip()
                    trees.append((gen, newick))
    return trees


def compare_values(val_a, val_b, tol=1e-10):
    """Compare two string values. Returns (match, diff_info)."""
    try:
        fa = float(val_a)
        fb = float(val_b)
        diff = abs(fa - fb)
        if diff <= tol:
            return True, None
        rel_diff = diff / max(abs(fa), abs(fb), 1e-300)
        return False, f"{fa:.10e} vs {fb:.10e} (diff={diff:.2e}, rel={rel_diff:.2e})"
    except ValueError:
        if val_a == val_b:
            return True, None
        return False, f"'{val_a}' vs '{val_b}'"


def compare_p_files(header_a, rows_a, header_b, rows_b, label_a="yvyra", label_b="mrbayes"):
    """Compare two .p file contents. Returns (n_match, n_differ, details)."""
    details = []

    # Compare headers (skip lnL columns which yvyra puts in .slk)
    common_headers = []
    cols_a = {}
    cols_b = {}
    for i, h in enumerate(header_a):
        if not h.startswith("lnL("):
            cols_a[h] = i
    for i, h in enumerate(header_b):
        if not h.startswith("lnL("):
            cols_b[h] = i
    common_headers = [h for h in cols_a if h in cols_b]

    if set(cols_a.keys()) - set(cols_b.keys()):
        details.append(f"Headers only in {label_a}: {set(cols_a.keys()) - set(cols_b.keys())}")
    if set(cols_b.keys()) - set(cols_a.keys()):
        details.append(f"Headers only in {label_b}: {set(cols_b.keys()) - set(cols_a.keys())}")

    n_match = 0
    n_differ = 0
    n_rows = min(len(rows_a), len(rows_b))
    if len(rows_a) != len(rows_b):
        details.append(f"Row count: {label_a}={len(rows_a)}, {label_b}={len(rows_b)}")

    for row_idx in range(n_rows):
        for h in common_headers:
            ca = cols_a[h]
            cb = cols_b[h]
            if ca < len(rows_a[row_idx]) and cb < len(rows_b[row_idx]):
                match, info = compare_values(rows_a[row_idx][ca], rows_b[row_idx][cb])
                if match:
                    n_match += 1
                else:
                    n_differ += 1
                    if n_differ <= 20:  # limit output
                        gen = rows_a[row_idx][0] if rows_a[row_idx] else "?"
                        details.append(f"  Row {row_idx} (gen {gen}), col '{h}': {info}")

    return n_match, n_differ, details


def compare_t_files(trees_a, trees_b, label_a="yvyra", label_b="mrbayes"):
    """Compare tree files. Returns (n_match, n_differ, details)."""
    details = []
    n_match = 0
    n_differ = 0
    n_trees = min(len(trees_a), len(trees_b))

    if len(trees_a) != len(trees_b):
        details.append(f"Tree count: {label_a}={len(trees_a)}, {label_b}={len(trees_b)}")

    for i in range(n_trees):
        gen_a, nwk_a = trees_a[i]
        gen_b, nwk_b = trees_b[i]
        if gen_a != gen_b:
            details.append(f"  Tree {i}: generation mismatch ({gen_a} vs {gen_b})")
            n_differ += 1
            continue
        if nwk_a == nwk_b:
            n_match += 1
        else:
            n_differ += 1
            if n_differ <= 5:
                details.append(f"  Tree {i} (gen {gen_a}): Newick strings differ")
                # Show first difference position
                for j in range(min(len(nwk_a), len(nwk_b))):
                    if nwk_a[j] != nwk_b[j]:
                        details.append(f"    First diff at pos {j}: ...{nwk_a[max(0,j-10):j+20]}...")
                        details.append(f"    vs                    : ...{nwk_b[max(0,j-10):j+20]}...")
                        break

    return n_match, n_differ, details


def main():
    parser = argparse.ArgumentParser(description="Cross-validate yvyra-c against upstream MrBayes")
    parser.add_argument("--yvyra", required=True, help="Path to yvyra binary")
    parser.add_argument("--mrbayes", required=True, help="Path to upstream mb binary")
    parser.add_argument("--nex", required=True, help="Path to nexus test file")
    parser.add_argument("--keep", action="store_true", help="Keep temp directories")
    args = parser.parse_args()

    yvyra_bin = os.path.abspath(args.yvyra)
    mb_bin = os.path.abspath(args.mrbayes)
    nex_path = os.path.abspath(args.nex)
    nex_name = os.path.basename(nex_path)

    print(f"=== Cross-validation: {nex_name} ===")
    print(f"  yvyra:   {yvyra_bin}")
    print(f"  mrbayes: {mb_bin}")
    print()

    # Create temp directories
    dir_yv = tempfile.mkdtemp(prefix="xval_yvyra_")
    dir_mb = tempfile.mkdtemp(prefix="xval_mb_")

    try:
        # Copy nexus file
        shutil.copy2(nex_path, dir_yv)
        shutil.copy2(nex_path, dir_mb)

        # Run both
        print("Running yvyra-c...")
        out_yv, rc_yv = run_program(yvyra_bin, os.path.join(dir_yv, nex_name), dir_yv)
        print(f"  Exit code: {rc_yv}")

        print("Running MrBayes...")
        out_mb, rc_mb = run_program(mb_bin, os.path.join(dir_mb, nex_name), dir_mb)
        print(f"  Exit code: {rc_mb}")

        if rc_yv != 0:
            print(f"\nERROR: yvyra failed with rc={rc_yv}")
            print(out_yv[-500:])
            return 1
        if rc_mb != 0:
            print(f"\nERROR: MrBayes failed with rc={rc_mb}")
            print(out_mb[-500:])
            return 1

        # Compare .p files
        print("\n--- .p file comparison ---")
        h_yv, r_yv = read_p_file(dir_yv)
        h_mb, r_mb = read_p_file(dir_mb)

        if h_yv is None:
            print("  ERROR: Could not read yvyra .p file")
        elif h_mb is None:
            print("  ERROR: Could not read MrBayes .p file")
        else:
            n_match, n_differ, details = compare_p_files(h_yv, r_yv, h_mb, r_mb)
            print(f"  Matching values: {n_match}")
            print(f"  Differing values: {n_differ}")
            if details:
                for d in details:
                    print(f"  {d}")
            if n_differ == 0:
                print("  RESULT: .p files IDENTICAL")
            else:
                print(f"  RESULT: .p files DIFFER ({n_differ} values)")

        # Compare .t files
        print("\n--- .t file comparison ---")
        t_yv = read_t_file(dir_yv)
        t_mb = read_t_file(dir_mb)

        if t_yv is None:
            print("  ERROR: Could not read yvyra .t file")
        elif t_mb is None:
            print("  ERROR: Could not read MrBayes .t file")
        else:
            n_match, n_differ, details = compare_t_files(t_yv, t_mb)
            print(f"  Matching trees: {n_match}")
            print(f"  Differing trees: {n_differ}")
            if details:
                for d in details:
                    print(f"  {d}")
            if n_differ == 0:
                print("  RESULT: .t files IDENTICAL")
            else:
                print(f"  RESULT: .t files DIFFER ({n_differ} trees)")

        # Extract best likelihoods from output
        print("\n--- Best likelihood comparison ---")
        best_yv = re.findall(r'Likelihood of best state for "cold" chain.*?was ([-\d.]+)', out_yv)
        best_mb = re.findall(r'Likelihood of best state for "cold" chain.*?was ([-\d.]+)', out_mb)
        if best_yv:
            print(f"  yvyra best lnL:   {best_yv[-1]}")
        if best_mb:
            print(f"  MrBayes best lnL: {best_mb[-1]}")
        if best_yv and best_mb:
            diff = abs(float(best_yv[-1]) - float(best_mb[-1]))
            print(f"  Difference: {diff:.6f}")

        print()

    finally:
        if not args.keep:
            shutil.rmtree(dir_yv, ignore_errors=True)
            shutil.rmtree(dir_mb, ignore_errors=True)
        else:
            print(f"Temp dirs kept: {dir_yv}, {dir_mb}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
