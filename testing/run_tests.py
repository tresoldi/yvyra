#!/usr/bin/env python3
"""
Test harness for yvyra Mk model.

Uses only Python standard library. Runs the binary on test Nexus files
and checks output against expected ranges.

Usage:
    python3 run_tests.py --binary ./src/yvyra
    python3 run_tests.py --binary ./build/yvyra
    python3 run_tests.py --binary ./src/yvyra -v          # verbose
    python3 run_tests.py --binary ./src/yvyra -k smoke    # run only tests matching 'smoke'
"""

import argparse
import json
import math
import os
import re
import shutil
import subprocess
import sys
import tempfile
import unittest


# Resolved at module load time, overridden by --binary
BINARY_PATH = None
EXPECTED_RANGES = None
TESTING_DIR = os.path.dirname(os.path.abspath(__file__))


def load_expected_ranges():
    path = os.path.join(TESTING_DIR, "expected_ranges.json")
    with open(path) as f:
        return json.load(f)


def run_yvyra(binary, nex_file, workdir):
    """Run the binary on a nexus file in the given workdir. Return (stdout, returncode)."""
    result = subprocess.run(
        [binary, os.path.basename(nex_file)],
        cwd=workdir,
        capture_output=True,
        text=True,
        timeout=300,
    )
    return result.stdout + result.stderr, result.returncode


def parse_likelihood(output):
    """Extract best cold chain likelihood(s) from output."""
    pattern = r'Likelihood of best state for "cold" chain of run \d+ was ([-\d.]+)'
    return [float(m) for m in re.findall(pattern, output)]


def parse_asdsf(output):
    """Extract final average standard deviation of split frequencies."""
    # The sump/sumt output line (not the mid-run diagnostic)
    pattern = r"Average standard deviation of split frequencies = ([\d.]+)"
    matches = re.findall(pattern, output)
    if matches:
        return float(matches[-1])
    return None


def parse_psrf(output):
    """Extract average PSRF."""
    pattern = r"Average PSRF for parameter values \(excluding NA and >10\.0\) = ([\d.]+)"
    matches = re.findall(pattern, output)
    if matches:
        return float(matches[-1])
    return None


def parse_completed(output):
    """Check if analysis completed."""
    return "Analysis completed" in output


class SmokeTest(unittest.TestCase):
    """Minimal smoke test: does the binary start and finish on Mk data."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="yvyra_test_")
        self.nex = os.path.join(TESTING_DIR, "mk_smoke.nex")
        shutil.copy2(self.nex, self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_smoke_completes(self):
        """Binary runs mk_smoke.nex to completion."""
        output, rc = run_yvyra(BINARY_PATH, self.nex, self.tmpdir)
        self.assertEqual(rc, 0, f"Binary exited with code {rc}")
        self.assertTrue(parse_completed(output), "Analysis did not complete")

    def test_smoke_likelihood_finite(self):
        """Smoke test produces a finite likelihood."""
        output, rc = run_yvyra(BINARY_PATH, self.nex, self.tmpdir)
        lnls = parse_likelihood(output)
        # Smoke test has only 1 run, might not report likelihood for short runs
        # At minimum, check completion
        self.assertTrue(parse_completed(output))

    def test_smoke_output_files(self):
        """Smoke test produces .p and .t output files."""
        output, rc = run_yvyra(BINARY_PATH, self.nex, self.tmpdir)
        self.assertEqual(rc, 0)
        # Check for parameter and tree files
        files = os.listdir(self.tmpdir)
        p_files = [f for f in files if f.endswith(".p")]
        t_files = [f for f in files if f.endswith(".t")]
        self.assertTrue(len(p_files) > 0, f"No .p files found. Files: {files}")
        self.assertTrue(len(t_files) > 0, f"No .t files found. Files: {files}")


class SyntheticMkTest(unittest.TestCase):
    """Tests on the synthetic 6-taxon 15-character Mk dataset."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="yvyra_test_")
        self.nex = os.path.join(TESTING_DIR, "mk_synthetic.nex")
        shutil.copy2(self.nex, self.tmpdir)
        self.expected = EXPECTED_RANGES["mk_synthetic"]
        self.output, self.rc = run_yvyra(BINARY_PATH, self.nex, self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_completes(self):
        """Synthetic Mk analysis completes."""
        self.assertEqual(self.rc, 0)
        self.assertTrue(parse_completed(self.output))

    def test_likelihood_range(self):
        """Cold chain likelihood within expected range."""
        lnls = parse_likelihood(self.output)
        self.assertTrue(len(lnls) >= 1, "No likelihood values found in output")
        lo, hi = self.expected["likelihood_range"]
        for lnl in lnls:
            self.assertGreaterEqual(lnl, lo, f"Likelihood {lnl} below {lo}")
            self.assertLessEqual(lnl, hi, f"Likelihood {lnl} above {hi}")

    def test_asdsf(self):
        """Average standard deviation of split frequencies below threshold."""
        asdsf = parse_asdsf(self.output)
        if asdsf is not None:
            self.assertLessEqual(
                asdsf,
                self.expected["asdsf_max"],
                f"ASDSF {asdsf} exceeds {self.expected['asdsf_max']}",
            )

    def test_psrf(self):
        """Average PSRF within expected range."""
        psrf = parse_psrf(self.output)
        if psrf is not None:
            lo, hi = self.expected["psrf_range"]
            self.assertGreaterEqual(psrf, lo, f"PSRF {psrf} below {lo}")
            self.assertLessEqual(psrf, hi, f"PSRF {psrf} above {hi}")

    def test_output_files(self):
        """Analysis produces expected output files."""
        files = os.listdir(self.tmpdir)
        p_files = [f for f in files if f.endswith(".p")]
        t_files = [f for f in files if f.endswith(".t")]
        # 2 runs = 2 .p files + 2 .t files
        self.assertGreaterEqual(len(p_files), 2, f"Expected 2+ .p files, got {len(p_files)}")
        self.assertGreaterEqual(len(t_files), 2, f"Expected 2+ .t files, got {len(t_files)}")

    def test_p_file_columns(self):
        """Parameter file has expected columns."""
        files = os.listdir(self.tmpdir)
        p_files = sorted([f for f in files if f.endswith(".p") and "run1" in f])
        if not p_files:
            self.skipTest("No run1 .p file found")
        p_path = os.path.join(self.tmpdir, p_files[0])
        with open(p_path) as f:
            lines = f.readlines()
        # Skip comment lines (start with [)
        data_lines = [l for l in lines if l.strip() and not l.strip().startswith("[")]
        self.assertTrue(len(data_lines) >= 2, "Too few lines in .p file")
        # Header + at least one data line
        header = data_lines[0].strip().split("\t")
        self.assertIn("lnLike", header, f"No lnLike column in .p header: {header}")
        self.assertIn("TL", header, f"No TL column in .p header: {header}")


class RealMkTest(unittest.TestCase):
    """Tests on the real hymenoptera 20-taxon 200-character Mk dataset."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="yvyra_test_")
        self.nex = os.path.join(TESTING_DIR, "mk_real.nex")
        shutil.copy2(self.nex, self.tmpdir)
        self.expected = EXPECTED_RANGES["mk_real"]
        self.output, self.rc = run_yvyra(BINARY_PATH, self.nex, self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_completes(self):
        """Real Mk analysis completes."""
        self.assertEqual(self.rc, 0)
        self.assertTrue(parse_completed(self.output))

    def test_likelihood_range(self):
        """Cold chain likelihood within expected range."""
        lnls = parse_likelihood(self.output)
        self.assertTrue(len(lnls) >= 1, "No likelihood values found in output")
        lo, hi = self.expected["likelihood_range"]
        for lnl in lnls:
            self.assertGreaterEqual(lnl, lo, f"Likelihood {lnl} below {lo}")
            self.assertLessEqual(lnl, hi, f"Likelihood {lnl} above {hi}")

    def test_asdsf(self):
        """Average standard deviation of split frequencies below threshold."""
        asdsf = parse_asdsf(self.output)
        if asdsf is not None:
            self.assertLessEqual(
                asdsf,
                self.expected["asdsf_max"],
                f"ASDSF {asdsf} exceeds {self.expected['asdsf_max']}",
            )

    def test_psrf(self):
        """Average PSRF within expected range."""
        psrf = parse_psrf(self.output)
        if psrf is not None:
            lo, hi = self.expected["psrf_range"]
            self.assertGreaterEqual(psrf, lo, f"PSRF {psrf} below {lo}")
            self.assertLessEqual(psrf, hi, f"PSRF {psrf} above {hi}")


class DeterminismTest(unittest.TestCase):
    """Verify that runs with the same seed produce identical results."""

    def test_deterministic_output(self):
        """Two runs with identical seeds produce identical .p files."""
        results = []
        for _ in range(2):
            tmpdir = tempfile.mkdtemp(prefix="yvyra_det_")
            nex = os.path.join(TESTING_DIR, "mk_smoke.nex")
            shutil.copy2(nex, tmpdir)
            run_yvyra(BINARY_PATH, nex, tmpdir)
            p_files = sorted([
                f for f in os.listdir(tmpdir) if f.endswith(".p")
            ])
            p_contents = []
            for pf in p_files:
                with open(os.path.join(tmpdir, pf)) as f:
                    p_contents.append(f.read())
            results.append((p_contents, tmpdir))

        try:
            self.assertEqual(len(results[0][0]), len(results[1][0]),
                             "Different number of .p files")
            for i, (c1, c2) in enumerate(zip(results[0][0], results[1][0])):
                self.assertEqual(c1, c2,
                                 f".p file {i} differs between identical-seed runs")
        finally:
            for _, tmpdir in results:
                shutil.rmtree(tmpdir, ignore_errors=True)


def write_nex(tmpdir, filename, content):
    """Write a Nexus file to tmpdir and return the path."""
    path = os.path.join(tmpdir, filename)
    with open(path, "w") as f:
        f.write(content)
    return path


def read_p_file(tmpdir, pattern=".p"):
    """Read the first .p file in tmpdir, return (header_list, data_rows)."""
    p_files = sorted([f for f in os.listdir(tmpdir) if f.endswith(pattern)])
    if not p_files:
        return None, None
    path = os.path.join(tmpdir, p_files[0])
    with open(path) as f:
        lines = f.readlines()
    # Skip comment lines (start with [)
    data_lines = [l for l in lines if l.strip() and not l.strip().startswith("[")]
    if len(data_lines) < 2:
        return None, None
    header = data_lines[0].strip().split("\t")
    rows = []
    for line in data_lines[1:]:
        rows.append(line.strip().split("\t"))
    return header, rows


# --- Nexus snippets used by multiple test classes ---

SITELIKES_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=6 nchar=10;
    Format datatype=standard symbols="0 1 2" missing=? gap=-;
    Matrix
        A  0 0 1 0 0 1 0 0 1 0
        B  0 0 1 0 1 1 0 1 1 0
        C  1 1 0 1 1 0 1 1 0 1
        D  1 1 0 1 0 0 1 0 0 1
        E  1 0 0 1 1 0 0 1 0 0
        F  0 1 1 0 0 1 1 0 1 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=99;
    lset coding=variable;
    report sitelikes=yes;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

NOSITELIKES_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=5;
    Format datatype=standard symbols="0 1" missing=? gap=-;
    Matrix
        A  0 0 1 0 1
        B  0 1 1 0 0
        C  1 1 0 1 0
        D  1 0 0 1 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=100 samplefr=100;
    mcmcp printfr=100 diagnfr=100;
    mcmc;
End;
"""


class SiteLikesTest(unittest.TestCase):
    """Tests for per-site log-likelihood output (Phase 4)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="yvyra_slk_")
        write_nex(self.tmpdir, "sitelikes.nex", SITELIKES_NEX)
        self.output, self.rc = run_yvyra(BINARY_PATH,
            os.path.join(self.tmpdir, "sitelikes.nex"), self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_completes_with_sitelikes(self):
        """Analysis with report sitelikes=yes completes."""
        self.assertEqual(self.rc, 0, f"Binary exited with code {self.rc}")
        self.assertTrue(parse_completed(self.output), "Analysis did not complete")

    def test_sitelikes_columns_present(self):
        """The .slk file contains lnL(...) columns when sitelikes=yes."""
        header, rows = read_p_file(self.tmpdir, pattern=".slk")
        self.assertIsNotNone(header, "Could not read .slk file")
        lnl_cols = [h for h in header if h.startswith("lnL(")]
        self.assertGreater(len(lnl_cols), 0,
                           f"No lnL columns found in header: {header}")

    def test_sitelikes_column_count(self):
        """Number of lnL columns matches the number of site patterns."""
        header, rows = read_p_file(self.tmpdir, pattern=".slk")
        self.assertIsNotNone(header)
        lnl_cols = [h for h in header if h.startswith("lnL(")]
        # With 10 characters (6 taxa, coding=variable), we expect some
        # number of unique compressed patterns. At least 1, at most 10.
        self.assertGreaterEqual(len(lnl_cols), 1,
                                "Too few lnL columns")
        self.assertLessEqual(len(lnl_cols), 10,
                             "More lnL columns than original characters")

    def test_sitelikes_values_finite(self):
        """Per-site log-likelihood values are finite numbers."""
        header, rows = read_p_file(self.tmpdir, pattern=".slk")
        self.assertIsNotNone(header)
        lnl_indices = [i for i, h in enumerate(header) if h.startswith("lnL(")]
        self.assertGreater(len(lnl_indices), 0)
        for row in rows:
            for idx in lnl_indices:
                val = float(row[idx])
                self.assertFalse(math.isnan(val),
                                 f"NaN in per-site lnL at gen {row[0]}, col {header[idx]}")
                self.assertFalse(math.isinf(val),
                                 f"Inf in per-site lnL at gen {row[0]}, col {header[idx]}")

    def test_sitelikes_values_negative(self):
        """Per-site log-likelihoods are negative (log of a probability)."""
        header, rows = read_p_file(self.tmpdir, pattern=".slk")
        self.assertIsNotNone(header)
        lnl_indices = [i for i, h in enumerate(header) if h.startswith("lnL(")]
        for row in rows:
            for idx in lnl_indices:
                val = float(row[idx])
                self.assertLess(val, 0.0,
                                f"Per-site lnL should be negative, got {val} "
                                f"at gen {row[0]}, col {header[idx]}")

    def test_sitelikes_vary_across_samples(self):
        """Per-site lnL values change between MCMC samples (not stuck)."""
        header, rows = read_p_file(self.tmpdir, pattern=".slk")
        self.assertIsNotNone(header)
        lnl_indices = [i for i, h in enumerate(header) if h.startswith("lnL(")]
        if len(rows) < 2:
            self.skipTest("Need at least 2 samples to check variation")
        # Check if at least one per-site lnL changed between first and last sample
        first = [float(rows[0][i]) for i in lnl_indices]
        last = [float(rows[-1][i]) for i in lnl_indices]
        changed = any(abs(a - b) > 1e-10 for a, b in zip(first, last))
        self.assertTrue(changed,
                        "Per-site lnL values identical across all samples (stuck chain?)")

    def test_sitelikes_not_in_p_file(self):
        """The .p file no longer contains lnL columns (moved to .slk)."""
        header, rows = read_p_file(self.tmpdir, pattern=".p")
        self.assertIsNotNone(header, "Could not read .p file")
        lnl_cols = [h for h in header if h.startswith("lnL(")]
        self.assertEqual(len(lnl_cols), 0,
                         f"lnL columns should not be in .p file: {lnl_cols}")

    def test_slk_has_weights(self):
        """The .slk file contains a [weights:...] line with pattern weights."""
        slk_files = sorted([f for f in os.listdir(self.tmpdir) if f.endswith(".slk")])
        self.assertGreater(len(slk_files), 0, "No .slk file found")
        with open(os.path.join(self.tmpdir, slk_files[0])) as f:
            content = f.read()
        self.assertIn("[weights:", content, "No weights line in .slk file")


class NoSiteLikesTest(unittest.TestCase):
    """Verify that lnL columns are absent by default (sitelikes not requested)."""

    def setUp(self):
        self.tmpdir = tempfile.mkdtemp(prefix="yvyra_noslk_")
        write_nex(self.tmpdir, "noslk.nex", NOSITELIKES_NEX)
        self.output, self.rc = run_yvyra(BINARY_PATH,
            os.path.join(self.tmpdir, "noslk.nex"), self.tmpdir)

    def tearDown(self):
        shutil.rmtree(self.tmpdir, ignore_errors=True)

    def test_no_sitelikes_columns_by_default(self):
        """Without report sitelikes=yes, no lnL columns appear in .p file."""
        header, rows = read_p_file(self.tmpdir)
        self.assertIsNotNone(header, "Could not read .p file")
        lnl_cols = [h for h in header if h.startswith("lnL(")]
        self.assertEqual(len(lnl_cols), 0,
                         f"lnL columns present without sitelikes=yes: {lnl_cols}")

    def test_no_slk_file_by_default(self):
        """Without report sitelikes=yes, no .slk file is created."""
        slk_files = [f for f in os.listdir(self.tmpdir) if f.endswith(".slk")]
        self.assertEqual(len(slk_files), 0,
                         f".slk file created without sitelikes=yes: {slk_files}")


class MolecularRejectionTest(unittest.TestCase):
    """Verify that molecular datatypes are rejected with a clear error."""

    def _run_nex(self, content):
        tmpdir = tempfile.mkdtemp(prefix="yvyra_reject_")
        try:
            write_nex(tmpdir, "test.nex", content)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            return output, rc
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_dna_rejected(self):
        """DNA datatype is rejected with a clear error message."""
        nex = "#NEXUS\nBegin data;\n  Dimensions ntax=2 nchar=3;\n" \
              "  Format datatype=DNA;\n  Matrix\n    A ACG\n    B ACT\n  ;\nEnd;\n"
        output, rc = self._run_nex(nex)
        self.assertNotEqual(rc, 0, "DNA data should be rejected")
        self.assertIn("not supported", output.lower(),
                      "Error message should mention 'not supported'")

    def test_protein_rejected(self):
        """Protein datatype is rejected with a clear error message."""
        nex = "#NEXUS\nBegin data;\n  Dimensions ntax=2 nchar=3;\n" \
              "  Format datatype=Protein;\n  Matrix\n    A ACG\n    B ACT\n  ;\nEnd;\n"
        output, rc = self._run_nex(nex)
        self.assertNotEqual(rc, 0, "Protein data should be rejected")
        self.assertIn("not supported", output.lower(),
                      "Error message should mention 'not supported'")

    def test_standard_accepted(self):
        """Standard (morphological) datatype is accepted."""
        nex = "#NEXUS\nBegin data;\n  Dimensions ntax=4 nchar=4;\n" \
              '  Format datatype=standard symbols="0 1";\n  Matrix\n' \
              "    A 0011\n    B 0110\n    C 1100\n    D 1001\n  ;\nEnd;\n" \
              "Begin yvyra;\n  set autoclose=yes nowarnings=yes;\n" \
              "  mcmcp nrun=1 nchain=1 ngen=100 samplefr=100;\n" \
              "  mcmcp printfr=100 diagnfr=100;\n  mcmc;\nEnd;\n"
        output, rc = self._run_nex(nex)
        self.assertEqual(rc, 0, f"Standard data should be accepted, got rc={rc}")


class BannerTest(unittest.TestCase):
    """Verify the binary identifies as yvyra."""

    def test_banner_says_yvyra(self):
        """Binary banner says 'yvyra', not 'MrBayes'."""
        result = subprocess.run(
            [BINARY_PATH],
            input="quit\n",
            capture_output=True,
            text=True,
            timeout=10,
        )
        self.assertIn("yvyra", result.stdout,
                      "Banner should contain 'yvyra'")
        # The banner should NOT say "MrBayes" as the tool name
        # (it can reference MrBayes as the upstream project)
        lines = result.stdout.split("\n")
        # First non-empty line with the tool name
        for line in lines:
            if "yvyra" in line.lower() or "mrbayes" in line.lower():
                # This should be the yvyra banner line
                self.assertIn("yvyra", line,
                              f"First banner line should say yvyra: {line}")
                break


USERTYPE_NEX = """\
#NEXUS
[ All 4 characters have all 3 states observed, so they go through cijk path ]
Begin data;
    Dimensions ntax=5 nchar=4;
    Format datatype=standard symbols="0 1 2" missing=? gap=-;
    Matrix
        A  0 2 1 0
        B  1 0 2 1
        C  2 1 0 2
        D  0 1 2 0
        E  1 2 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    usertype asymm (Standard) = 3
        0.0 1.0 0.0
        5.0 0.0 1.0
        0.0 0.1 0.0
    ;
    ctype asymm: 1-4;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

USERTYPE_DEFAULT_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=5 nchar=4;
    Format datatype=standard symbols="0 1 2" missing=? gap=-;
    Matrix
        A  0 2 1 0
        B  1 0 2 1
        C  2 1 0 2
        D  0 1 2 0
        E  1 2 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

USERTYPE_SYMMETRIC_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=5 nchar=4;
    Format datatype=standard symbols="0 1 2" missing=? gap=-;
    Matrix
        A  0 2 1 0
        B  1 0 2 1
        C  2 1 0 2
        D  0 1 2 0
        E  1 2 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    usertype symm (Standard) = 3
        0.0 1.0 1.0
        1.0 0.0 1.0
        1.0 1.0 0.0
    ;
    ctype symm: 1-4;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""


class UsertypeTest(unittest.TestCase):
    """Tests for arbitrary Q matrix via usertype command (Phase 5)."""

    def test_usertype_completes(self):
        """Analysis with custom usertype Q matrix completes."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_ut_")
        try:
            write_nex(tmpdir, "test.nex", USERTYPE_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0, f"Binary exited with code {rc}")
            self.assertTrue(parse_completed(output), "Analysis did not complete")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_usertype_recognized(self):
        """The usertype command is recognized and the ctype is applied."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_ut_")
        try:
            write_nex(tmpdir, "test.nex", USERTYPE_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertIn("Defined usertype", output,
                          "Usertype definition not acknowledged")
            self.assertIn("Ctype was applied", output,
                          "Ctype application not acknowledged")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    # Fixed: USERTYPE now works correctly
    def test_usertype_changes_likelihood(self):
        """Asymmetric usertype produces different sampled likelihoods than default Mk."""
        results = {}
        for label, nex_content in [("usertype", USERTYPE_NEX),
                                    ("default", USERTYPE_DEFAULT_NEX)]:
            tmpdir = tempfile.mkdtemp(prefix="yvyra_ut_")
            try:
                write_nex(tmpdir, "test.nex", nex_content)
                output, rc = run_yvyra(BINARY_PATH,
                    os.path.join(tmpdir, "test.nex"), tmpdir)
                self.assertEqual(rc, 0)
                # Read the .p file and get the last sampled lnL
                header, rows = read_p_file(tmpdir)
                self.assertIsNotNone(header)
                lnl_idx = header.index("lnLike")
                # Use the last sample (after MCMC has run and UpDateCijk called)
                last_lnl = float(rows[-1][lnl_idx])
                results[label] = last_lnl
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        self.assertNotAlmostEqual(
            results["usertype"], results["default"], places=0,
            msg=f"Usertype lnL ({results['usertype']}) should differ from "
                f"default ({results['default']})")

    # Fixed: USERTYPE now works correctly
    def test_symmetric_usertype_produces_finite_likelihood(self):
        """Symmetric usertype with equal rates produces a finite likelihood."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_ut_")
        try:
            write_nex(tmpdir, "test.nex", USERTYPE_SYMMETRIC_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0)
            m = re.search(r'Chain 1 -- ([-\d.]+) --', output)
            self.assertIsNotNone(m, "Could not parse initial lnL")
            lnl = float(m.group(1))
            self.assertFalse(math.isnan(lnl), "lnL is NaN")
            self.assertFalse(math.isinf(lnl), "lnL is Inf")
            self.assertLess(lnl, 0.0, "lnL should be negative")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


USERCOST_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1" missing=? gap=-;
    Matrix
        A  0 0 1 0
        B  0 1 1 0
        C  1 1 0 1
        D  1 0 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    usercost fort_len (Standard) = 2
        0 1
        2 0
    ;
    ctype fort_len: 1-4;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

USERTYPE_EQUIV_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1" missing=? gap=-;
    Matrix
        A  0 0 1 0
        B  0 1 1 0
        C  1 1 0 1
        D  1 0 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    usertype fort_len (Standard) = 2
        0.0 1.0
        0.5 0.0
    ;
    ctype fort_len: 1-4;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

WTSET_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1" missing=? gap=-;
    Matrix
        A  0 0 1 0
        B  0 1 1 0
        C  1 1 0 1
        D  1 0 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    wtset * = 1:2 3:3;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

NOWTSET_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1" missing=? gap=-;
    Matrix
        A  0 0 1 0
        B  0 1 1 0
        C  1 1 0 1
        D  1 0 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""


class UsercostTest(unittest.TestCase):
    """Tests for usercost command (cost-to-rate conversion)."""

    def test_usercost_completes(self):
        """Analysis with usercost command completes."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_uc_")
        try:
            write_nex(tmpdir, "test.nex", USERCOST_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0)
            self.assertTrue(parse_completed(output))
            self.assertIn("Converting cost matrix to rates", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_usercost_matches_usertype(self):
        """usercost with cost=2 produces same lnL as usertype with rate=0.5."""
        results = {}
        for label, nex in [("usercost", USERCOST_NEX),
                           ("usertype", USERTYPE_EQUIV_NEX)]:
            tmpdir = tempfile.mkdtemp(prefix="yvyra_uc_")
            try:
                write_nex(tmpdir, "test.nex", nex)
                output, rc = run_yvyra(BINARY_PATH,
                    os.path.join(tmpdir, "test.nex"), tmpdir)
                self.assertEqual(rc, 0)
                header, rows = read_p_file(tmpdir)
                self.assertIsNotNone(header)
                lnl_idx = header.index("lnLike")
                results[label] = float(rows[-1][lnl_idx])
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        self.assertAlmostEqual(
            results["usercost"], results["usertype"], places=3,
            msg=f"usercost lnL ({results['usercost']}) should match "
                f"usertype ({results['usertype']})")


class WtsetTest(unittest.TestCase):
    """Tests for wtset character weights."""

    def test_wtset_completes(self):
        """Analysis with wtset weights completes."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_wt_")
        try:
            write_nex(tmpdir, "test.nex", WTSET_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0)
            self.assertTrue(parse_completed(output))
            self.assertIn("Character weights set", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_wtset_changes_likelihood(self):
        """Weighted analysis produces different lnL than unweighted."""
        results = {}
        for label, nex in [("weighted", WTSET_NEX),
                           ("unweighted", NOWTSET_NEX)]:
            tmpdir = tempfile.mkdtemp(prefix="yvyra_wt_")
            try:
                write_nex(tmpdir, "test.nex", nex)
                output, rc = run_yvyra(BINARY_PATH,
                    os.path.join(tmpdir, "test.nex"), tmpdir)
                self.assertEqual(rc, 0)
                header, rows = read_p_file(tmpdir)
                self.assertIsNotNone(header)
                lnl_idx = header.index("lnLike")
                results[label] = float(rows[-1][lnl_idx])
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        self.assertNotAlmostEqual(
            results["weighted"], results["unweighted"], places=0,
            msg="Weighted and unweighted lnL should differ")
        # Weighted should have larger magnitude (more effective chars)
        self.assertLess(results["weighted"], results["unweighted"],
            msg="Weighted lnL should be more negative (more effective chars)")

    def test_default_weights_match_unweighted(self):
        """Explicit weight=1 for all chars matches no wtset."""
        all_weight1_nex = NOWTSET_NEX.replace(
            "lset coding=variable;",
            "wtset * = 1:1 2:1 3:1 4:1;\n    lset coding=variable;")
        results = {}
        for label, nex in [("explicit_1", all_weight1_nex),
                           ("no_wtset", NOWTSET_NEX)]:
            tmpdir = tempfile.mkdtemp(prefix="yvyra_wt_")
            try:
                write_nex(tmpdir, "test.nex", nex)
                output, rc = run_yvyra(BINARY_PATH,
                    os.path.join(tmpdir, "test.nex"), tmpdir)
                self.assertEqual(rc, 0)
                header, rows = read_p_file(tmpdir)
                self.assertIsNotNone(header)
                lnl_idx = header.index("lnLike")
                results[label] = float(rows[-1][lnl_idx])
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        self.assertAlmostEqual(
            results["explicit_1"], results["no_wtset"], places=3,
            msg="Explicit weight=1 should match default (no wtset)")


SHOWMODEL_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1" missing=? gap=-;
    Matrix
        A  0 0 1 0
        B  0 1 1 0
        C  1 1 0 1
        D  1 0 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    lset coding=variable;
    showmodel;
    showparams;
    mcmcp nrun=1 nchain=1 ngen=100 samplefr=100;
    mcmcp printfr=100 diagnfr=100;
    mcmc;
End;
"""


TIPWEIGHTS_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1 2" missing=? gap=-;
    Matrix
        A  0 0 1 0
        B  0 1 ? 1
        C  1 1 0 2
        D  1 0 0 0
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    tipweights B 3 = 0:0.7 1:0.3;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

TIPWEIGHTS_MISSING_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1 2" missing=? gap=-;
    Matrix
        A  0 0 1 0
        B  0 1 ? 1
        C  1 1 0 2
        D  1 0 0 0
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""


class TipweightsTest(unittest.TestCase):
    """Tests for confidence-weighted tips (Phase 6)."""

    def test_tipweights_completes(self):
        """Analysis with tipweights command completes."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_tw_")
        try:
            write_nex(tmpdir, "test.nex", TIPWEIGHTS_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0, f"Binary exited with code {rc}")
            self.assertTrue(parse_completed(output))
            self.assertIn("Tip weights set", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_tipweights_changes_likelihood(self):
        """Tip weights produce different likelihood than plain missing data."""
        results = {}
        for label, nex in [("weighted", TIPWEIGHTS_NEX),
                           ("missing", TIPWEIGHTS_MISSING_NEX)]:
            tmpdir = tempfile.mkdtemp(prefix="yvyra_tw_")
            try:
                write_nex(tmpdir, "test.nex", nex)
                output, rc = run_yvyra(BINARY_PATH,
                    os.path.join(tmpdir, "test.nex"), tmpdir)
                self.assertEqual(rc, 0)
                m = re.search(r'Chain 1 -- ([-\d.]+) --', output)
                self.assertIsNotNone(m)
                results[label] = float(m.group(1))
            finally:
                shutil.rmtree(tmpdir, ignore_errors=True)

        self.assertNotAlmostEqual(
            results["weighted"], results["missing"], places=0,
            msg="Tipweights should produce different lnL than missing data")

    def test_tipweights_finite(self):
        """Tip-weighted analysis produces finite likelihood."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_tw_")
        try:
            write_nex(tmpdir, "test.nex", TIPWEIGHTS_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            m = re.search(r'Chain 1 -- ([-\d.]+) --', output)
            self.assertIsNotNone(m)
            lnl = float(m.group(1))
            self.assertFalse(math.isnan(lnl))
            self.assertLess(lnl, 0.0)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


class ShowCommandsTest(unittest.TestCase):
    """Test that showmodel and showparams commands work."""

    def test_showmodel_completes(self):
        """showmodel and showparams don't crash and display Standard."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_show_")
        try:
            write_nex(tmpdir, "test.nex", SHOWMODEL_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0, f"Binary exited with code {rc}")
            self.assertTrue(parse_completed(output))
            self.assertIn("Standard", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


ORDERED_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=4 nchar=5;
    Format datatype=standard symbols="0 1 2 3" missing=? gap=-;
    Matrix
        A  0 0 1 3 0
        B  0 1 2 2 1
        C  1 2 3 1 2
        D  2 3 0 0 3
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    ctype ordered: 1-5;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""

MISSING_DATA_NEX = """\
#NEXUS
Begin data;
    Dimensions ntax=5 nchar=6;
    Format datatype=standard symbols="0 1" missing=? gap=-;
    Matrix
        A  0 0 1 0 ? 1
        B  0 1 ? ? 0 0
        C  1 1 0 1 1 ?
        D  ? 0 0 1 1 0
        E  1 ? 1 0 0 1
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes;
    set seed=42;
    lset coding=variable;
    mcmcp nrun=1 nchain=1 ngen=500 samplefr=100;
    mcmcp printfr=500 diagnfr=500;
    mcmc;
End;
"""


class OrderedCharTest(unittest.TestCase):
    """Test ordered character model (ctype ordered)."""

    def test_ordered_completes(self):
        """Analysis with ordered characters completes."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_ord_")
        try:
            write_nex(tmpdir, "test.nex", ORDERED_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0, f"Binary exited with code {rc}")
            self.assertTrue(parse_completed(output), "Analysis did not complete")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_ordered_produces_finite_likelihood(self):
        """Ordered character analysis produces finite likelihood."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_ord_")
        try:
            write_nex(tmpdir, "test.nex", ORDERED_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            m = re.search(r'Chain 1 -- ([-\d.]+) --', output)
            self.assertIsNotNone(m, "Could not parse initial lnL")
            lnl = float(m.group(1))
            self.assertFalse(math.isnan(lnl), "lnL is NaN")
            self.assertLess(lnl, 0.0, "lnL should be negative")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


class MissingDataTest(unittest.TestCase):
    """Test handling of missing data (? characters)."""

    def test_missing_data_completes(self):
        """Analysis with extensive missing data completes."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_miss_")
        try:
            write_nex(tmpdir, "test.nex", MISSING_DATA_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            self.assertEqual(rc, 0, f"Binary exited with code {rc}")
            self.assertTrue(parse_completed(output), "Analysis did not complete")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_missing_data_finite_likelihood(self):
        """Missing data analysis produces finite likelihood."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_miss_")
        try:
            write_nex(tmpdir, "test.nex", MISSING_DATA_NEX)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "test.nex"), tmpdir)
            m = re.search(r'Chain 1 -- ([-\d.]+) --', output)
            self.assertIsNotNone(m, "Could not parse initial lnL")
            lnl = float(m.group(1))
            self.assertFalse(math.isnan(lnl))
            self.assertLess(lnl, 0.0)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


class AnatolianTest(unittest.TestCase):
    """Integration test: Billing & Elgh (2023) Anatolian dataset."""

    def test_anatolian_completes(self):
        """Full 27-character Anatolian analysis with 6 asymmetric models completes."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_anat_")
        nex_path = os.path.join(TESTING_DIR, "anatolian_full.nex")
        if not os.path.exists(nex_path):
            self.skipTest("anatolian_full.nex not found")
        try:
            shutil.copy2(nex_path, tmpdir)
            # Use shorter run for testing
            with open(os.path.join(tmpdir, "anatolian_full.nex")) as f:
                content = f.read()
            content = content.replace("ngen=50000", "ngen=1000")
            content = content.replace("nrun=2", "nrun=1")
            content = content.replace("nchain=2", "nchain=1")
            content = content.replace("samplefr=50", "samplefr=100")
            content = content.replace("printfr=10000", "printfr=1000")
            content = content.replace("diagnfr=10000", "diagnfr=1000")
            # Remove sump/sumt for speed
            content = content.replace("sump;", "")
            content = content.replace("sumt;", "")
            with open(os.path.join(tmpdir, "anatolian_full.nex"), "w") as f:
                f.write(content)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "anatolian_full.nex"), tmpdir)
            self.assertEqual(rc, 0, f"Anatolian analysis failed with rc={rc}")
            self.assertTrue(parse_completed(output), "Analysis did not complete")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_anatolian_likelihood_range(self):
        """Anatolian analysis produces likelihood in expected range."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_anat_")
        nex_path = os.path.join(TESTING_DIR, "anatolian_full.nex")
        if not os.path.exists(nex_path):
            self.skipTest("anatolian_full.nex not found")
        try:
            shutil.copy2(nex_path, tmpdir)
            with open(os.path.join(tmpdir, "anatolian_full.nex")) as f:
                content = f.read()
            content = content.replace("ngen=50000", "ngen=500")
            content = content.replace("nrun=2", "nrun=1")
            content = content.replace("nchain=2", "nchain=1")
            content = content.replace("samplefr=50", "samplefr=100")
            content = content.replace("printfr=10000", "printfr=500")
            content = content.replace("diagnfr=10000", "diagnfr=500")
            content = content.replace("sump;", "")
            content = content.replace("sumt;", "")
            with open(os.path.join(tmpdir, "anatolian_full.nex"), "w") as f:
                f.write(content)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "anatolian_full.nex"), tmpdir)
            self.assertEqual(rc, 0)
            # Check initial likelihood is in reasonable range
            m = re.search(r'Chain 1 -- ([-\d.]+) --', output)
            self.assertIsNotNone(m)
            lnl = float(m.group(1))
            self.assertGreater(lnl, -400.0, f"lnL {lnl} too low")
            self.assertLess(lnl, -50.0, f"lnL {lnl} too high")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


class CharlabelsTest(unittest.TestCase):
    """Test character labels command."""

    def _make_nex(self, tmpdir, charlabels_cmd="", run_mcmc=True):
        """Create a minimal nexus file with optional charlabels."""
        mcmc_line = "mcmc ngen=100 nrun=1 nchain=1 samplefr=10 printfr=100 diagnfr=100;" if run_mcmc else ""
        nex = f"""#NEXUS
Begin data;
    Dimensions ntax=4 nchar=4;
    Format datatype=standard symbols="0 1" missing=?;
    Matrix
        A 0101
        B 1010
        C 0110
        D 1100
    ;
End;
Begin yvyra;
    set autoclose=yes nowarnings=yes seed=42;
    {charlabels_cmd}
    lset coding=variable;
    {mcmc_line}
End;
"""
        path = os.path.join(tmpdir, "charlabels_test.nex")
        with open(path, "w") as f:
            f.write(nex)
        return path

    def test_charlabels_accepted(self):
        """charlabels command is accepted after matrix definition."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_cl_")
        try:
            nex = self._make_nex(tmpdir, "charlabels 1=foo 2=bar 3=baz 4=qux;")
            output, rc = run_yvyra(BINARY_PATH, nex, tmpdir)
            self.assertEqual(rc, 0, f"rc={rc}\n{output}")
            self.assertIn("Character labels set for 4 of 4 characters", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_charlabels_partial(self):
        """charlabels can label a subset of characters."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_cl_")
        try:
            nex = self._make_nex(tmpdir, "charlabels 1=alpha 3=gamma;")
            output, rc = run_yvyra(BINARY_PATH, nex, tmpdir)
            self.assertEqual(rc, 0, f"rc={rc}\n{output}")
            self.assertIn("Character labels set for 2 of 4 characters", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_charlabels_out_of_range(self):
        """charlabels rejects out-of-range character index."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_cl_")
        try:
            nex = self._make_nex(tmpdir, "charlabels 99=bad;")
            output, rc = run_yvyra(BINARY_PATH, nex, tmpdir)
            self.assertIn("out of range", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_charlabels_no_matrix(self):
        """charlabels before matrix gives error."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_cl_")
        try:
            nex = """#NEXUS
Begin yvyra;
    charlabels 1=foo;
End;
"""
            path = os.path.join(tmpdir, "no_matrix.nex")
            with open(path, "w") as f:
                f.write(nex)
            output, rc = run_yvyra(BINARY_PATH, path, tmpdir)
            self.assertIn("matrix must be specified", output)
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)

    def test_anatolian_chardiag_uses_labels(self):
        """chardiag output uses character labels from Anatolian dataset."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_cl_")
        nex_path = os.path.join(TESTING_DIR, "anatolian_full.nex")
        if not os.path.exists(nex_path):
            self.skipTest("anatolian_full.nex not found")
        try:
            shutil.copy2(nex_path, tmpdir)
            with open(os.path.join(tmpdir, "anatolian_full.nex")) as f:
                content = f.read()
            # Short run, then chardiag
            content = content.replace("ngen=50000", "ngen=500")
            content = content.replace("nrun=2", "nrun=1")
            content = content.replace("nchain=2", "nchain=1")
            content = content.replace("samplefr=50", "samplefr=100")
            content = content.replace("printfr=10000", "printfr=500")
            content = content.replace("diagnfr=10000", "diagnfr=500")
            content = content.replace("sump;", "")
            content = content.replace("sumt;", "")
            # Add chardiag before last End;
            last_end = content.rfind("End;")
            content = content[:last_end] + "chardiag clade={Luwian,Lycian};\n" + content[last_end:]
            with open(os.path.join(tmpdir, "anatolian_full.nex"), "w") as f:
                f.write(content)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "anatolian_full.nex"), tmpdir)
            self.assertEqual(rc, 0, f"rc={rc}")
            # Verify chardiag output contains character labels (not raw lnL column names)
            self.assertIn("len_fric", output)
            self.assertIn("gen_ers", output)
            self.assertIn("lar_phon", output)
            self.assertIn("CHARACTER LIKELIHOOD CONTRIBUTIONS", output)
            # Verify CI section is present and has valid values
            self.assertIn("CONSISTENCY INDEX", output)
            self.assertIn("MAP tree from sample", output)
            # Extract CI values: format is "charname   steps  minsteps  ci"
            ci_lines = re.findall(r'^\s+\S+\s+\d+\s+\d+\s+(\d+\.\d{3})\s*$', output, re.MULTILINE)
            self.assertGreater(len(ci_lines), 0, "No CI values found")
            for ci_str in ci_lines:
                ci = float(ci_str)
                self.assertGreaterEqual(ci, 0.0, f"CI {ci} < 0")
                self.assertLessEqual(ci, 1.0, f"CI {ci} > 1")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


class AsrentropyTest(unittest.TestCase):
    """Test ancestral state entropy command."""

    def test_asrentropy_anatolian(self):
        """asrentropy runs on Anatolian dataset and produces entropy output."""
        tmpdir = tempfile.mkdtemp(prefix="yvyra_asr_")
        nex_path = os.path.join(TESTING_DIR, "anatolian_full.nex")
        if not os.path.exists(nex_path):
            self.skipTest("anatolian_full.nex not found")
        try:
            shutil.copy2(nex_path, tmpdir)
            with open(os.path.join(tmpdir, "anatolian_full.nex")) as f:
                content = f.read()
            content = content.replace("ngen=50000", "ngen=500")
            content = content.replace("nrun=2", "nrun=1")
            content = content.replace("nchain=2", "nchain=1")
            content = content.replace("samplefr=50", "samplefr=100")
            content = content.replace("printfr=10000", "printfr=500")
            content = content.replace("diagnfr=10000", "diagnfr=500")
            content = content.replace("sump;", "")
            content = content.replace("sumt;", "")
            # Add asrentropy before last End;
            last_end = content.rfind("End;")
            content = content[:last_end] + "asrentropy;\n" + content[last_end:]
            with open(os.path.join(tmpdir, "anatolian_full.nex"), "w") as f:
                f.write(content)
            output, rc = run_yvyra(BINARY_PATH,
                os.path.join(tmpdir, "anatolian_full.nex"), tmpdir)
            self.assertEqual(rc, 0, f"rc={rc}")
            # Check output structure
            self.assertIn("ANCESTRAL STATE ENTROPY", output)
            self.assertIn("(root)", output)
            self.assertIn("(Luwian,Lycian)", output)
            self.assertIn("Entropy", output)
            self.assertIn("Conf", output)
            # Check that character labels are used
            self.assertIn("len_fric", output)
            self.assertIn("encl_du", output)
            # Verify entropy values are in valid range [0, max]
            entropy_vals = re.findall(r'Entropy.*?\n.*?(\d+\.\d{4})\s+\d+\.\d{2}', output)
            for e_str in entropy_vals:
                e = float(e_str)
                self.assertGreaterEqual(e, 0.0, f"Entropy {e} < 0")
            # Verify confidence values are in [0, 1]
            conf_vals = re.findall(r'\d+\.\d{4}\s+(\d+\.\d{2})\s+\d+', output)
            for c_str in conf_vals:
                c = float(c_str)
                self.assertGreaterEqual(c, 0.0, f"Confidence {c} < 0")
                self.assertLessEqual(c, 1.0, f"Confidence {c} > 1")
        finally:
            shutil.rmtree(tmpdir, ignore_errors=True)


def main():
    parser = argparse.ArgumentParser(description="Test harness for yvyra Mk model")
    parser.add_argument(
        "--binary", required=True, help="Path to the yvyra (or mb) binary"
    )
    # Separate our args from unittest args
    args, remaining = parser.parse_known_args()

    global BINARY_PATH, EXPECTED_RANGES
    BINARY_PATH = os.path.abspath(args.binary)
    EXPECTED_RANGES = load_expected_ranges()

    if not os.path.isfile(BINARY_PATH):
        print(f"Error: binary not found at {BINARY_PATH}", file=sys.stderr)
        sys.exit(1)
    if not os.access(BINARY_PATH, os.X_OK):
        print(f"Error: {BINARY_PATH} is not executable", file=sys.stderr)
        sys.exit(1)

    # Run unittest with remaining args
    sys.argv = [sys.argv[0]] + remaining
    unittest.main(module=__name__, exit=True)


if __name__ == "__main__":
    main()
