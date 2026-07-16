import csv
import os
import subprocess
import sys

import pytest

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

# Space-delimited fixture (shared with test_fetch.py).
MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")

# Tab-delimited fixture: same alignment content as example.maf, but with
# tab-delimited 's' lines. Regression fixture for a bug where gc's line
# parser assumed a literal single-space delimiter (`line.startswith("s ")`)
# and silently found zero sequence lines on tab-delimited real-world MAFs.
TAB_MAF_FILE = os.path.join(TEST_DIR, "example-tabs.maf")
TAB_INDEX_FILE = os.path.join(TEST_DIR, "example-tabs.maf.block.idx")

EXPECTED_GC = {
    "human": 0.44827586,
    "chimp": 0.56097561,
    "gorilla": 0.25000000,
}
EXPECTED_MEAN_GC = sum(EXPECTED_GC.values()) / len(EXPECTED_GC)


def read_gc_csv(path):
    with open(path, "r", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        return {row["species"]: float(row["gc"]) for row in reader}


def run_gc(args, tmp_path):
    cmd = [sys.executable, "-m", "mafutils", "gc"] + args
    result = subprocess.run(cmd, check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, f"gc failed (exit {result.returncode}):\n{result.stderr}"


def assert_gc_matches_expected(gc_csv_path):
    got = read_gc_csv(gc_csv_path)
    assert got.keys() == EXPECTED_GC.keys()
    for species, expected_value in EXPECTED_GC.items():
        assert got[species] == pytest.approx(expected_value, abs=1e-6), (
            f"{species}: expected {expected_value}, got {got[species]}"
        )


def test_gc_sequential_space_delimited(tmp_path):
    """Sequential path (no index) on the standard space-delimited fixture."""
    out_prefix = os.path.join(str(tmp_path), "out")
    run_gc([MAF_FILE, "-o", out_prefix], tmp_path)
    assert_gc_matches_expected(out_prefix + ".gc.csv")


def test_gc_parallel_space_delimited(tmp_path):
    """Parallel path (index + processes>1) on the standard space-delimited fixture."""
    out_prefix = os.path.join(str(tmp_path), "out")
    run_gc([MAF_FILE, INDEX_FILE, "-o", out_prefix, "-p", "2"], tmp_path)
    assert_gc_matches_expected(out_prefix + ".gc.csv")


def test_gc_mean_file(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    run_gc([MAF_FILE, "-o", out_prefix], tmp_path)
    with open(out_prefix + ".gc.mean.txt", "r", encoding="utf-8") as fp:
        mean_value = float(fp.read().strip())
    assert mean_value == pytest.approx(EXPECTED_MEAN_GC, abs=1e-6)


def test_gc_sequential_tab_delimited(tmp_path):
    """
    Regression test: sequential path on a tab-delimited MAF (no index) must
    find the same species/GC values as the space-delimited fixture, since
    the underlying alignment content is identical.
    """
    out_prefix = os.path.join(str(tmp_path), "out")
    run_gc([TAB_MAF_FILE, "-o", out_prefix], tmp_path)
    assert_gc_matches_expected(out_prefix + ".gc.csv")


def test_gc_parallel_tab_delimited_falls_back_and_matches(tmp_path):
    """
    Regression test: even when a tab-delimited MAF is uncompressed and an
    index + --processes>1 are given (so the parallel path runs), results
    must still match, since both paths share the same line-parsing logic.
    """
    out_prefix = os.path.join(str(tmp_path), "out")
    run_gc([TAB_MAF_FILE, TAB_INDEX_FILE, "-o", out_prefix, "-p", "2"], tmp_path)
    assert_gc_matches_expected(out_prefix + ".gc.csv")
