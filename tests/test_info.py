import os
import re
import subprocess
import sys

import pytest

from mafutils.lib import common as COMMON

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "example.maf")
TABS_MAF = os.path.join(TEST_DIR, "example-tabs.maf")

# Hand-verifiable against tests/example.maf, and cross-checked against the
# independent values already pinned in test_stats.py (10 blocks, 61 alignment
# columns, 20 sequence lines).
EXPECTED = {
    "blocks": 10,
    "scaffolds": 5,      # chr1..chr4, chrX -- chr4 appears in two runs, one scaffold
    "ref_bases": 59,
    "aln_cols": 61,
    "seq_lines": 20,
    "max_seqs": 2,
}


def run_info(args):
    return subprocess.run(
        [sys.executable, "-m", "mafutils", "info"] + args,
        check=False, cwd=REPO_ROOT, capture_output=True, text=True,
    )


def build_index(tmp_path, maf=MAF_FILE, name="t.maf"):
    """A freshly built index, i.e. one whose header carries the aggregates."""
    import shutil
    copy = os.path.join(str(tmp_path), name)
    shutil.copyfile(maf, copy)
    r = subprocess.run([sys.executable, "-m", "mafutils", "index", copy],
                       check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return copy


def value(output, label):
    m = re.search(rf"^{re.escape(label)}\s*:\s*([0-9,]+)", output, re.M)
    assert m, f"{label!r} not found in:\n{output}"
    return int(m.group(1).replace(",", ""))


# ---------------------------------------------------------------------------
# Index header aggregates
# ---------------------------------------------------------------------------


def test_index_records_aggregates_in_the_header(tmp_path):
    maf = build_index(tmp_path)
    header = COMMON.readIndexHeader(COMMON.deriveBlockIndexPath(maf))
    for key, expected in EXPECTED.items():
        assert int(header[key]) == expected, key


def test_both_index_headers_are_identical(tmp_path):
    """
    validate cross-checks the two headers with exact equality, so aggregates
    present in only one would make every future validate report MISMATCH.
    """
    maf = build_index(tmp_path)
    block = open(COMMON.deriveBlockIndexPath(maf)).readline()
    scaffold = open(COMMON.deriveScaffoldIndexPath(maf)).readline()
    assert block == scaffold


def test_validate_still_verifies_with_aggregates(tmp_path):
    maf = build_index(tmp_path)
    r = subprocess.run([sys.executable, "-m", "mafutils", "validate", maf],
                       check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stdout + r.stderr
    assert "VERIFIED" in (r.stdout + r.stderr)


def test_format_version_unchanged(tmp_path):
    """Additive keys only -- bumping the format would imply an incompatibility
    that does not exist, and old readers already ignore unknown keys."""
    maf = build_index(tmp_path)
    header = COMMON.readIndexHeader(COMMON.deriveBlockIndexPath(maf))
    assert header["format"] == "2"


# ---------------------------------------------------------------------------
# info: the three tiers
# ---------------------------------------------------------------------------


def test_info_reports_exact_values_from_the_header(tmp_path):
    maf = build_index(tmp_path)
    r = run_info([maf])
    assert r.returncode == 0, r.stderr
    out = r.stdout
    assert value(out, "blocks") == EXPECTED["blocks"]
    assert value(out, "ref scaffolds") == EXPECTED["scaffolds"]
    assert value(out, "ref bases") == EXPECTED["ref_bases"]
    assert value(out, "aln columns") == EXPECTED["aln_cols"]
    assert value(out, "seq lines") == EXPECTED["seq_lines"]
    assert value(out, "max seqs/block") == EXPECTED["max_seqs"]
    # no fallback warning when the header already has them
    assert "predates" not in (r.stdout + r.stderr)


def test_info_falls_back_to_scanning_an_older_index(tmp_path):
    """
    An index without the aggregate keys must still produce the right numbers,
    warn that it is scanning, and name the index size so a long wait on a large
    index is never a surprise.
    """
    maf = build_index(tmp_path)
    index_path = COMMON.deriveBlockIndexPath(maf)
    lines = open(index_path).readlines()
    # strip the aggregate keys, leaving a v0.6.0-shaped header
    lines[0] = re.sub(r"\s+(blocks|scaffolds|ref_bases|aln_cols|seq_lines|max_seqs)=\d+", "", lines[0])
    open(index_path, "w").writelines(lines)
    assert "blocks=" not in open(index_path).readline()

    r = run_info([maf])
    assert r.returncode == 0, r.stderr
    combined = r.stdout + r.stderr
    assert "predates the aggregate header fields" in combined
    assert "Rebuild with `mafutils index`" in combined
    assert re.search(r"scanning the index \([\d.]+ \w+\)", combined), combined
    # and the scanned numbers equal what the header would have said
    for label, key in (("blocks", "blocks"), ("ref bases", "ref_bases"),
                       ("aln columns", "aln_cols"), ("seq lines", "seq_lines")):
        assert value(r.stdout, label) == EXPECTED[key], label


def test_info_without_an_index_says_so_and_does_not_scan(tmp_path):
    import shutil
    copy = os.path.join(str(tmp_path), "lonely.maf")
    shutil.copyfile(MAF_FILE, copy)

    r = run_info([copy])
    assert r.returncode == 0, r.stderr
    out = r.stdout
    assert "index" in out and "not found" in out
    assert "unavailable without an index" in out
    assert "run `mafutils index` first" in out
    # species still works -- it reads the MAF directly, not the index
    assert value(out, "species") == 3


# ---------------------------------------------------------------------------
# Species sampling
# ---------------------------------------------------------------------------


def test_species_are_labelled_as_a_sample(tmp_path):
    maf = build_index(tmp_path)
    out = run_info([maf]).stdout
    assert value(out, "species") == 3
    assert "possibly incomplete" in out
    assert "sampled from the first" in out
    for name in ("human", "chimp", "gorilla"):
        assert name in out
    assert "mafutils stats" in out


def test_sample_blocks_zero_skips_the_species_section(tmp_path):
    maf = build_index(tmp_path)
    out = run_info([maf, "--sample-blocks", "0"]).stdout
    # "species" still appears in the seq-lines note, so assert on the section
    assert "possibly incomplete" not in out
    assert not re.search(r"^species\s*:", out, re.M), out
    assert "mafutils stats" not in out
    assert value(out, "blocks") == EXPECTED["blocks"]


def test_tab_delimited_maf_yields_the_same_species(tmp_path):
    """
    Real MAFs are tab-delimited; detecting sequence lines with
    startswith("s ") silently matches nothing on them, a false negative this
    project has hit before.
    """
    os.makedirs(tmp_path / "tabs", exist_ok=True)
    tabs = build_index(tmp_path / "tabs", maf=TABS_MAF, name="tabs.maf")
    out = run_info([tabs]).stdout
    assert value(out, "species") == 3
    for name in ("human", "chimp", "gorilla"):
        assert name in out


def test_negative_sample_blocks_is_rejected(tmp_path):
    maf = build_index(tmp_path)
    r = run_info([maf, "--sample-blocks", "-1"])
    assert r.returncode != 0
