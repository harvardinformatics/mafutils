import os
import subprocess
import sys

import pytest

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))
MAF_FILE = os.path.join(TEST_DIR, "example.maf")

HEADER = "BED check failed"


def run_fetch(args):
    return subprocess.run([sys.executable, "-m", "mafutils", "fetch"] + args,
                          check=False, cwd=REPO_ROOT, capture_output=True, text=True)


def fetch_with_bed(tmp_path, lines, *extra):
    bed = os.path.join(str(tmp_path), "r.bed")
    with open(bed, "w") as fp:
        fp.write("".join(line + "\n" for line in lines))
    out = os.path.join(str(tmp_path), "out")
    r = run_fetch([MAF_FILE, bed, "-o", out] + list(extra))
    return r, r.stdout + r.stderr, out


def assert_failed_before_any_work(r, combined, out):
    assert r.returncode == 1, combined
    assert combined.count(HEADER) == 1               # one error, not one per region
    assert "Traceback" not in combined
    assert sorted(os.listdir(out)) == ["maf_fetch.log"]


# ---------------------------------------------------------------------------
# Unusable lines -- these used to crash (the skip path passed file= to a logger)
# or, for bad intervals, pass silently as "no overlapping blocks"
# ---------------------------------------------------------------------------


def test_non_integer_coordinates_are_reported_with_line_number(tmp_path):
    r, combined, out = fetch_with_bed(tmp_path, ["chr1\t2\t6\tok", "chr1\tabc\t12\tbad"])
    assert_failed_before_any_work(r, combined, out)
    assert "e.g. line 2: start/end not integers" in combined


def test_too_few_columns_is_reported(tmp_path):
    r, combined, out = fetch_with_bed(tmp_path, ["chr1\t2\t6\tok", "chr1\t8"])
    assert_failed_before_any_work(r, combined, out)
    assert "e.g. line 2: fewer than 3 columns" in combined


@pytest.mark.parametrize("line", ["chr1\t20\t6\tinverted", "chr1\t5\t5\tempty", "chr1\t-1\t5\tnegative"])
def test_invalid_intervals_are_reported(line, tmp_path):
    r, combined, out = fetch_with_bed(tmp_path, ["chr2\t5\t10\tok", line])
    assert_failed_before_any_work(r, combined, out)
    assert "e.g. line 2: invalid interval" in combined


def test_track_and_browser_header_lines_are_skipped(tmp_path):
    r, combined, _out = fetch_with_bed(
        tmp_path, ["browser position chr1:1-100", "track name=peaks", "# comment", "chr2\t5\t10\tok"], "-b", "id")
    assert r.returncode == 0, combined
    assert HEADER not in combined


# ---------------------------------------------------------------------------
# -b id: every region needs a name
# ---------------------------------------------------------------------------


def test_basename_id_without_names_fails_once_before_any_work(tmp_path):
    """
    Regression: each unnamed region used to log its own error, then fetch ran
    anyway and a worker crashed with KeyError('output_basename'), leaving a
    half-written output directory.
    """
    r, combined, out = fetch_with_bed(
        tmp_path, ["chr1\t2\t6", "chr1\t8\t12\tnamed", "chr2\t5\t10"], "-b", "id", "-p", "2")
    assert_failed_before_any_work(r, combined, out)
    assert "2 of 3 regions have no name (4th column), which -b id needs, e.g. chr1:2-6." in combined


def test_bed_with_no_names_at_all_gets_one_short_line(tmp_path):
    r, combined, out = fetch_with_bed(tmp_path, ["chr1\t2\t6", "chr2\t5\t10"], "-b", "id")
    assert_failed_before_any_work(r, combined, out)
    assert "r.bed: no region has a name (4th column), which -b id needs. Nothing was fetched." in combined


def test_unnamed_regions_are_fine_without_basename_id(tmp_path):
    r, combined, _out = fetch_with_bed(tmp_path, ["chr1\t2\t6", "chr2\t5\t10"], "-b", "coords")
    assert r.returncode == 0, combined


# ---------------------------------------------------------------------------
# Output-file collisions: the later region would silently replace the earlier
# ---------------------------------------------------------------------------


def test_duplicate_names_with_basename_id_are_reported(tmp_path):
    r, combined, out = fetch_with_bed(
        tmp_path, ["chr1\t2\t6\tsame", "chr2\t5\t10\tsame", "chr3\t4\t8\tsame", "chr4\t4\t8\tunique"], "-b", "id")
    assert_failed_before_any_work(r, combined, out)
    assert "1 name is shared by several regions" in combined
    assert "same (3 regions)" in combined


def test_same_name_on_different_scaffolds_is_fine_with_scaffold_subdirs(tmp_path):
    lines = ["chr1\t2\t6\tsame", "chr2\t5\t10\tsame"]
    r, combined, out = fetch_with_bed(tmp_path, lines, "-b", "id", "--scaffold-subdirs")
    assert r.returncode == 0, combined
    assert os.path.isfile(os.path.join(out, "chr1", "same.maf"))
    assert os.path.isfile(os.path.join(out, "chr2", "same.maf"))


def test_same_name_on_the_same_scaffold_collides_even_with_scaffold_subdirs(tmp_path):
    r, combined, out = fetch_with_bed(
        tmp_path, ["chr1\t2\t6\tsame", "chr1\t8\t12\tsame"], "-b", "id", "--scaffold-subdirs")
    assert_failed_before_any_work(r, combined, out)
    assert "same (2 regions)" in combined


def test_identical_regions_with_basename_coords_are_reported(tmp_path):
    r, combined, out = fetch_with_bed(tmp_path, ["chr1\t2\t6", "chr2\t5\t10", "chr1\t2\t6"], "-b", "coords")
    assert_failed_before_any_work(r, combined, out)
    assert "1 region is listed more than once" in combined
    assert "chr1:2-6 (2 times)" in combined


def test_identical_regions_are_fine_with_basename_count(tmp_path):
    """count numbers every region, so duplicates get distinct files."""
    r, combined, _out = fetch_with_bed(tmp_path, ["chr1\t2\t6", "chr1\t2\t6"], "-b", "count")
    assert r.returncode == 0, combined


# ---------------------------------------------------------------------------
# Everything at once: one error listing every problem
# ---------------------------------------------------------------------------


def test_all_problems_are_reported_together(tmp_path):
    r, combined, out = fetch_with_bed(
        tmp_path,
        ["chr1\tabc\t6\tbad", "chr1\t9\t3\tinverted", "chr1\t2\t6", "chr2\t5\t10\tdup", "chr3\t4\t8\tdup"],
        "-b", "id")
    assert_failed_before_any_work(r, combined, out)
    assert "nothing was fetched:\n  - 2 unusable lines, e.g. line 1: start/end not integers" in combined
    assert "\n  - 1 of 3 regions have no name" in combined
    assert "\n  - 1 name is shared by several regions" in combined
    assert "dup (2 regions)" in combined


def test_clean_bed_logs_that_the_check_passed(tmp_path):
    r, combined, _out = fetch_with_bed(tmp_path, ["chr1\t2\t6\ta", "chr2\t5\t10\tb"], "-b", "id")
    assert r.returncode == 0, combined
    assert "BED check passed: 2 regions." in combined


def test_bed_with_no_data_lines_is_reported_not_crashed(tmp_path):
    """An empty BED used to crash on batching (range() with a zero step)."""
    r, combined, out = fetch_with_bed(tmp_path, ["# only a comment", "track name=x"])
    assert_failed_before_any_work(r, combined, out)
    assert "no usable regions" in combined
