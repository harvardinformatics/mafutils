"""
Index format 3: the per-block sequence count (block.idx column 6) and the
seq_lines / max_seqs header totals count 's' lines only. Formats 1-2 also
counted i/e/q lines, which inflated the count on MAFs that have them.

Fixture: multiz45way-excerpt.maf -- the `##maf` line plus 14 consecutive real
blocks (galGal6 chr6:20883553-20891747) of a 45-way multiz MAF, with 72 s, 58 i
and 12 e lines. It includes the block at chr6:20884758 (2 s + 1 i + 1 e lines)
that a format-2 index reported as 4 sequences.

Malformed inputs for the error tests are built from copies of the fixture in
tmp_path; nothing here is invented sequence data.
"""

import csv
import gzip
import os
import subprocess
import sys

import pytest

from mafutils.lib import bgzf
from mafutils.lib import common as COMMON

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(TEST_DIR)
MAF_FILE = os.path.join(TEST_DIR, "multiz45way-excerpt.maf")

REPORTED_BLOCK_START = "20884758"


def run(args):
    return subprocess.run([sys.executable, "-m", "mafutils"] + args,
                          check=False, cwd=REPO_ROOT, capture_output=True, text=True)


def maf_text():
    with open(MAF_FILE, "r", encoding="utf-8") as fp:
        return fp.read()


def s_counts_per_block(text):
    """Independent count: lines whose first field is 's', per 'a' block."""
    counts = []
    for line in text.splitlines():
        fields = line.split()
        if not fields:
            continue
        if fields[0] == "a":
            counts.append(0)
        elif fields[0] == "s":
            counts[-1] += 1
    return counts


def read_index(path):
    with open(path, "r", encoding="utf-8") as fp:
        lines = fp.read().splitlines()
    return COMMON.readIndexHeader(path), [line.split("\t") for line in lines[1:]]


def build_index(maf_path):
    result = run(["index", maf_path, "--progress-interval", "0"])
    assert result.returncode == 0, result.stdout + result.stderr
    return COMMON.deriveBlockIndexPath(maf_path), COMMON.deriveScaffoldIndexPath(maf_path)


def copy_as(tmp_path, compression):
    data = maf_text().encode("ascii")
    if compression == "none":
        path = os.path.join(str(tmp_path), "bird.maf")
        with open(path, "wb") as fp:
            fp.write(data)
    elif compression == "gz":
        path = os.path.join(str(tmp_path), "bird.maf.gz")
        with gzip.open(path, "wb") as fp:
            fp.write(data)
    else:
        path = os.path.join(str(tmp_path), "bird.maf.bgz")
        with bgzf.BgzfWriter(path, "wb") as fp:
            fp.write(data)
    assert COMMON.detectCompression(path) == compression
    return path


# ---------------------------------------------------------------------------
# Counts
# ---------------------------------------------------------------------------


def test_fixture_has_the_line_types_this_file_is_about():
    kinds = [line[:1] for line in maf_text().splitlines() if line]
    assert kinds.count("s") == 72 and kinds.count("i") == 58 and kinds.count("e") == 12


@pytest.mark.parametrize("compression", ["none", "gz", "bgzip"])
def test_count_is_s_lines_only(tmp_path, compression):
    block_idx, _ = build_index(copy_as(tmp_path, compression))
    header, rows = read_index(block_idx)
    expected = s_counts_per_block(maf_text())

    assert [int(r[5]) for r in rows] == expected
    assert header["format"] == "3"
    assert int(header["seq_lines"]) == sum(expected) == 72
    assert int(header["max_seqs"]) == max(expected)
    assert int(header["blocks"]) == len(expected) == 14
    assert [int(r[5]) for r in rows if r[1] == REPORTED_BLOCK_START] == [2]


def test_offsets_bracket_blocks_with_that_many_s_lines(tmp_path):
    maf = copy_as(tmp_path, "none")
    block_idx, _ = build_index(maf)
    _, rows = read_index(block_idx)
    with open(maf, "rb") as fp:
        data = fp.read()
    for row in rows:
        block = data[int(row[6]):int(row[7])]
        assert block.startswith(b"a ")
        assert sum(line.startswith(b"s") for line in block.split(b"\n")) == int(row[5])


def test_comment_and_blank_lines_do_not_change_the_counts(tmp_path):
    plain = copy_as(tmp_path, "none")
    noisy = os.path.join(str(tmp_path), "noisy.maf")
    out = []
    for line in maf_text().splitlines(keepends=True):
        out.append(line)
        if line.startswith("s "):
            out.append("# a comment inside a block\n")
        elif not line.strip():
            out.append("\n\n")
    with open(noisy, "w", encoding="utf-8") as fp:
        fp.writelines(out)

    _, plain_rows = read_index(build_index(plain)[0])
    _, noisy_rows = read_index(build_index(noisy)[0])
    # Offsets (columns 7-8) move with the inserted bytes; everything else must not.
    assert [r[:6] for r in noisy_rows] == [r[:6] for r in plain_rows]


# ---------------------------------------------------------------------------
# Malformed MAFs: one clear error, and no index left behind
# ---------------------------------------------------------------------------


def _track_line_on_top(lines):
    return lines[:1] + ["track name=birds\n"] + lines[1:]


def _indented_s_line(lines):
    i = [n for n, line in enumerate(lines) if line.startswith("s ")][1]
    return lines[:i] + [" " + lines[i]] + lines[i + 1:]


def _s_line_before_first_block(lines):
    first_s = next(line for line in lines if line.startswith("s "))
    return lines[:1] + [first_s] + lines[1:]


def _block_without_s_lines(lines):
    a_lines = [n for n, line in enumerate(lines) if line.startswith("a ")]
    start, end = a_lines[2], a_lines[3]
    return lines[:start] + [l for l in lines[start:end] if not l.startswith("s ")] + lines[end:]


@pytest.mark.parametrize("mutate, message", [
    (_track_line_on_top, "unrecognized line before the first block: 'track name=birds'"),
    (_indented_s_line, "unrecognized line in block 1:"),
    (_s_line_before_first_block, "'s' line before the first 'a' line"),
    (_block_without_s_lines, "block 3 has no 's' line"),
], ids=["track-line", "indented-s-line", "s-before-first-a", "block-without-s"])
def test_malformed_maf_errors_and_writes_no_index(tmp_path, mutate, message):
    maf = os.path.join(str(tmp_path), "bad.maf")
    with open(maf, "w", encoding="utf-8") as fp:
        fp.writelines(mutate(maf_text().splitlines(keepends=True)))

    result = run(["index", maf, "--progress-interval", "0"])
    combined = result.stdout + result.stderr
    assert result.returncode == 1, combined
    assert message in combined
    assert "No index was written" in combined
    assert "Traceback" not in combined
    assert not os.path.exists(COMMON.deriveBlockIndexPath(maf))
    assert not os.path.exists(COMMON.deriveScaffoldIndexPath(maf))


# ---------------------------------------------------------------------------
# Older indexes (format 2 / headerless): refused only where the count is shown
# ---------------------------------------------------------------------------


def old_index_maf(tmp_path, kind):
    """A MAF copy whose fresh index pair is rewritten as format 2, or headerless."""
    maf = copy_as(tmp_path, "none")
    for path in build_index(maf):
        with open(path, "r", encoding="utf-8") as fp:
            lines = fp.readlines()
        if kind == "format2":
            assert " format=3 " in lines[0]
            lines[0] = lines[0].replace(" format=3 ", " format=2 ")
        else:
            lines = lines[1:]
        with open(path, "w", encoding="utf-8") as fp:
            fp.writelines(lines)
    return maf


@pytest.mark.parametrize("kind, described", [
    ("format2", "index format 2"),
    ("headerless", "no index header"),
])
def test_stats_refuses_an_old_index_unless_the_block_table_is_off(tmp_path, kind, described):
    maf = old_index_maf(tmp_path, kind)
    out = os.path.join(str(tmp_path), "out")

    refused = run(["stats", maf, "-o", out, "--progress-interval", "0"])
    combined = refused.stdout + refused.stderr
    assert refused.returncode == 1, combined
    assert described in combined and "--no-block-table" in combined
    assert "Rebuild it with `mafutils index" in combined
    assert not os.path.exists(out + ".overall.tsv")

    allowed = run(["stats", maf, "-o", out, "--no-block-table", "--progress-interval", "0"])
    assert allowed.returncode == 0, allowed.stdout + allowed.stderr
    assert os.path.exists(out + ".overall.tsv")


@pytest.mark.parametrize("kind, described", [
    ("format2", "index format 2"),
    ("headerless", "a headerless index"),
])
def test_info_marks_old_counts_unavailable_and_exits_1(tmp_path, kind, described):
    maf = old_index_maf(tmp_path, kind)
    result = run(["info", maf, "--sample-blocks", "0"])
    assert result.returncode == 1, result.stdout + result.stderr
    assert "blocks         : 14" in result.stdout
    assert f"seq lines      : unavailable ({described} may count i/e/q lines)" in result.stdout
    assert f"max seqs/block : unavailable ({described} may count i/e/q lines)" in result.stdout
    assert "Rebuild it with `mafutils index" in result.stderr


def test_validate_reports_an_old_format_without_changing_the_verdict(tmp_path):
    maf = old_index_maf(tmp_path, "format2")
    result = run(["validate", maf])
    combined = result.stdout + result.stderr
    assert result.returncode == 0, combined
    assert "VERIFIED" in combined
    assert "index format 2" in combined


def test_validate_reports_a_current_format(tmp_path):
    maf = copy_as(tmp_path, "none")
    build_index(maf)
    result = run(["validate", maf])
    assert result.returncode == 0, result.stdout + result.stderr
    assert "format     : 3" in result.stdout + result.stderr


def test_fetch_and_gc_accept_a_format_2_index(tmp_path):
    maf = old_index_maf(tmp_path, "format2")
    bed = os.path.join(str(tmp_path), "region.bed")
    with open(bed, "w", encoding="utf-8") as fp:
        fp.write("chr6\t20884000\t20885000\tr1\n")

    fetched = run(["fetch", maf, bed, "-b", "id", "-o", os.path.join(str(tmp_path), "fetch"),
                   "--progress-interval", "0"])
    assert fetched.returncode == 0, fetched.stdout + fetched.stderr
    assert os.path.isfile(os.path.join(str(tmp_path), "fetch", "r1.maf"))

    gc = run(["gc", maf, "-o", os.path.join(str(tmp_path), "gc"), "--progress-interval", "0"])
    assert gc.returncode == 0, gc.stdout + gc.stderr

    for result in (fetched, gc):
        assert "i/e/q" not in result.stdout + result.stderr


# ---------------------------------------------------------------------------
# stats: with a format-3 index, its count and stats' own parse agree
# ---------------------------------------------------------------------------


def test_stats_index_and_parsed_counts_agree(tmp_path):
    maf = copy_as(tmp_path, "none")
    build_index(maf)
    out = os.path.join(str(tmp_path), "out")
    result = run(["stats", maf, "-o", out, "--progress-interval", "0"])
    assert result.returncode == 0, result.stdout + result.stderr

    with open(out + ".block.tsv", "r", encoding="utf-8") as fp:
        rows = list(csv.DictReader(fp, delimiter="\t"))
    assert len(rows) == 14
    assert [int(r["n_seq_lines_index"]) for r in rows] == [int(r["n_seq_lines_parsed"]) for r in rows]
    assert sum(int(r["n_seq_lines_parsed"]) for r in rows) == 72
