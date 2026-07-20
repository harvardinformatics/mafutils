import os
import shutil
import subprocess
import sys
import time

import pytest

from mafutils.lib import common as COMMON

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")

# Deliberately preserved pre-hash-header fixture (see DEVELOPMENT.md) -- used
# directly as a "headerless for these fields" case rather than constructing
# one by hand.
HEADERLESS_SCAFFOLD_INDEX = os.path.join(TEST_DIR, "example.maf.scaffold.idx")


def run_module(args, cwd=REPO_ROOT):
    cmd = [sys.executable, "-m", "mafutils"] + args
    return subprocess.run(cmd, check=False, cwd=cwd, capture_output=True, text=True)


def combined_output(result):
    return result.stdout + result.stderr


# ---------------------------------------------------------------------------
# Header contents
# ---------------------------------------------------------------------------


def test_index_header_has_size_mtime_hash():
    header = COMMON.readIndexHeader(INDEX_FILE)
    assert header is not None
    assert header["format"] == "2"
    assert int(header["size"]) == os.path.getsize(MAF_FILE)
    assert header["hash"] == COMMON.computeFileHash(MAF_FILE)


def test_hash_is_reproducible(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "rebuild")
    result = run_module(["index", MAF_FILE, out_prefix + ".block.idx", out_prefix + ".scaffold.idx"])
    assert result.returncode == 0, result.stderr

    rebuilt_header = COMMON.readIndexHeader(out_prefix + ".block.idx")
    original_header = COMMON.readIndexHeader(INDEX_FILE)
    assert rebuilt_header["hash"] == original_header["hash"]


# ---------------------------------------------------------------------------
# mafutils validate
# ---------------------------------------------------------------------------


def test_validate_clean_match_exits_0(tmp_path):
    # A freshly-built block+scaffold pair, unlike the checked-in MAF_FILE
    # fixture whose default-path scaffold index is deliberately kept
    # headerless (see HEADERLESS_SCAFFOLD_INDEX above) to exercise that path
    # elsewhere -- that fixture would (correctly) make this case
    # UNVERIFIABLE rather than VERIFIED, since it can't cross-check a
    # headerless scaffold index.
    maf_copy = os.path.join(str(tmp_path), "clean.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    result = run_module(["validate", maf_copy])
    assert result.returncode == 0, combined_output(result)
    output = combined_output(result)
    assert "VERIFIED" in output
    assert "block and scaffold index headers match" in output


def test_validate_scaffold_header_mismatch_exits_1(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "pair_mismatch.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    # Simulate the scaffold index having been built from a different file
    # state than the block index (e.g. one was rebuilt, the other wasn't).
    scaffold_index = COMMON.deriveScaffoldIndexPath(maf_copy)
    lines = open(scaffold_index).readlines()
    lines[0] = lines[0].replace("hash=md5:", "hash=md5:deadbeef")
    open(scaffold_index, "w").writelines(lines)

    result = run_module(["validate", maf_copy])
    assert result.returncode == 1, combined_output(result)
    output = combined_output(result)
    assert "MISMATCH" in output
    assert "block index" in output and "scaffold index" in output


def test_validate_missing_scaffold_index_exits_2(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "no_scaffold.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    os.remove(COMMON.deriveScaffoldIndexPath(maf_copy))

    result = run_module(["validate", maf_copy])
    assert result.returncode == 2, combined_output(result)
    output = combined_output(result)
    assert "UNVERIFIABLE" in output
    assert "Scaffold index not found" in output


def test_validate_mtime_drift_only_still_exits_0(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "touched.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    time.sleep(COMMON.MTIME_TOLERANCE_SECONDS + 1)
    os.utime(maf_copy, None)  # bump mtime without touching content

    result = run_module(["validate", maf_copy])
    assert result.returncode == 0, combined_output(result)
    output = combined_output(result)
    assert "VERIFIED" in output
    assert "mtime" in output and "differs" in output  # noted, but not fatal


def test_validate_content_change_exits_1(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "modified.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    content = open(maf_copy, "rb").read().replace(b"ACCTG", b"TTTTT", 1)
    open(maf_copy, "wb").write(content)

    result = run_module(["validate", maf_copy])
    assert result.returncode == 1, combined_output(result)
    assert "MISMATCH" in combined_output(result)


def test_validate_headerless_index_exits_2(tmp_path):
    result = run_module(["validate", MAF_FILE, HEADERLESS_SCAFFOLD_INDEX])
    assert result.returncode == 2, combined_output(result)
    assert "UNVERIFIABLE" in combined_output(result)


# ---------------------------------------------------------------------------
# --verify-hash on gc/stats/fetch
# ---------------------------------------------------------------------------


def test_gc_verify_hash_passes_on_match(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    result = run_module(["gc", MAF_FILE, INDEX_FILE, "-p", "2", "--verify-hash", "-o", out_prefix])
    assert result.returncode == 0, combined_output(result)
    assert "Verified" in combined_output(result)


def test_stats_verify_hash_fails_on_tampered_file(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "tampered.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    content = open(maf_copy, "rb").read().replace(b"ACCTG", b"TTTTT", 1)
    open(maf_copy, "wb").write(content)

    result = run_module(["stats", maf_copy, "--verify-hash", "-o", os.path.join(str(tmp_path), "out")])
    assert result.returncode != 0
    assert "content mismatch" in combined_output(result)


def test_fetch_verify_hash_errors_on_fully_headerless_index(tmp_path):
    bed_path = os.path.join(str(tmp_path), "regions.bed")
    with open(bed_path, "w") as fp:
        fp.write("chr1\t0\t10\tregion1\n")

    result = run_module([
        "fetch", MAF_FILE, bed_path,
        "--index", HEADERLESS_SCAFFOLD_INDEX,
        "--verify-hash",
        "-o", os.path.join(str(tmp_path), "out"),
    ])
    assert result.returncode != 0
    assert "no mafutils header" in combined_output(result)


def test_verify_hash_errors_when_header_present_but_no_hash_field(tmp_path):
    # Hand-crafted: a format=2-shaped header that's missing the hash= field
    # specifically (doesn't occur via normal `mafutils index` usage, which
    # always writes size/mtime/hash together, but exercises the other
    # "nothing to verify against" branch in validateIndexHeader).
    index_no_hash = os.path.join(str(tmp_path), "no_hash.block.idx")
    with open(INDEX_FILE) as src, open(index_no_hash, "w") as dst:
        lines = src.readlines()
        header = lines[0].rsplit(" hash=", 1)[0] + "\n"
        dst.write(header)
        dst.writelines(lines[1:])

    result = run_module(["gc", MAF_FILE, index_no_hash, "-p", "2", "--verify-hash", "-o", os.path.join(str(tmp_path), "out")])
    assert result.returncode != 0
    assert "no stored hash" in combined_output(result)


# ---------------------------------------------------------------------------
# Default (non-strict) validation: size mismatch errors, mtime-only warns
# ---------------------------------------------------------------------------


def test_default_validation_errors_on_size_mismatch(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "resized.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    with open(maf_copy, "a") as fp:
        fp.write("\n")  # append a byte, changing size without rebuilding the index

    result = run_module(["gc", maf_copy, "-p", "2", "-o", os.path.join(str(tmp_path), "out")])
    assert result.returncode != 0
    assert "size mismatch" in combined_output(result)


def test_default_validation_warns_only_on_mtime_mismatch(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "touched2.maf")
    shutil.copyfile(MAF_FILE, maf_copy)
    index_result = run_module(["index", maf_copy])
    assert index_result.returncode == 0, index_result.stderr

    time.sleep(COMMON.MTIME_TOLERANCE_SECONDS + 1)
    os.utime(maf_copy, None)

    result = run_module(["gc", maf_copy, "-p", "2", "-o", os.path.join(str(tmp_path), "out")])
    assert result.returncode == 0, combined_output(result)
    assert "modification time differs" in combined_output(result)
