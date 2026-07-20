import os
import subprocess
import sys

import pytest

from mafutils.lib import common as COMMON

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_NONE = os.path.join(TEST_DIR, "example.maf")
INDEX_NONE = os.path.join(TEST_DIR, "example.maf.block.idx")

MAF_GZ = os.path.join(TEST_DIR, "example.maf.gz")
INDEX_GZ = os.path.join(TEST_DIR, "example.maf.gz.block.idx")

MAF_BGZ = os.path.join(TEST_DIR, "example.maf.bgz")
INDEX_BGZ = os.path.join(TEST_DIR, "example.maf.bgz.block.idx")

BED_FILE = os.path.join(TEST_DIR, "example.bed")


def run_module(args, cwd=REPO_ROOT):
    cmd = [sys.executable, "-m", "mafutils"] + args
    return subprocess.run(cmd, check=False, cwd=cwd, capture_output=True, text=True)


def combined_output(result):
    # loginit.py routes INFO to stdout and WARNING/ERROR to stderr -- check
    # both together so these tests don't depend on exactly which level a
    # given log line is at.
    return result.stdout + result.stderr


# ---------------------------------------------------------------------------
# detectCompression / index header
# ---------------------------------------------------------------------------


def test_detect_compression_distinguishes_bgzip_from_gzip():
    assert COMMON.detectCompression(MAF_NONE) == "none"
    assert COMMON.detectCompression(MAF_GZ) == "gz"
    assert COMMON.detectCompression(MAF_BGZ) == "bgzip"


def test_index_files_carry_a_matching_header():
    for maf, index, expected_compression in [
        (MAF_NONE, INDEX_NONE, "none"),
        (MAF_GZ, INDEX_GZ, "gz"),
        (MAF_BGZ, INDEX_BGZ, "bgzip"),
    ]:
        header = COMMON.readIndexHeader(index)
        assert header is not None
        assert header["compression"] == expected_compression
        assert header["maf"] == os.path.basename(maf)


# ---------------------------------------------------------------------------
# gc: bgzip parallel path must match the uncompressed/sequential baseline
# ---------------------------------------------------------------------------


def test_gc_bgzip_parallel_matches_uncompressed(tmp_path):
    baseline_prefix = os.path.join(str(tmp_path), "baseline")
    bgz_prefix = os.path.join(str(tmp_path), "bgz")

    result = run_module(["gc", MAF_NONE, "-o", baseline_prefix])
    assert result.returncode == 0, result.stderr

    result = run_module(["gc", MAF_BGZ, INDEX_BGZ, "-o", bgz_prefix, "-p", "2"])
    assert result.returncode == 0, result.stderr
    assert "PARALLEL" in combined_output(result)

    baseline_csv = open(baseline_prefix + ".gc.csv").read()
    bgz_csv = open(bgz_prefix + ".gc.csv").read()
    assert baseline_csv == bgz_csv


# ---------------------------------------------------------------------------
# stats: bgzip parallel path must match; gzip must force a single process
# ---------------------------------------------------------------------------


def test_stats_bgzip_parallel_matches_uncompressed(tmp_path):
    baseline_prefix = os.path.join(str(tmp_path), "baseline")
    bgz_prefix = os.path.join(str(tmp_path), "bgz")

    result = run_module(["stats", MAF_NONE, INDEX_NONE, "-o", baseline_prefix])
    assert result.returncode == 0, result.stderr

    result = run_module(["stats", MAF_BGZ, INDEX_BGZ, "-o", bgz_prefix, "-p", "3", "--chunk-size", "2"])
    assert result.returncode == 0, result.stderr
    assert "tasks" in combined_output(result)

    assert open(baseline_prefix + ".overall.tsv").read() == open(bgz_prefix + ".overall.tsv").read()
    assert open(baseline_prefix + ".species.tsv").read() == open(bgz_prefix + ".species.tsv").read()


def test_stats_gzip_forces_single_process_and_matches(tmp_path):
    baseline_prefix = os.path.join(str(tmp_path), "baseline")
    gz_prefix = os.path.join(str(tmp_path), "gz")

    result = run_module(["stats", MAF_NONE, INDEX_NONE, "-o", baseline_prefix])
    assert result.returncode == 0, result.stderr

    result = run_module(["stats", MAF_GZ, INDEX_GZ, "-o", gz_prefix, "-p", "4"])
    assert result.returncode == 0, result.stderr
    output = combined_output(result)
    assert "inefficient on" in output
    assert "across 1 process" in output

    assert open(baseline_prefix + ".overall.tsv").read() == open(gz_prefix + ".overall.tsv").read()


# ---------------------------------------------------------------------------
# fetch: bgzip parallel path, gzip fallback, and the nested/overlapping-region
# case that motivated the gzip redesign (a single sort key can't order this
# correctly; the prefetch-cache approach sidesteps the problem entirely).
# ---------------------------------------------------------------------------


def _normalize_source_line(text):
    return "\n".join(
        line if not line.startswith("## Source MAF:") else "## Source MAF: X"
        for line in text.splitlines()
    )


def test_fetch_bgzip_parallel_matches_uncompressed(tmp_path):
    baseline_dir = os.path.join(str(tmp_path), "baseline")
    bgz_dir = os.path.join(str(tmp_path), "bgz")

    result = run_module(["fetch", MAF_NONE, BED_FILE, "--index", INDEX_NONE, "-b", "id", "-o", baseline_dir])
    assert result.returncode == 0, result.stderr

    result = run_module(["fetch", MAF_BGZ, BED_FILE, "--index", INDEX_BGZ, "-b", "id", "-o", bgz_dir, "-p", "4"])
    assert result.returncode == 0, result.stderr

    for name in os.listdir(baseline_dir):
        if not name.endswith(".maf"):
            continue
        baseline_text = _normalize_source_line(open(os.path.join(baseline_dir, name)).read())
        bgz_text = _normalize_source_line(open(os.path.join(bgz_dir, name)).read())
        assert baseline_text == bgz_text, f"mismatch for {name}"


def test_fetch_gzip_falls_back_and_matches(tmp_path):
    baseline_dir = os.path.join(str(tmp_path), "baseline")
    gz_dir = os.path.join(str(tmp_path), "gz")

    result = run_module(["fetch", MAF_NONE, BED_FILE, "--index", INDEX_NONE, "-b", "id", "-o", baseline_dir])
    assert result.returncode == 0, result.stderr

    result = run_module(["fetch", MAF_GZ, BED_FILE, "--index", INDEX_GZ, "-b", "id", "-o", gz_dir, "-p", "4"])
    assert result.returncode == 0, result.stderr
    output = combined_output(result)
    assert "inefficient on" in output
    assert "Prefetching" in output

    for name in os.listdir(baseline_dir):
        if not name.endswith(".maf"):
            continue
        baseline_text = _normalize_source_line(open(os.path.join(baseline_dir, name)).read())
        gz_text = _normalize_source_line(open(os.path.join(gz_dir, name)).read())
        assert baseline_text == gz_text, f"mismatch for {name}"


def test_fetch_gzip_nested_overlapping_regions(tmp_path):
    """
    Regression test for the gzip redesign: one region spans multiple blocks,
    another is fully nested inside that span. No single ascending sort key
    over regions can order both correctly (confirmed during design review),
    which is why fetch precomputes+dedupes the needed blocks globally instead
    of relying on region processing order at all.

    "outer-spanning" (chr1:0-40) is itself a pre-existing no-output case
    (a trim-mismatch on a reference gap, identical in both compression
    paths -- not something this test is checking) -- padded with guaranteed
    single-block hits so that one no-output region doesn't trip fetch's
    unrelated zero-overlap-fraction fail-fast threshold.
    """
    bed_path = os.path.join(str(tmp_path), "nested.bed")
    with open(bed_path, "w") as fp:
        fp.write("chr1\t18\t26\tinner-nested\n")
        fp.write("chr1\t0\t40\touter-spanning\n")
        for i in range(10):
            fp.write(f"chr2\t5\t10\tpadding-{i}\n")

    baseline_dir = os.path.join(str(tmp_path), "baseline")
    gz_dir = os.path.join(str(tmp_path), "gz")

    result = run_module(["fetch", MAF_NONE, bed_path, "--index", INDEX_NONE, "-b", "id", "-o", baseline_dir])
    assert result.returncode == 0, result.stderr

    result = run_module(["fetch", MAF_GZ, bed_path, "--index", INDEX_GZ, "-b", "id", "-o", gz_dir])
    assert result.returncode == 0, result.stderr

    for name in ["inner-nested.maf", "outer-spanning.maf"]:
        baseline_path = os.path.join(baseline_dir, name)
        gz_path = os.path.join(gz_dir, name)
        baseline_exists = os.path.exists(baseline_path)
        assert baseline_exists == os.path.exists(gz_path), f"existence mismatch for {name}"
        if baseline_exists:
            baseline_text = _normalize_source_line(open(baseline_path).read())
            gz_text = _normalize_source_line(open(gz_path).read())
            assert baseline_text == gz_text, f"mismatch for {name}"


# ---------------------------------------------------------------------------
# Index/MAF compression mismatch must error; a missing header must only warn.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("command", ["gc", "stats"])
def test_mismatched_index_compression_errors(command, tmp_path):
    args = [command, MAF_NONE, INDEX_BGZ, "-o", os.path.join(str(tmp_path), "out")]
    if command == "gc":
        args = args + ["-p", "2"]
    result = run_module(args)
    assert result.returncode != 0
    assert "compression mismatch" in combined_output(result)


def test_missing_header_warns_not_errors(tmp_path):
    """A block index with no mafutils header (as if built by an older
    version) must still work, with a warning rather than a hard failure."""
    headerless_index = os.path.join(str(tmp_path), "headerless.block.idx")
    with open(INDEX_NONE) as src, open(headerless_index, "w") as dst:
        lines = src.readlines()
        assert lines[0].startswith("#"), "expected a header line to strip"
        dst.writelines(lines[1:])
    assert COMMON.readIndexHeader(headerless_index) is None

    result = run_module(["gc", MAF_NONE, headerless_index, "-o", os.path.join(str(tmp_path), "out"), "-p", "2"])
    assert result.returncode == 0, result.stderr
    assert "no mafutils header" in combined_output(result)


# ---------------------------------------------------------------------------
# Regression: an 'a' line straddling a real BGZF block boundary
# ---------------------------------------------------------------------------


def test_itermafblocks_handles_bgzf_block_boundary_straddling_a_line(tmp_path):
    """
    Regression test for a real bug found on a ~114k-BGZF-block production
    file: iterMafBlocks() used to compute block boundaries via
    `stream.tell() - len(line)`, which is invalid arithmetic on BGZF's
    packed virtual offsets whenever the 'a' line straddles a real BGZF
    block boundary -- it can "borrow" across the packed representation,
    landing on a byte position that was never a valid block start. This
    engineers that exact scenario using the vendored BgzfWriter's
    known-precise 65536-byte block flushing, rather than relying on a
    multi-gigabyte fixture to hit it by chance.
    """
    from mafutils.lib import bgzf

    content = "##maf version=1\n"
    i = 0
    # Coarse approach: bulk filler blocks to get within ~200 bytes of a
    # 65536-byte (one BGZF block) boundary.
    while (65536 - (len(content) % 65536)) > 200:
        content += f"a score=1\ns human.chr1 {i * 4} 4 + 1000 ACGT\n\n"
        i += 1

    straddle_line = "a score=999\n"
    # Fine adjustment: pad with a comment line (skipped by iterMafBlocks,
    # but still consumes bytes/position) to land with a small gap -- one
    # smaller than the straddle line -- before the boundary.
    target_gap = 6
    gap = 65536 - (len(content) % 65536)
    pad_len = gap - target_gap
    assert pad_len > 0, "coarse filler overshot the target; adjust filler size"
    content += "#" * (pad_len - 1) + "\n"
    assert 65536 - (len(content) % 65536) == target_gap

    straddle_start = len(content)
    content += straddle_line
    straddle_end = len(content)
    assert straddle_start // 65536 != (straddle_end - 1) // 65536, (
        "test construction failed to straddle a 65536-byte boundary"
    )
    content += f"s human.chr1 {i * 4} 4 + 1000 TTTT\n\n"
    for _ in range(5):
        i += 1
        content += f"a score=1\ns human.chr1 {i * 4} 4 + 1000 ACGT\n\n"

    maf_bytes = content.encode("ascii")

    bgz_path = os.path.join(str(tmp_path), "boundary.maf.bgz")
    with bgzf.BgzfWriter(bgz_path, "wb") as f:
        f.write(maf_bytes)

    none_path = os.path.join(str(tmp_path), "boundary.maf")
    with open(none_path, "wb") as f:
        f.write(maf_bytes)

    # Ground truth: confirm a real BGZF block boundary falls inside the
    # straddle line (i.e. the test actually engineered the scenario it
    # claims to).
    real_block_data_starts = []
    with open(bgz_path, "rb") as f:
        for _start_offset, _block_length, data_start, _data_len in bgzf.BgzfBlocks(f):
            real_block_data_starts.append(data_start)
    assert any(straddle_start < b < straddle_end for b in real_block_data_starts), (
        "test construction failed: no real BGZF block boundary falls inside the straddle line"
    )

    # The actual regression check: offsets computed against the bgzip copy
    # must match those computed against the uncompressed copy in count, and
    # every block must be independently readable via a real seek.
    with COMMON.openMaf(none_path, "none", "rt") as f:
        none_offsets = [(s, e) for _, s, e in COMMON.iterMafBlocks(f)]

    with COMMON.openMaf(bgz_path, "bgzip", "rt") as f:
        bgz_offsets = [(s, e) for _, s, e in COMMON.iterMafBlocks(f)]

    assert len(none_offsets) == len(bgz_offsets) == i + 1

    with COMMON.openMaf(bgz_path, "bgzip", "rb") as maf_fp:
        for start, end in bgz_offsets:
            block_bytes = COMMON.readMafBlockBytes(maf_fp, "bgzip", start, end)
            assert block_bytes.startswith(b"a"), f"unreadable/misaligned block at {start}-{end}"


# ---------------------------------------------------------------------------
# Index path auto-derivation
# ---------------------------------------------------------------------------


def test_index_auto_derivation_for_gc(tmp_path):
    maf_copy = os.path.join(str(tmp_path), "auto.maf")
    with open(MAF_NONE) as src, open(maf_copy, "w") as dst:
        dst.write(src.read())

    result = run_module(["index", maf_copy])
    assert result.returncode == 0, result.stderr
    assert os.path.isfile(maf_copy + ".block.idx")
    assert os.path.isfile(maf_copy + ".scaffold.idx")

    result = run_module(["gc", maf_copy, "-o", os.path.join(str(tmp_path), "out"), "-p", "2"])
    assert result.returncode == 0, result.stderr
    assert "Found index at default location" in combined_output(result)
