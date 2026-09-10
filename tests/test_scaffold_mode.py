import io
import os
import subprocess
import sys

import pytest

from mafutils.lib import common as COMMON

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "example.maf")
SCAFFOLD_INDEX = os.path.join(TEST_DIR, "example.maf.scaffold.idx")

HEADER = "# mafutils-index format=2 maf=t.maf compression=none size=776 mtime=1.0 hash=md5:x\n"


def run_fetch(args):
    return subprocess.run(
        [sys.executable, "-m", "mafutils", "fetch"] + args,
        check=False, cwd=REPO_ROOT, capture_output=True, text=True,
    )


def write_bed(path, lines):
    with open(path, "w") as fp:
        for line in lines:
            fp.write(line + "\n")


# ---------------------------------------------------------------------------
# copyMafRangeToStream: the streaming replacement for a single giant read
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("compression,fixture", [
    ("none", "example.maf"),
    ("gz", "example.maf.gz"),
    ("bgzip", "example.maf.bgz"),
])
@pytest.mark.parametrize("chunk_size", [1, 7, 64, 8 * 1024 * 1024])
def test_copy_range_matches_read_bytes_at_any_chunk_size(compression, fixture, chunk_size):
    """
    Streaming a range in bounded chunks must produce exactly what the
    single-read version produces, including when the chunk size is far
    smaller than the range (so the copy loop runs many times) -- the real
    fixtures are small enough that a default chunk size would only ever
    loop once, which would leave the multi-chunk path untested.
    """
    maf = os.path.join(TEST_DIR, fixture)
    index = maf + ".block.idx"

    with open(index) as fp:
        rows = [ln.split("\t") for ln in fp if ln.strip() and not ln.startswith("#")]
    start, end = int(rows[0][6]), int(rows[-1][7])  # span many blocks at once

    with COMMON.openMaf(maf, compression, "rb") as handle:
        expected = COMMON.readMafBlockBytes(handle, compression, start, end)

    buf = io.BytesIO()
    with COMMON.openMaf(maf, compression, "rb") as handle:
        written = COMMON.copyMafRangeToStream(handle, compression, start, end, buf, chunk_size=chunk_size)

    assert buf.getvalue() == expected
    assert written == len(expected)


# ---------------------------------------------------------------------------
# fetch -m scaffold: must fail loudly, never report success on failure
# ---------------------------------------------------------------------------


def test_scaffold_mode_succeeds_normally(tmp_path):
    bed = os.path.join(str(tmp_path), "s.bed")
    write_bed(bed, ["chr1\t0\t50\tchr1", "chr4\t0\t50\tchr4"])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([MAF_FILE, bed, "-i", SCAFFOLD_INDEX, "-m", "scaffold", "-o", out])
    assert result.returncode == 0, result.stderr

    for scaffold in ("chr1", "chr4"):
        path = os.path.join(out, f"{scaffold}.maf")
        assert os.path.getsize(path) > 0, f"{scaffold} output is empty"
        with open(path) as fp:
            assert fp.read().lstrip().startswith("#") or "a" in fp.read(), path


def test_scaffold_mode_exits_nonzero_when_extraction_yields_nothing(tmp_path):
    """
    Regression test for a silent-failure bug: a scaffold that produced no
    output was logged as an error but still reported as "Wrote <scaffold>",
    and fetch exited 0 -- so a downstream pipeline consumed a 0-byte file
    believing it had succeeded. Simulated here with a stale index whose
    offsets point past EOF; the real trigger was a MemoryError from trying
    to read a 597GB scaffold in one allocation.
    """
    maf = os.path.join(str(tmp_path), "t.maf")
    with open(MAF_FILE) as src, open(maf, "w") as dst:
        dst.write(src.read())

    stale_index = os.path.join(str(tmp_path), "stale.scaffold.idx")
    with open(stale_index, "w") as fp:
        fp.write(HEADER)
        fp.write("chr1\t900000\t999999\n")

    bed = os.path.join(str(tmp_path), "s.bed")
    write_bed(bed, ["chr1\t0\t50\tchr1"])

    result = run_fetch([maf, bed, "-i", stale_index, "-m", "scaffold",
                        "-o", os.path.join(str(tmp_path), "out")])
    combined = result.stdout + result.stderr

    assert result.returncode != 0, "a scaffold that extracted nothing must fail the run"
    assert "Wrote chr1" not in combined, "must not claim it wrote a scaffold it failed to extract"
    assert "0 alignment bytes" in combined
