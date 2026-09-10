import os
import subprocess
import sys

import pytest

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")
SCAFFOLD_INDEX = os.path.join(TEST_DIR, "example.maf.scaffold.idx")
BED_FILE = os.path.join(TEST_DIR, "example.bed")

# 8 real blocks on one scaffold (CM000994.3), extracted from production data.
REAL_MAF = os.path.join(TEST_DIR, "real-excerpt.maf")
REAL_INDEX = os.path.join(TEST_DIR, "real-excerpt.maf.block.idx")


def run_fetch(args):
    return subprocess.run(
        [sys.executable, "-m", "mafutils", "fetch"] + args,
        check=False, cwd=REPO_ROOT, capture_output=True, text=True,
    )


def outputs(directory, ext=".maf"):
    if not os.path.isdir(directory):
        return {}
    return {
        name: open(os.path.join(directory, name), "rb").read()
        for name in os.listdir(directory)
        if name.endswith(ext)
    }


# ---------------------------------------------------------------------------
# Block mode streams the index, so results must not depend on how regions are
# partitioned across workers, nor on the order they appear in the BED.
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("processes", [1, 2, 3, 8])
def test_output_is_invariant_to_process_count(processes, tmp_path):
    baseline = os.path.join(str(tmp_path), "p1")
    assert run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", baseline, "-p", "1"]).returncode == 0

    out = os.path.join(str(tmp_path), f"p{processes}")
    result = run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", out, "-p", str(processes)])
    assert result.returncode == 0, result.stderr
    assert outputs(out) == outputs(baseline)


def test_output_is_invariant_to_bed_order(tmp_path):
    """
    Regions are sorted into index order internally, so a shuffled BED must
    produce byte-identical per-region output. (Guards the sort/merge-join: if
    the streaming assumed BED order, a reversed BED would seek backwards and
    silently drop blocks.)
    """
    with open(BED_FILE) as fp:
        lines = [ln for ln in fp if ln.strip()]

    forward = os.path.join(str(tmp_path), "forward")
    assert run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", forward]).returncode == 0

    reversed_bed = os.path.join(str(tmp_path), "reversed.bed")
    with open(reversed_bed, "w") as fp:
        # Normalize line endings: example.bed's last line has no trailing
        # newline, so writing the reversed list verbatim would splice two
        # records together.
        fp.write("".join(ln.rstrip("\n") + "\n" for ln in reversed(lines)))

    out = os.path.join(str(tmp_path), "reversed")
    result = run_fetch([MAF_FILE, reversed_bed, "--index", INDEX_FILE, "-b", "id", "-o", out])
    assert result.returncode == 0, result.stderr
    assert outputs(out) == outputs(forward)


def test_scaffold_appearing_in_multiple_runs_is_fully_covered(tmp_path):
    """
    example.maf interleaves scaffolds (chr4 ... chrX ... chr4), so chr4 occupies
    two separate runs of the index. A streaming implementation that assumes one
    contiguous run per scaffold silently loses the later blocks -- 08-span-
    multiple-gaps lives in chr4's SECOND run, so its output is the canary.
    """
    with open(SCAFFOLD_INDEX) as fp:
        scaffolds = [ln.split("\t")[0] for ln in fp if not ln.startswith("#") and ln.strip()]
    assert len(scaffolds) > len(set(scaffolds)), "fixture no longer exercises multi-run scaffolds"

    bed = os.path.join(str(tmp_path), "chr4.bed")
    with open(bed, "w") as fp:
        fp.write("chr4\t4\t8\tfirst-run\n")     # chr4's first run
        fp.write("chr4\t31\t33\tsecond-run\n")  # chr4's second run, after chrX

    out = os.path.join(str(tmp_path), "out")
    result = run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-b", "id", "-o", out])
    assert result.returncode == 0, result.stderr

    for name in ("first-run.maf", "second-run.maf"):
        path = os.path.join(out, name)
        assert os.path.isfile(path) and os.path.getsize(path) > 0, f"{name} missing/empty"
        assert "chr4" in open(path).read()


# ---------------------------------------------------------------------------
# Single large scaffold: the degenerate shape for a streaming design.
# ---------------------------------------------------------------------------


def test_many_regions_on_one_scaffold_use_all_workers(tmp_path):
    """All 8 blocks of real-excerpt.maf are on one scaffold; splitting by region
    count (not by scaffold) must still spread them across workers."""
    bed = os.path.join(str(tmp_path), "many.bed")
    with open(bed, "w") as fp:
        for i in range(8):
            start = 125888838 + i * 100
            fp.write(f"CM000994.3\t{start}\t{start + 50}\tr{i}\n")

    single = os.path.join(str(tmp_path), "p1")
    assert run_fetch([REAL_MAF, bed, "--index", REAL_INDEX, "-b", "id", "-o", single, "-p", "1"]).returncode == 0

    out = os.path.join(str(tmp_path), "p4")
    result = run_fetch([REAL_MAF, bed, "--index", REAL_INDEX, "-b", "id", "-o", out, "-p", "4"])
    assert result.returncode == 0, result.stderr
    assert "into 4 batches" in (result.stdout + result.stderr)
    assert outputs(out) == outputs(single)


def test_single_region_spanning_whole_scaffold(tmp_path):
    """
    One BED region covering an entire scaffold. The blocks for it are streamed
    rather than gathered into a list, so this must work without holding every
    overlapping block in memory.
    """
    bed = os.path.join(str(tmp_path), "whole.bed")
    with open(bed, "w") as fp:
        fp.write("CM000994.3\t0\t250000000\twhole\n")

    out = os.path.join(str(tmp_path), "out")
    result = run_fetch([REAL_MAF, bed, "--index", REAL_INDEX, "-b", "id", "-o", out])
    assert result.returncode == 0, result.stderr

    content = open(os.path.join(out, "whole.maf")).read()
    # All 8 blocks of the excerpt fall inside this range.
    assert content.count("\na") + content.startswith("a") == 8, "expected every block in the scaffold"


# ---------------------------------------------------------------------------
# The block/scaffold index pair must be consistent, and say so if not.
# ---------------------------------------------------------------------------


def test_mismatched_index_pair_fails_loudly(tmp_path):
    """
    Block mode locates blocks via byte ranges recorded in the scaffold index, so
    a pair from different `mafutils index` runs cannot be trusted together. An
    off-by-one boundary alone is enough to skip every run's first block, which
    would silently drop regions -- so it must be an error, not a warning.
    """
    maf = os.path.join(str(tmp_path), "t.maf")
    open(maf, "wb").write(open(MAF_FILE, "rb").read())
    open(maf + ".block.idx", "wb").write(open(INDEX_FILE, "rb").read())

    # Scaffold index whose header disagrees with the block index's.
    with open(SCAFFOLD_INDEX) as src, open(maf + ".scaffold.idx", "w") as dst:
        for line in src:
            dst.write(line.replace("hash=md5:", "hash=md5:deadbeef") if line.startswith("#") else line)

    bed = os.path.join(str(tmp_path), "r.bed")
    with open(bed, "w") as fp:
        fp.write("chr1\t2\t6\tr1\n")

    result = run_fetch([maf, bed, "--index", maf + ".block.idx", "-b", "id",
                        "-o", os.path.join(str(tmp_path), "out")])
    assert result.returncode != 0
    assert "do not come from the same" in (result.stdout + result.stderr)
