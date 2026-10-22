import os
import subprocess
import sys

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))
MAF_FILE = os.path.join(TEST_DIR, "example.maf")
BED_FILE = os.path.join(TEST_DIR, "example.bed")


def run_fetch(args):
    return subprocess.run([sys.executable, "-m", "mafutils", "fetch"] + args,
                          check=False, cwd=REPO_ROOT, capture_output=True, text=True)


# ---------------------------------------------------------------------------
# maf_fetch_summary.tsv: batch order, not worker-completion order
# ---------------------------------------------------------------------------


def read_summary(out_dir):
    with open(os.path.join(out_dir, "maf_fetch_summary.tsv")) as fp:
        return fp.read().splitlines()


def test_summary_rows_identical_for_any_process_count(tmp_path):
    """
    Rows used to be written as each batch finished, so their order followed
    worker timing. They are now concatenated in batch order, and batches
    partition the index-ordered regions, so the file is the same however many
    workers there are (and however they are scheduled).
    """
    tables = {}
    for p in (1, 2, 3, 4):
        out = os.path.join(str(tmp_path), f"p{p}")
        r = run_fetch([MAF_FILE, BED_FILE, "-b", "id", "-o", out, "-p", str(p)])
        assert r.returncode == 0, r.stderr
        tables[p] = read_summary(out)
    assert tables[1] == tables[2] == tables[3] == tables[4]
    assert len(tables[1]) == 1 + 10                                  # header + every region


def test_summary_rows_follow_index_order(tmp_path):
    """Regions are sorted into index order (MAF scaffold order, then start)."""
    r = run_fetch([MAF_FILE, BED_FILE, "-b", "id", "-o", str(tmp_path), "-p", "3"])
    assert r.returncode == 0, r.stderr
    with open(MAF_FILE + ".scaffold.idx") as fp:
        first_run = {}
        for ln in fp:
            if ln.startswith("#") or not ln.strip():
                continue
            first_run.setdefault(ln.split("\t")[0], len(first_run))
    rows = [ln.split("\t") for ln in read_summary(str(tmp_path))[1:]]
    keys = [(first_run.get(r[0], len(first_run)), int(r[1]), int(r[2])) for r in rows]
    assert keys == sorted(keys)
