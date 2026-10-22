import csv
import os
import re
import subprocess
import sys

import pytest

from mafutils.stats import SPECIES_SCAFFOLD_HEADERS, writeSpeciesScaffolds

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "example.maf")
TABS_MAF = os.path.join(TEST_DIR, "example-tabs.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")

# Hand-derived from tests/example.maf's s lines:
#   (species, scaffold): (src_size, blocks, aligned_bases)
EXPECTED_SCAFFOLDS = {
    ("human", "chr1"): (1000, 3, 25),
    ("human", "chr2"): (2000, 1, 6),
    ("human", "chr3"): (1000, 2, 9),
    ("human", "chr4"): (1000, 3, 14),
    ("human", "chrX"): (500, 1, 4),
    ("chimp", "chr1"): (1000, 3, 20),
    ("chimp", "chr2"): (2000, 1, 6),
    ("chimp", "chr3"): (1000, 2, 7),
    ("chimp", "chr4"): (1000, 2, 8),
    ("gorilla", "chr4"): (1000, 1, 4),
    ("gorilla", "chrX"): (800, 1, 4),
}
# species -> (scaffolds_in_maf, scaffold_bases_in_maf)
EXPECTED_ROLLUP = {"human": (5, 5500), "chimp": (4, 5000), "gorilla": (2, 1800)}

ORIGINAL_SPECIES_COLUMNS = [
    "species", "blocks_present", "blocks_missing", "presence_fraction", "copy_lines",
    "duplicated_blocks", "duplicate_copies_total", "gaps_total", "nongaps_total",
    "gap_fraction_cells", "avg_gaps_per_present_block",
    "avg_gap_fraction_per_present_block", "all_gap_blocks",
]


def run_stats(args):
    return subprocess.run(
        [sys.executable, "-m", "mafutils", "stats"] + args,
        check=False, cwd=REPO_ROOT, capture_output=True, text=True,
    )


def read_tsv(path):
    with open(path, encoding="utf-8") as fp:
        return list(csv.DictReader(fp, delimiter="\t"))


def stats_outputs(tmp_path, maf=MAF_FILE, extra=()):
    prefix = os.path.join(str(tmp_path), "out")
    args = [maf, "-o", prefix] + list(extra)
    if maf == MAF_FILE:
        args[1:1] = [INDEX_FILE]
    r = run_stats(args)
    assert r.returncode == 0, r.stdout + r.stderr
    return prefix


# ---------------------------------------------------------------------------
# species_scaffolds.tsv
# ---------------------------------------------------------------------------


def test_scaffold_table_matches_hand_derived_values(tmp_path):
    prefix = stats_outputs(tmp_path)
    rows = read_tsv(prefix + ".species_scaffolds.tsv")
    got = {(r["species"], r["scaffold"]): (int(r["src_size"]), int(r["blocks"]), int(r["aligned_bases"]))
           for r in rows}
    assert got == EXPECTED_SCAFFOLDS


def test_same_named_scaffolds_in_different_species_are_kept_apart(tmp_path):
    """
    human.chr4 and chimp.chr4 are different sequences that share a name; keying
    by scaffold alone would merge their blocks and aligned bases.
    """
    prefix = stats_outputs(tmp_path)
    rows = {(r["species"], r["scaffold"]): r for r in read_tsv(prefix + ".species_scaffolds.tsv")}
    assert rows[("human", "chr4")]["aligned_bases"] == "14"
    assert rows[("chimp", "chr4")]["aligned_bases"] == "8"
    assert rows[("gorilla", "chr4")]["aligned_bases"] == "4"


def test_scaffold_table_sorted_by_species_then_largest_first(tmp_path):
    prefix = stats_outputs(tmp_path)
    rows = read_tsv(prefix + ".species_scaffolds.tsv")
    keys = [(r["species"], -int(r["src_size"]), r["scaffold"]) for r in rows]
    assert keys == sorted(keys)
    human = [r["scaffold"] for r in rows if r["species"] == "human"]
    assert human[0] == "chr2"            # 2000, the largest


def test_aligned_fraction_is_aligned_over_src_size(tmp_path):
    prefix = stats_outputs(tmp_path)
    for r in read_tsv(prefix + ".species_scaffolds.tsv"):
        assert float(r["aligned_fraction"]) == pytest.approx(int(r["aligned_bases"]) / int(r["src_size"]))


def test_aligned_fraction_is_not_clamped(tmp_path):
    """A region aligned into several blocks can exceed its scaffold length;
    that is information, so it must not be capped at 1.0."""
    out = os.path.join(str(tmp_path), "s.tsv")
    writeSpeciesScaffolds(out, {("sp", "s1"): {"src_size": 10, "blocks": 3, "aligned_bases": 25}})
    row = read_tsv(out)[0]
    assert list(row) == SPECIES_SCAFFOLD_HEADERS
    assert row["aligned_fraction"] == "2.50000000"


def test_missing_scaffold_name_and_size_are_written_as_dots(tmp_path):
    out = os.path.join(str(tmp_path), "s.tsv")
    writeSpeciesScaffolds(out, {("sp", ""): {"src_size": None, "blocks": 1, "aligned_bases": 4}})
    row = read_tsv(out)[0]
    assert (row["scaffold"], row["src_size"], row["aligned_fraction"]) == (".", ".", ".")


def test_tab_delimited_maf_gives_identical_scaffold_table(tmp_path):
    os.makedirs(tmp_path / "space")
    os.makedirs(tmp_path / "tabs")
    space = stats_outputs(tmp_path / "space")
    tabs = stats_outputs(tmp_path / "tabs", maf=TABS_MAF)
    assert read_tsv(space + ".species_scaffolds.tsv") == read_tsv(tabs + ".species_scaffolds.tsv")


# ---------------------------------------------------------------------------
# species.tsv roll-up columns
# ---------------------------------------------------------------------------


def test_species_rollup_columns(tmp_path):
    prefix = stats_outputs(tmp_path)
    rows = {r["species"]: r for r in read_tsv(prefix + ".species.tsv")}
    for sp, (n, bases) in EXPECTED_ROLLUP.items():
        assert int(rows[sp]["scaffolds_in_maf"]) == n, sp
        assert int(rows[sp]["scaffold_bases_in_maf"]) == bases, sp


def test_new_species_columns_are_appended_not_inserted(tmp_path):
    """Positional parsers of species.tsv must keep working."""
    prefix = stats_outputs(tmp_path)
    header = open(prefix + ".species.tsv").readline().rstrip("\n").split("\t")
    assert header[:len(ORIGINAL_SPECIES_COLUMNS)] == ORIGINAL_SPECIES_COLUMNS
    assert header[len(ORIGINAL_SPECIES_COLUMNS):] == ["scaffolds_in_maf", "scaffold_bases_in_maf"]


def test_expected_but_absent_species_gets_zero_scaffolds(tmp_path):
    species_file = os.path.join(str(tmp_path), "sp.txt")
    open(species_file, "w").write("human\nchimp\ngorilla\norangutan\n")
    prefix = stats_outputs(tmp_path, extra=["--expected-species-file", species_file])
    rows = {r["species"]: r for r in read_tsv(prefix + ".species.tsv")}
    assert rows["orangutan"]["scaffolds_in_maf"] == "0"
    assert rows["orangutan"]["scaffold_bases_in_maf"] == "0"


# ---------------------------------------------------------------------------
# Invariants against existing, independently pinned statistics
# ---------------------------------------------------------------------------


def test_per_scaffold_aligned_bases_sum_to_nongaps_total(tmp_path):
    """
    Every non-gap base of a species comes from exactly one of its scaffolds, so
    the per-scaffold sums must equal the species' existing nongaps_total --
    an independent check on the new accounting.
    """
    prefix = stats_outputs(tmp_path)
    species = {r["species"]: r for r in read_tsv(prefix + ".species.tsv")}
    totals = {}
    for r in read_tsv(prefix + ".species_scaffolds.tsv"):
        totals[r["species"]] = totals.get(r["species"], 0) + int(r["aligned_bases"])
    for sp, total in totals.items():
        assert total == int(species[sp]["nongaps_total"]), sp


def test_per_scaffold_blocks_cover_blocks_present(tmp_path):
    """
    sum(blocks) over a species' scaffolds is >= its blocks_present, with
    equality when no block carries two copies of the species (true here).
    """
    prefix = stats_outputs(tmp_path)
    species = {r["species"]: r for r in read_tsv(prefix + ".species.tsv")}
    totals = {}
    for r in read_tsv(prefix + ".species_scaffolds.tsv"):
        totals[r["species"]] = totals.get(r["species"], 0) + int(r["blocks"])
    for sp, total in totals.items():
        assert total == int(species[sp]["blocks_present"]), sp


@pytest.mark.parametrize("processes,chunk_size", [(1, 1), (1, 3), (2, 1), (2, 2), (4, 3)])
def test_scaffold_outputs_invariant_to_chunking(processes, chunk_size, tmp_path):
    os.makedirs(tmp_path / "base")
    os.makedirs(tmp_path / "run")
    base = stats_outputs(tmp_path / "base", extra=["-p", "1", "--chunk-size", "10000"])
    run = stats_outputs(tmp_path / "run", extra=["-p", str(processes), "--chunk-size", str(chunk_size)])
    assert read_tsv(run + ".species_scaffolds.tsv") == read_tsv(base + ".species_scaffolds.tsv")
    cols = ["species", "scaffolds_in_maf", "scaffold_bases_in_maf"]
    pick = lambda path: [{c: r[c] for c in cols} for r in read_tsv(path)]
    assert pick(run + ".species.tsv") == pick(base + ".species.tsv")


# ---------------------------------------------------------------------------
# srcSize integrity
# ---------------------------------------------------------------------------


def _altered_copy(tmp_path, old, new):
    """A copy of example.maf with one s line changed, freshly indexed."""
    text = open(MAF_FILE).read()
    assert text.count(old) == 1
    maf = os.path.join(str(tmp_path), "altered.maf")
    open(maf, "w").write(text.replace(old, new))
    r = subprocess.run([sys.executable, "-m", "mafutils", "index", maf],
                       check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    return maf


@pytest.mark.parametrize("processes,chunk_size", [(1, 10000), (2, 3)])
def test_conflicting_src_size_fails_cleanly(processes, chunk_size, tmp_path):
    """
    One species.scaffold reporting two lengths is an inconsistent MAF. It must
    fail -- both within a chunk and when the two lines land in different
    chunks merged across workers -- naming the scaffold, without a traceback.
    """
    maf = _altered_copy(tmp_path, "s chimp.chr4 28 5 + 1000 TGG-C-T", "s chimp.chr4 28 5 + 998 TGG-C-T")
    r = run_stats([maf, "-o", os.path.join(str(tmp_path), "out"),
                   "-p", str(processes), "--chunk-size", str(chunk_size)])
    combined = r.stdout + r.stderr
    assert r.returncode != 0
    assert re.search(r"chimp\.chr4 has conflicting srcSize values \((1000 and 998|998 and 1000)\)", combined), combined
    assert "Traceback" not in combined


def test_unparseable_src_size_affects_only_the_scaffold_table(tmp_path):
    maf = _altered_copy(tmp_path, "s chimp.chr4 28 5 + 1000 TGG-C-T", "s chimp.chr4 28 5 + abc TGG-C-T")
    prefix = os.path.join(str(tmp_path), "out")
    r = run_stats([maf, "-o", prefix])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "did not parse as an integer" in (r.stdout + r.stderr)

    overall = {row["metric"]: row["value"] for row in read_tsv(prefix + ".overall.tsv")}
    assert overall["scaffold_size_parse_errors"] == "1"
    # the line still counts toward existing statistics, unchanged
    assert overall["total_sequence_lines"] == "20"
    species = {r["species"]: r for r in read_tsv(prefix + ".species.tsv")}
    assert species["chimp"]["nongaps_total"] == "41"
    # the other chimp.chr4 line (block 8) still supplies its size
    rows = {(r["species"], r["scaffold"]): r for r in read_tsv(prefix + ".species_scaffolds.tsv")}
    assert rows[("chimp", "chr4")]["src_size"] == "1000"


# ---------------------------------------------------------------------------
# Dashboard
# ---------------------------------------------------------------------------


def test_dashboard_states_truncation_and_shows_largest_first(tmp_path):
    pytest.importorskip("matplotlib")
    prefix = stats_outputs(tmp_path, extra=["--html-dashboard", "--dashboard-top-scaffolds", "2"])
    page = open(prefix + ".dashboard.html").read()
    summaries = [re.sub(r"<[^>]+>", "", m) for m in re.findall(r"<summary>.*?</summary>", page)]
    assert any("human" in s and "showing 2 of 5 scaffolds (largest by length)" in s for s in summaries), summaries
    assert any("gorilla" in s and "showing 2 of 2 scaffolds" in s and "largest" not in s for s in summaries)
    human = re.search(r"<summary><b>human</b>.*?</details>", page, re.S).group(0)
    first = re.search(r"<tbody><tr><td>([^<]+)</td>", human).group(1)
    assert first == "chr2"
    assert "present in the" in page and "MAF</b>" in page
