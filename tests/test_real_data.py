import csv
import os
import subprocess
import sys

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "real-excerpt.maf")
INDEX_FILE = os.path.join(TEST_DIR, "real-excerpt.maf.block.idx")

# 8 real, complete alignment blocks, extracted read-only from gwct's real
# ~42GB production MAF (data/hamsters/uncompressed/...), a couple GB in to
# skip past the initial all-N region. Unlike example.maf (hand-crafted, used
# to independently verify exact output values in test_stats.py/test_gc.py),
# this fixture exists to catch real-world-data surprises a hand-crafted
# fixture wouldn't think to include (real species-naming conventions, real
# gap patterns, real block-size distribution) -- so these tests check
# plausibility/robustness/"does it run cleanly on real data", not exact
# hand-computed values.
EXPECTED_SPECIES = {
    "Phodopus_sungorus_new",
    "Phodopus_campbelli",
    "Phodopus_roborovskii",
    "Cricetulus_griseus",
    "Cricetus_cricetus",
    "Dicrostonyx_torquatus",
    "Peromyscus_maniculatus_North",
    "Microtus_arvalis",
    "Microtus_pennsylvanicus",
    "Arvicola_amphibius",
    "Neotoma_floridana",
    "Mus_musculus",
    "Rattus_norvegicus",
    "Peromyscus_californicus",
    "Onychomys_torridus",
}

# All 8 blocks are on the same real reference scaffold/coordinate range.
REAL_SCAFFOLD = "CM000994.3"
REAL_START = 125888838
REAL_END = 125891358  # last block's start (125891324) + size (34)


def run(cmd_args):
    cmd = [sys.executable, "-m", "mafutils"] + cmd_args
    return subprocess.run(cmd, check=False, cwd=REPO_ROOT, capture_output=True, text=True)


def read_overall(path):
    with open(path, "r", encoding="utf-8") as fp:
        lines = fp.readlines()[1:]
    return dict(line.rstrip("\n").split("\t") for line in lines)


def test_stats_runs_cleanly_on_real_data(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    result = run(["stats", MAF_FILE, INDEX_FILE, "-o", out_prefix])
    assert result.returncode == 0, result.stderr

    overall = read_overall(out_prefix + ".overall.tsv")
    assert int(overall["total_blocks"]) == 8
    assert int(overall["total_species"]) == 15
    assert int(overall["total_parse_errors"]) == 0
    assert int(overall["blocks_with_zero_parsed_seq_lines_but_index_nonzero"]) == 0


def test_stats_species_list_matches_known_species(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    result = run(["stats", MAF_FILE, INDEX_FILE, "-o", out_prefix])
    assert result.returncode == 0, result.stderr

    with open(out_prefix + ".species.tsv", "r", encoding="utf-8") as fp:
        reader = csv.DictReader(fp, delimiter="\t")
        got_species = {row["species"] for row in reader}
    assert got_species == EXPECTED_SPECIES


def test_gc_runs_cleanly_and_plausibly_on_real_data(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    result = run(["gc", MAF_FILE, INDEX_FILE, "-o", out_prefix])
    assert result.returncode == 0, result.stderr

    with open(out_prefix + ".gc.csv", "r", encoding="utf-8") as fp:
        reader = csv.DictReader(fp)
        gc_values = {row["species"]: float(row["gc"]) for row in reader}

    assert gc_values.keys() == EXPECTED_SPECIES
    for species, gc in gc_values.items():
        # Real rodent genomic GC content clusters tightly around 0.40-0.45 in
        # this fixture; 0.30-0.55 gives real margin without being so loose it
        # would miss a genuine encoding/parsing bug (e.g. counting gaps as GC).
        assert 0.30 <= gc <= 0.55, f"{species}: implausible GC content {gc}"


def test_fetch_runs_cleanly_on_real_data(tmp_path):
    bed_path = os.path.join(str(tmp_path), "region.bed")
    with open(bed_path, "w", encoding="utf-8") as fp:
        fp.write(f"{REAL_SCAFFOLD}\t{REAL_START}\t{REAL_END}\tregion1\n")

    out_dir = str(tmp_path)
    result = run([
        "fetch", MAF_FILE, bed_path,
        "--index", INDEX_FILE,
        "-b", "id",
        "-o", out_dir,
    ])
    assert result.returncode == 0, result.stderr

    output_file = os.path.join(out_dir, "region1.maf")
    assert os.path.isfile(output_file), "Expected fetch output file was not created"
    with open(output_file, "r", encoding="utf-8") as fp:
        content = fp.read()
    assert content.count("\na") + content.startswith("a") >= 1, "Expected at least one alignment block in fetch output"
    assert "Mus_musculus" in content
