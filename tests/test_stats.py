import csv
import os
import subprocess
import sys

import pytest

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")

# Hand-verified directly against tests/example.maf's raw content (10 blocks;
# human present in all 10, chimp in 8 [blocks 1,2,3,4,5,6,8,10], gorilla in 2
# [blocks 7,9]). Every block has exactly 2 taxa, so parsimony-informative
# sites (requiring >=2 states each seen >=2 times) are structurally
# impossible; variable sites occur only in block 3 (6 columns) and block 4
# (1 column), the only two blocks where both taxa have different non-gap
# bases in the same column.
EXPECTED_OVERALL = {
    "total_blocks": 10,
    "total_species": 3,
    "total_alignment_columns": 61,
    "total_sequence_lines": 20,
    "all_gap_sites": 0,
    "variable_sites": 7,
    "invariant_sites": 54,
    "parsimony_informative_sites": 0,
    "total_parse_errors": 0,
    "blocks_with_zero_parsed_seq_lines_but_index_nonzero": 0,
}

# blocks_present / gaps_total / nongaps_total per species, hand-counted
# directly from each block's raw sequence strings.
EXPECTED_SPECIES_RAW = {
    "human": {"blocks_present": 10, "gaps_total": 3, "nongaps_total": 58},
    "chimp": {"blocks_present": 8, "gaps_total": 11, "nongaps_total": 41},
    "gorilla": {"blocks_present": 2, "gaps_total": 1, "nongaps_total": 8},
}


def run_stats(args):
    cmd = [sys.executable, "-m", "mafutils", "stats"] + args
    result = subprocess.run(cmd, check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, f"stats failed (exit {result.returncode}):\n{result.stderr}"
    return result


def read_overall(path):
    with open(path, "r", encoding="utf-8") as fp:
        lines = fp.readlines()[1:]
    return dict(line.rstrip("\n").split("\t") for line in lines)


def read_species(path):
    with open(path, "r", encoding="utf-8") as fp:
        reader = csv.DictReader(fp, delimiter="\t")
        return {row["species"]: row for row in reader}


def assert_overall_matches(overall):
    for key, expected in EXPECTED_OVERALL.items():
        assert int(overall[key]) == expected, f"{key}: expected {expected}, got {overall[key]}"


def assert_species_matches(species):
    assert species.keys() == EXPECTED_SPECIES_RAW.keys()
    for sp, expected in EXPECTED_SPECIES_RAW.items():
        row = species[sp]
        assert int(row["blocks_present"]) == expected["blocks_present"], sp
        assert int(row["gaps_total"]) == expected["gaps_total"], sp
        assert int(row["nongaps_total"]) == expected["nongaps_total"], sp

        # Derived ratios computed from the hand-verified raw counts above,
        # not hardcoded floats -- still an independent check of the formula
        # application, without risking manual 8-decimal arithmetic errors.
        blocks_present = expected["blocks_present"]
        blocks_missing = 10 - blocks_present
        total_cells = expected["gaps_total"] + expected["nongaps_total"]
        assert int(row["blocks_missing"]) == blocks_missing, sp
        assert float(row["presence_fraction"]) == pytest.approx(blocks_present / 10, abs=1e-6), sp
        assert float(row["gap_fraction_cells"]) == pytest.approx(expected["gaps_total"] / total_cells, abs=1e-6), sp
        assert float(row["avg_gaps_per_present_block"]) == pytest.approx(
            expected["gaps_total"] / blocks_present, abs=1e-6
        ), sp


def test_overall_stats_match_hand_verified_values(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    run_stats([MAF_FILE, INDEX_FILE, "-o", out_prefix])
    assert_overall_matches(read_overall(out_prefix + ".overall.tsv"))


def test_species_stats_match_hand_verified_values(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    run_stats([MAF_FILE, INDEX_FILE, "-o", out_prefix])
    assert_species_matches(read_species(out_prefix + ".species.tsv"))


def test_no_block_table_omits_block_file(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    run_stats([MAF_FILE, INDEX_FILE, "-o", out_prefix, "--no-block-table"])
    assert not os.path.exists(out_prefix + ".block.tsv")


def test_block_table_present_by_default(tmp_path):
    out_prefix = os.path.join(str(tmp_path), "out")
    run_stats([MAF_FILE, INDEX_FILE, "-o", out_prefix])
    assert os.path.isfile(out_prefix + ".block.tsv")


def test_sequential_and_parallel_paths_produce_identical_output(tmp_path):
    """
    Regression guard for the processes=1 ProcessPoolExecutor-skip fix: the
    direct in-process call (-p 1) and the real multi-worker pool path
    (-p 2, forced into multiple chunks via a small --chunk-size) must
    produce identical overall/species output, not just "doesn't crash".
    """
    seq_prefix = os.path.join(str(tmp_path), "seq")
    par_prefix = os.path.join(str(tmp_path), "par")
    run_stats([MAF_FILE, INDEX_FILE, "-o", seq_prefix, "-p", "1"])
    run_stats([MAF_FILE, INDEX_FILE, "-o", par_prefix, "-p", "2", "--chunk-size", "3"])

    assert read_overall(seq_prefix + ".overall.tsv") == read_overall(par_prefix + ".overall.tsv")
    assert read_species(seq_prefix + ".species.tsv") == read_species(par_prefix + ".species.tsv")


def test_expected_species_file_reports_missing_species(tmp_path):
    species_file = os.path.join(str(tmp_path), "species.txt")
    with open(species_file, "w", encoding="utf-8") as fp:
        fp.write("human\nchimp\ngorilla\norangutan\n")

    out_prefix = os.path.join(str(tmp_path), "out")
    run_stats([MAF_FILE, INDEX_FILE, "-o", out_prefix, "--expected-species-file", species_file])

    overall = read_overall(out_prefix + ".overall.tsv")
    assert int(overall["total_species"]) == 4  # 3 observed + 1 expected-but-absent

    species = read_species(out_prefix + ".species.tsv")
    assert "orangutan" in species
    assert int(species["orangutan"]["blocks_present"]) == 0
    assert int(species["orangutan"]["blocks_missing"]) == 10
