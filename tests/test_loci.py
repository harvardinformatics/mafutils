import csv
import os
import subprocess
import sys
from collections import defaultdict

import pytest

from mafutils.fetch import (ELEMENT_LOCI_HEADERS, chunkGaps, classifyChunks,
                            fastaScaffoldField, summarizeElementLoci)

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))

MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")
BED_FILE = os.path.join(TEST_DIR, "example.bed")

# 6 real blocks (15 species) covering two real regions chosen because they
# exhibit the rare classes that tests/example.maf cannot produce:
# multi_strand and multi_scaffold. Extracted read-only from production data.
CLASSES_MAF = os.path.join(TEST_DIR, "loci-classes.maf")
CLASSES_INDEX = CLASSES_MAF + ".block.idx"
CLASSES_REGIONS = [
    ("CM000994.3", 3541891, 3541950, "multistrand"),
    ("CM000994.3", 3602239, 3602958, "multiscaffold"),
]


def run_fetch(args):
    return subprocess.run(
        [sys.executable, "-m", "mafutils", "fetch"] + args,
        check=False, cwd=REPO_ROOT, capture_output=True, text=True,
    )


def write_bed(path, rows):
    with open(path, "w") as fp:
        for r in rows:
            fp.write("\t".join(str(x) for x in r) + "\n")


def read_loci(path):
    with open(path) as fp:
        return list(csv.DictReader(fp, delimiter="\t"))


def headers(fa_path):
    with open(fa_path) as fp:
        return [ln.strip() for ln in fp if ln.startswith(">")]


# ---------------------------------------------------------------------------
# classifyChunks / fastaScaffoldField (unit level)
# ---------------------------------------------------------------------------


def chunk(scaffold="s1", start=0, size=10, strand="+", src_size=1000, bi=1):
    """bi = block_index, the block ordinal this chunk came from."""
    return {"scaffold": scaffold, "start": start, "size": size,
            "end": start + size, "strand": strand, "src_size": src_size,
            "block_index": bi}


def test_single_chunk_is_single():
    assert classifyChunks([chunk()]) == ("single", 0)


def test_adjacent_chunks_are_contiguous():
    chunks = [chunk(start=0, size=10), chunk(start=10, size=5)]
    assert classifyChunks(chunks) == ("contiguous", 0)


def test_gap_makes_it_split_and_reports_the_gap():
    chunks = [chunk(start=0, size=10), chunk(start=60, size=5)]
    assert classifyChunks(chunks) == ("split", 50)


def test_largest_magnitude_gap_wins_and_sign_is_kept():
    """
    Overlapping chunks (negative gap) mean something different from a plain
    gap, so the sign must survive rather than being abs()'d away.
    """
    chunks = [chunk(start=0, size=10), chunk(start=5, size=10), chunk(start=100, size=5)]
    cls, gap = classifyChunks(chunks)
    assert cls == "split"
    assert gap == 85                      # 100 - (5+10), the largest magnitude

    overlapping = [chunk(start=0, size=50), chunk(start=10, size=5)]
    cls, gap = classifyChunks(overlapping)
    assert cls == "split"
    assert gap == -40                     # signed, not 40


def test_differing_scaffolds_beat_other_signals():
    chunks = [chunk(scaffold="s1", start=0, size=10), chunk(scaffold="s2", start=10, size=5)]
    assert classifyChunks(chunks) == ("multi_scaffold", None)


def test_differing_strands():
    chunks = [chunk(start=0, size=10, strand="+"), chunk(start=10, size=5, strand="-")]
    assert classifyChunks(chunks) == ("multi_strand", None)


def test_scaffold_field_dedupes_preserving_block_order():
    chunks = [chunk(scaffold="b"), chunk(scaffold="a"), chunk(scaffold="b")]
    assert fastaScaffoldField(chunks) == "b,a"
    assert fastaScaffoldField([chunk(scaffold="")]) == ""


# ---------------------------------------------------------------------------
# The scaffold must be in the header regardless of dedupe / expected-species
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("extra", [
    [],
    ["--fasta-dedupe", "most-seq"],
    ["--expected-species", "human,chimp,gorilla"],
    ["--fasta-dedupe", "most-seq", "--expected-species", "human,chimp,gorilla"],
])
def test_scaffold_present_in_header_for_every_dedupe_mode(extra, tmp_path):
    """
    The scaffold used to live only in the FASTA dict key, so keying by species
    (which --fasta-dedupe and --expected-species both force) discarded it and
    the header could not name a locus.
    """
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [("chr2", 5, 10, "r1")])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id"] + extra)
    assert result.returncode == 0, result.stderr

    hdrs = headers(os.path.join(out, "r1.fa"))
    present = [h for h in hdrs if h.startswith((">human", ">chimp"))]
    assert present, hdrs
    for h in present:
        assert h.startswith((">human.chr2", ">chimp.chr2")), h


def test_species_only_stays_bare(tmp_path):
    """species-only exists so headers match tree tip labels; it never gains a scaffold."""
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [("chr2", 5, 10, "r1")])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-b", "id", "-o", out,
                        "-f", "-fh", "species-only", "--fasta-dedupe", "most-seq"])
    assert result.returncode == 0, result.stderr
    assert sorted(headers(os.path.join(out, "r1.fa"))) == [">chimp", ">human"]


def test_species_coords_id_survives_a_bed_with_no_name_column(tmp_path):
    """Regression: id_str was assigned conditionally but used unconditionally."""
    bed = os.path.join(str(tmp_path), "noid.bed")
    write_bed(bed, [("chr2", 5, 10)])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-o", out,
                        "-f", "-fh", "species-coords-id"])
    assert result.returncode == 0, result.stderr
    # default --basename is 'coords', so the file is named scaffold-start-end
    for h in headers(os.path.join(out, "chr2-5-10.fa")):
        assert not h.endswith(" "), f"trailing space with no region id: {h!r}"
        assert "id:" not in h


# ---------------------------------------------------------------------------
# The loci table
# ---------------------------------------------------------------------------


def test_loci_table_requires_fasta(tmp_path):
    result = run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id",
                        "-o", os.path.join(str(tmp_path), "out"), "--loci-table"])
    assert result.returncode != 0
    assert "--loci-table requires --fasta" in (result.stdout + result.stderr)


def test_loci_table_grain_and_class_consistency(tmp_path):
    out = os.path.join(str(tmp_path), "out")
    result = run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id", "--loci-table"])
    assert result.returncode == 0, result.stderr

    rows = read_loci(os.path.join(out, "maf_fetch_loci.tsv"))
    assert rows

    groups = defaultdict(list)
    for r in rows:
        groups[(r["region_scaffold"], r["region_start"], r["region_end"], r["species"])].append(r)

    for key, g in groups.items():
        # one row per contributing chunk, indexed 1..n, class constant per group
        assert len({r["class"] for r in g}) == 1, key
        assert len({r["max_gap"] for r in g}) == 1, key
        assert int(g[0]["n_chunks"]) == len(g), key
        assert sorted(int(r["chunk_index"]) for r in g) == list(range(1, len(g) + 1)), key
        if len(g) == 1:
            assert g[0]["class"] == "single", key


def test_loci_table_matches_the_examples_known_shape(tmp_path):
    """
    03-crossblocks-missing is a hand-checkable case: human's two chunks are 6bp
    apart (split) while chimp's are adjacent (contiguous) in the same region.
    """
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [("chr1", 24, 33, "cross")])
    out = os.path.join(str(tmp_path), "out")
    result = run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id", "--loci-table"])
    assert result.returncode == 0, result.stderr

    rows = read_loci(os.path.join(out, "maf_fetch_loci.tsv"))
    by_sp = defaultdict(list)
    for r in rows:
        by_sp[r["species"]].append(r)

    assert by_sp["human"][0]["class"] == "split"
    assert by_sp["human"][0]["max_gap"] == "6"
    assert by_sp["chimp"][0]["class"] == "contiguous"
    assert by_sp["chimp"][0]["max_gap"] == "0"
    # src_size is emitted so forward coordinates are one subtraction away
    assert all(int(r["src_size"]) > 0 for r in rows)


@pytest.mark.parametrize("processes", [1, 4])
def test_loci_table_is_deterministic_across_process_counts(processes, tmp_path):
    base = os.path.join(str(tmp_path), "p1")
    assert run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", base,
                      "-f", "-fh", "species-coords-id", "--loci-table", "-p", "1"]).returncode == 0

    out = os.path.join(str(tmp_path), f"p{processes}")
    result = run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id", "--loci-table", "-p", str(processes)])
    assert result.returncode == 0, result.stderr
    with open(os.path.join(out, "maf_fetch_loci.tsv")) as a, open(os.path.join(base, "maf_fetch_loci.tsv")) as b:
        assert a.read() == b.read()


def test_loci_table_does_not_change_the_fasta(tmp_path):
    without = os.path.join(str(tmp_path), "without")
    with_ = os.path.join(str(tmp_path), "with")
    common = [MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-f", "-fh", "species-coords-id"]
    assert run_fetch(common + ["-o", without]).returncode == 0
    assert run_fetch(common + ["-o", with_, "--loci-table"]).returncode == 0

    for name in os.listdir(without):
        if name.endswith(".fa"):
            with open(os.path.join(without, name)) as a, open(os.path.join(with_, name)) as b:
                assert a.read() == b.read(), name


# ---------------------------------------------------------------------------
# The rare classes, on real data
# ---------------------------------------------------------------------------


def test_real_fixture_produces_the_rare_classes(tmp_path):
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, CLASSES_REGIONS)
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([CLASSES_MAF, bed, "--index", CLASSES_INDEX, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id", "--fasta-dedupe", "most-seq",
                        "--loci-table"])
    assert result.returncode == 0, result.stderr

    rows = read_loci(os.path.join(out, "maf_fetch_loci.tsv"))
    seen = {r["class"] for r in rows}
    assert {"multi_scaffold", "multi_strand", "split", "contiguous"} <= seen, seen

    # max_gap is undefined across different source sequences
    for r in rows:
        if r["class"] in ("multi_scaffold", "multi_strand"):
            assert r["max_gap"] == ".", r


def test_multi_scaffold_header_comma_joins_and_warns(tmp_path):
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [CLASSES_REGIONS[1]])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([CLASSES_MAF, bed, "--index", CLASSES_INDEX, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id", "--fasta-dedupe", "most-seq",
                        "--loci-table", "--verbose"])
    assert result.returncode == 0, result.stderr

    rows = read_loci(os.path.join(out, "maf_fetch_loci.tsv"))
    multi = {r["species"] for r in rows if r["class"] == "multi_scaffold"}
    assert multi, "fixture no longer exercises multi_scaffold"

    hdrs = headers(os.path.join(out, "multiscaffold.fa"))
    for sp in multi:
        h = next(x for x in hdrs if x.startswith(f">{sp}."))
        scaffolds = {r["src_scaffold"] for r in rows if r["species"] == sp}
        assert "," in h, h
        for scaf in scaffolds:
            assert scaf in h, (scaf, h)

    assert "multi-scaffold" in (result.stdout + result.stderr)


# ---------------------------------------------------------------------------
# Zero-size chunks and the no_bases class
# ---------------------------------------------------------------------------


def test_all_zero_size_chunks_are_no_bases_not_single():
    """
    A species can have an s-line in a block whose aligned sequence falls
    entirely outside the region once trimmed, leaving an all-gap row and a
    zero-width header span. That used to classify as a clean `single`
    (1,648 such rows on real data).
    """
    assert classifyChunks([chunk(size=0)]) == ("no_bases", 0)
    assert classifyChunks([chunk(start=0, size=0), chunk(start=50, size=0)]) == ("no_bases", 0)


def test_zero_size_chunks_do_not_create_phantom_gaps():
    """
    Real case (Dicrostonyx_torquatus, CM000994.3:10038193): a zero-size chunk
    in one block plus 27 real bases in a later block. Measuring a gap from a
    chunk containing no bases invented a 718bp `split`; the species really has
    one aligned piece.
    """
    chunks = [chunk(start=9949990, size=0), chunk(start=9950708, size=27)]
    assert classifyChunks(chunks) == ("single", 0)
    assert chunkGaps(chunks) == []


def test_zero_size_chunk_is_still_emitted_as_a_row(tmp_path):
    """Excluded from the arithmetic, but never hidden from the table."""
    from mafutils.fetch import writeLociRows
    import io
    buf = io.StringIO()
    chunks = [chunk(start=100, size=0, bi=1), chunk(start=500, size=27, bi=2)]
    writeLociRows(buf, {"scaffold": "chr1", "start": 0, "end": 1000, "id": "r1"},
                  "sp", chunks, "single", 0)
    lines = [l for l in buf.getvalue().splitlines() if l]
    assert len(lines) == 2, lines
    assert "\t0\t" in lines[0]          # the zero-size chunk survives as a row


# ---------------------------------------------------------------------------
# Per-element summary columns
# ---------------------------------------------------------------------------


def test_element_summary_counts_and_concentration():
    per_species = {
        # two species split at the SAME boundary with near-identical gaps
        "A": [chunk(start=0, size=10, bi=1), chunk(start=43, size=5, bi=2)],
        "B": [chunk(start=0, size=10, bi=1), chunk(start=44, size=5, bi=2)],
        "C": [chunk(start=0, size=10, bi=1), chunk(start=10, size=5, bi=2)],
        "D": [chunk(size=0, bi=1)],
        "E": [chunk(start=0, size=10, bi=1)],
    }
    vals = dict(zip(ELEMENT_LOCI_HEADERS, summarizeElementLoci(per_species)))
    assert vals["n.species"] == 5
    assert vals["n.split"] == 2
    assert vals["n.contiguous"] == 1
    assert vals["n.no.bases"] == 1
    assert vals["n.single"] == 1
    assert vals["split.bases.max"] == 34            # B's total
    assert vals["split.gap.median"] == 33           # median of [33, 34]
    assert vals["split.boundaries"] == 1
    assert vals["split.max.boundary.n.species"] == 2
    assert vals["split.max.boundary.gap.spread"] == 1


def test_element_summary_separates_shared_from_diffuse():
    """
    max.boundary.n.species vs n.split is the shared-vs-independent signal:
    equal means one event, much lower means scattered.
    """
    shared = {f"sp{i}": [chunk(start=0, size=10, bi=1), chunk(start=43, size=5, bi=2)] for i in range(4)}
    v = dict(zip(ELEMENT_LOCI_HEADERS, summarizeElementLoci(shared)))
    assert v["n.split"] == 4 and v["split.max.boundary.n.species"] == 4 and v["split.boundaries"] == 1

    diffuse = {
        "a": [chunk(start=0, size=10, bi=1), chunk(start=43, size=5, bi=2)],
        "b": [chunk(start=0, size=10, bi=2), chunk(start=43, size=5, bi=3)],
        "c": [chunk(start=0, size=10, bi=3), chunk(start=43, size=5, bi=4)],
    }
    v = dict(zip(ELEMENT_LOCI_HEADERS, summarizeElementLoci(diffuse)))
    assert v["n.split"] == 3 and v["split.max.boundary.n.species"] == 1 and v["split.boundaries"] == 3


def test_summary_has_the_element_columns_and_matches_the_loci_table(tmp_path):
    out = os.path.join(str(tmp_path), "out")
    result = run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id", "--loci-table"])
    assert result.returncode == 0, result.stderr

    summary = list(csv.DictReader(open(os.path.join(out, "maf_fetch_summary.tsv")), delimiter="\t"))
    for col in ELEMENT_LOCI_HEADERS:
        assert col in summary[0], col

    # the per-element counts must equal the loci table's own grouping
    rows = read_loci(os.path.join(out, "maf_fetch_loci.tsv"))
    seen, counts = set(), defaultdict(lambda: defaultdict(int))
    for r in rows:
        key = (r["region_scaffold"], r["region_start"], r["region_end"])
        if (key, r["species"]) in seen:
            continue
        seen.add((key, r["species"]))
        counts[key][r["class"]] += 1

    checked = 0
    for r in summary:
        key = (r["scaffold"], r["start"], r["end"])
        if key not in counts:
            continue
        checked += 1
        for cls, col in (("single", "n.single"), ("contiguous", "n.contiguous"),
                         ("split", "n.split"), ("multi_scaffold", "n.multi.scaffold"),
                         ("multi_strand", "n.multi.strand"), ("no_bases", "n.no.bases")):
            assert counts[key][cls] == int(r[col]), (key, col)
    assert checked > 0


def test_element_columns_identical_in_maf_and_fasta_mode(tmp_path):
    """
    The summary columns come from trimMafBlock's chunks on the MAF path and
    from mafBlockToFasta's on the FASTA path; for a MAF with no duplicate
    copies of a species per block they must agree exactly.
    """
    maf_out = os.path.join(str(tmp_path), "maf")
    fa_out = os.path.join(str(tmp_path), "fa")
    common = [MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id"]
    assert run_fetch(common + ["-o", maf_out]).returncode == 0
    assert run_fetch(common + ["-o", fa_out, "-f", "-fh", "species-coords-id"]).returncode == 0

    def element_cols(path):
        out = {}
        for r in csv.DictReader(open(os.path.join(path, "maf_fetch_summary.tsv")), delimiter="\t"):
            out[r["basename"]] = {c: r[c] for c in ELEMENT_LOCI_HEADERS}
        return out

    assert element_cols(maf_out) == element_cols(fa_out)


@pytest.mark.parametrize("processes", [1, 4])
def test_element_columns_deterministic_across_process_counts(processes, tmp_path):
    base = os.path.join(str(tmp_path), "p1")
    assert run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", base, "-p", "1"]).returncode == 0
    out = os.path.join(str(tmp_path), f"p{processes}")
    assert run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", out, "-p", str(processes)]).returncode == 0

    def rows(path):
        return sorted(
            tuple(r[c] for c in ["basename"] + ELEMENT_LOCI_HEADERS)
            for r in csv.DictReader(open(os.path.join(path, "maf_fetch_summary.tsv")), delimiter="\t")
        )
    assert rows(out) == rows(base)


# ---------------------------------------------------------------------------
# block_index / gap_to_next
# ---------------------------------------------------------------------------


def test_loci_table_block_index_and_gap_to_next(tmp_path):
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [("chr1", 24, 33, "cross")])
    out = os.path.join(str(tmp_path), "out")
    assert run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-b", "id", "-o", out,
                      "-f", "-fh", "species-coords-id", "--loci-table"]).returncode == 0

    rows = read_loci(os.path.join(out, "maf_fetch_loci.tsv"))
    by_sp = defaultdict(list)
    for r in rows:
        by_sp[r["species"]].append(r)

    # human is split by 6bp between the region's two blocks; chimp is contiguous
    human = sorted(by_sp["human"], key=lambda r: int(r["chunk_start"]))
    assert [r["block_index"] for r in human] == ["1", "2"]
    assert human[0]["gap_to_next"] == "6"      # the gap sits at boundary 1->2
    assert human[1]["gap_to_next"] == "."      # last chunk has no next

    chimp = sorted(by_sp["chimp"], key=lambda r: int(r["chunk_start"]))
    assert chimp[0]["gap_to_next"] == "0"


def test_real_fixture_concentration_columns(tmp_path):
    """
    tests/loci-classes.maf, region 3602239-3602958: 13 species split across 3
    boundaries with 13 sharing the busiest one -- a mostly-shared interruption.
    """
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [CLASSES_REGIONS[1]])
    out = os.path.join(str(tmp_path), "out")
    assert run_fetch([CLASSES_MAF, bed, "--index", CLASSES_INDEX, "-b", "id", "-o", out,
                      "-f", "-fh", "species-coords-id", "--fasta-dedupe", "most-seq"]).returncode == 0

    r = next(csv.DictReader(open(os.path.join(out, "maf_fetch_summary.tsv")), delimiter="\t"))
    assert int(r["n.split"]) == 13
    assert int(r["split.boundaries"]) == 3
    assert int(r["split.max.boundary.n.species"]) == 13
    assert int(r["split.max.boundary.gap.spread"]) == 1654
    assert int(r["n.multi.scaffold"]) == 1


def test_no_bases_end_to_end_on_real_data(tmp_path):
    """
    End-to-end cover for `no_bases`, the class most likely to silently corrupt
    a downstream BED: filtering `single|contiguous` on the old behaviour would
    have admitted these as clean, with zero-width intervals.

    A 5bp window at the very start of the fixture's first block leaves four
    species with an s-line whose aligned bases all fall outside the trimmed
    columns -- a real all-gap row, not a constructed one.
    """
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [("CM000994.3", 3541841, 3541846, "nobases")])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([CLASSES_MAF, bed, "--index", CLASSES_INDEX, "-b", "id", "-o", out,
                        "-f", "-fh", "species-coords-id", "--fasta-dedupe", "most-seq",
                        "--loci-table"])
    assert result.returncode == 0, result.stderr

    rows = read_loci(os.path.join(out, "maf_fetch_loci.tsv"))
    no_bases = {r["species"] for r in rows if r["class"] == "no_bases"}
    assert len(no_bases) == 4, sorted(no_bases)

    # every no_bases row really is a zero-size chunk, and none is called clean
    for r in rows:
        if r["species"] in no_bases:
            assert int(r["chunk_size"]) == 0, r
            assert r["class"] not in ("single", "contiguous"), r

    # the summary count agrees with the loci table
    summary = next(csv.DictReader(open(os.path.join(out, "maf_fetch_summary.tsv")), delimiter="\t"))
    assert int(summary["n.no.bases"]) == 4
    assert int(summary["n.single"]) == 10

    # and the FASTA rows are genuinely all-gap with a zero-width header span
    seqs = {}
    name = None
    for line in open(os.path.join(out, "nobases.fa")):
        if line.startswith(">"):
            name = line.strip()
            seqs[name] = ""
        else:
            seqs[name] += line.strip()
    for header, seq in seqs.items():
        species = header[1:].split(".")[0]
        if species in no_bases:
            assert set(seq) == {"-"}, header
            start, end = header.split(":")[1].split("(")[0].split("-")
            assert start == end, header


# ---------------------------------------------------------------------------
# Reference-side summary columns
# ---------------------------------------------------------------------------


def test_summary_schema_marks_reference_columns(tmp_path):
    out = os.path.join(str(tmp_path), "out")
    assert run_fetch([MAF_FILE, BED_FILE, "--index", INDEX_FILE, "-b", "id", "-o", out]).returncode == 0

    header = open(os.path.join(out, "maf_fetch_summary.tsv")).readline().rstrip("\n").split("\t")
    assert "ref.n.overlapping.blocks" in header
    assert "ref.block.bases" in header
    # removed: it was structurally all-zeros and invited being read as a
    # per-species contiguity check, which it never was
    assert "interblock.distances" not in header
    assert "block.lengths" not in header


def test_ref_coverage_gap_warns_and_shortens_ref_block_bases(tmp_path):
    """
    03-crossblocks-missing spans two blocks 6bp apart in the reference, so the
    region holds reference positions no block covers. That used to surface only
    as an `interblock.distances` value of [6]; it is now a warning, and
    ref.block.bases falls short of the region width by exactly that much.
    """
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [("chr1", 24, 33, "cross")])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-b", "id", "-o", out, "--verbose"])
    assert result.returncode == 0, result.stderr
    combined = result.stdout + result.stderr
    assert "ref-coverage-gap" in combined
    assert "do not tile the reference" in combined

    r = next(csv.DictReader(open(os.path.join(out, "maf_fetch_summary.tsv")), delimiter="\t"))
    width = int(r["end"]) - int(r["start"])
    assert int(r["ref.block.bases"]) == width - 6


def test_no_ref_coverage_gap_warning_when_blocks_tile(tmp_path):
    """
    A region inside a single block cannot have a reference coverage gap, so the
    warning must not fire -- it would be noise on every normal run.
    """
    bed = os.path.join(str(tmp_path), "r.bed")
    write_bed(bed, [("chr2", 5, 10, "single")])
    out = os.path.join(str(tmp_path), "out")

    result = run_fetch([MAF_FILE, bed, "--index", INDEX_FILE, "-b", "id", "-o", out, "--verbose"])
    assert result.returncode == 0, result.stderr
    assert "ref-coverage-gap" not in (result.stdout + result.stderr)

    r = next(csv.DictReader(open(os.path.join(out, "maf_fetch_summary.tsv")), delimiter="\t"))
    assert int(r["ref.block.bases"]) == int(r["end"]) - int(r["start"])
