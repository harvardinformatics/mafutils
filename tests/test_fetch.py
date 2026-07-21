import os
import subprocess
import sys
import pytest

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))
BED_FILE = os.path.join(TEST_DIR, "example.bed")
#BED_FILE = os.path.join(TEST_DIR, "crossblocks-missing-species.bed")
MAF_FILE = os.path.join(TEST_DIR, "example.maf")
INDEX_FILE = os.path.join(TEST_DIR, "example.maf.block.idx")
EXPECTED_MAF_DIR = os.path.join(TEST_DIR, "expected-maf")
EXPECTED_FASTA_DIR = os.path.join(TEST_DIR, "expected-fasta")

# Regions with no overlapping alignment block are a normal, expected outcome
# (fetch typically processes many regions per run) -- fetch warns and moves
# on rather than failing. This is distinct from 09-missingchrom, where the
# BED references a scaffold that isn't in the index at all, which is a
# harder error and does exit non-zero (see the empty-expected-file branch
# below).
WARN_ONLY_NO_OVERLAP_REGIONS = {"06-noalign"}

def get_regions():
    """Yields (chrom, start, end, id) for each region in the BED file."""
    with open(BED_FILE) as bed:
        for line in bed:
            if line.strip() and not line.startswith("#"):
                fields = line.strip().split()
                yield fields[0], int(fields[1]), int(fields[2]), fields[3]

@pytest.mark.parametrize("chrom,start,end,region_id", list(get_regions()))
def test_maf_fetch_region(chrom, start, end, region_id, tmp_path):
    """Run maf_fetch for just one region and compare output to expected."""
    out_dir = str(tmp_path)
    # Write a one-line BED for this region
    per_test_bed = os.path.join(out_dir, f"{region_id}.bed")
    os.makedirs(out_dir, exist_ok=True)
    with open(per_test_bed, "w") as out:
        out.write(f"{chrom}\t{start}\t{end}\t{region_id}\n")

    # Run maf_fetch
    output_file = os.path.join(out_dir, f"{region_id}.maf")
    cmd = [
        sys.executable, "-m", "mafutils", "fetch",
        MAF_FILE,
        per_test_bed,
        "--index", INDEX_FILE,
        "-b", "id",
        "-o", out_dir
    ]
    expected_file = os.path.join(EXPECTED_MAF_DIR, f"{region_id}.maf")
    result = subprocess.run(cmd, check=False, cwd=REPO_ROOT, capture_output=True, text=True)

    if region_id in WARN_ONLY_NO_OVERLAP_REGIONS:
        assert result.returncode == 0, f"Expected exit 0 (warn, not fail) for region {region_id}: {result.stderr}"
        assert not os.path.exists(output_file), f"Did not expect output file for region {region_id}"
        assert "no-overlap" in result.stderr, f"Expected a no-overlap warning for region {region_id}"
        return

    if os.path.getsize(expected_file) == 0:
        assert result.returncode != 0, f"Expected non-zero exit for region {region_id}"
        assert not os.path.exists(output_file), f"Did not expect output file for region {region_id}"
        return

    assert result.returncode == 0, f"Unexpected non-zero exit for region {region_id}: {result.returncode}"

    got = open(output_file).read()
    expected = open(expected_file).read()
    
    # Strip trailing blank lines from both outputs before comparison
    got = got.rstrip('\n')
    expected = expected.rstrip('\n')

    assert got == expected, f"Output for region {region_id} does not match expected"

@pytest.mark.parametrize("chrom,start,end,region_id", list(get_regions()))
def test_maf_fetch_fasta_region(chrom, start, end, region_id, tmp_path):
    """Run maf_fetch for just one region (FASTA output) and compare to expected."""
    out_dir = str(tmp_path)
    per_test_bed = os.path.join(out_dir, f"{region_id}.bed")
    os.makedirs(out_dir, exist_ok=True)
    with open(per_test_bed, "w") as out:
        out.write(f"{chrom}\t{start}\t{end}\t{region_id}\n")

    # Run maf_fetch with --fasta/-f
    output_file = os.path.join(out_dir, f"{region_id}.fa")
    cmd = [
        sys.executable, "-m", "mafutils", "fetch",
        MAF_FILE,
        per_test_bed,
        "--index", INDEX_FILE,
        "-b", "id",
        "-o", out_dir,
        "-f",
        "-fh", "species-coords-id"
    ]
    expected_file = os.path.join(EXPECTED_FASTA_DIR, f"{region_id}.fa")
    result = subprocess.run(cmd, check=False, cwd=REPO_ROOT, capture_output=True, text=True)

    if region_id in WARN_ONLY_NO_OVERLAP_REGIONS:
        assert result.returncode == 0, f"Expected exit 0 (warn, not fail) for region {region_id}: {result.stderr}"
        assert not os.path.exists(output_file), f"Did not expect FASTA output file for region {region_id}"
        assert "no-overlap" in result.stderr, f"Expected a no-overlap warning for region {region_id}"
        return

    if os.path.getsize(expected_file) == 0:
        assert result.returncode != 0, f"Expected non-zero FASTA exit for region {region_id}"
        assert not os.path.exists(output_file), f"Did not expect FASTA output file for region {region_id}"
        return

    assert result.returncode == 0, f"Unexpected non-zero FASTA exit for region {region_id}: {result.returncode}"

    got = open(output_file).read()
    expected = open(expected_file).read()

    # Strip trailing blank lines from both outputs before comparison
    got = got.rstrip('\n')
    expected = expected.rstrip('\n')

    assert got == expected, f"FASTA output for region {region_id} does not match expected"


def test_scaffold_subdirs_groups_output_by_scaffold(tmp_path):
    """
    --scaffold-subdirs groups per-region output files into
    <outdir>/<scaffold>/<basename>.maf instead of a flat <outdir>/. Runs
    with -p 2 (real multi-worker ProcessPoolExecutor dispatch, matching how
    block-mode writes actually happen) across regions spanning multiple
    scaffolds, including two scaffolds (chr1, chr4) with multiple regions
    each, to confirm they land in the same subdirectory without any
    cross-worker race (directories are pre-created upfront, see fetch.py).
    """
    out_dir = str(tmp_path)
    cmd = [
        sys.executable, "-m", "mafutils", "fetch",
        MAF_FILE, BED_FILE,
        "--index", INDEX_FILE,
        "-b", "id",
        "-o", out_dir,
        "--scaffold-subdirs",
        "-p", "2",
    ]
    result = subprocess.run(cmd, check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode == 0, f"Unexpected non-zero exit: {result.stderr}"

    # chr1 has 3 regions with real output (06-noalign has no overlap, so no
    # file); chr4 has 2; chr2/chr3/chrX have 1 each. chrZ (09-missingchrom)
    # never gets that far (unknown scaffold errors before any writes).
    assert sorted(os.listdir(os.path.join(out_dir, "chr1"))) == [
        "01-singleblock-gap.maf",
        "02-truncatedblock.maf",
        "03-crossblocks-missing.maf",
    ]
    assert os.listdir(os.path.join(out_dir, "chr2")) == ["04-singleblock.maf"]
    assert os.listdir(os.path.join(out_dir, "chr3")) == ["05-crossblocks.maf"]
    assert os.listdir(os.path.join(out_dir, "chrX")) == ["07-negstrand.maf"]
    assert sorted(os.listdir(os.path.join(out_dir, "chr4"))) == [
        "08-span-multiple-gaps.maf",
        "10-crossblocks-missing-species.maf",
    ]

    # Content should be identical to the flat-layout output for the same region.
    with open(os.path.join(out_dir, "chr2", "04-singleblock.maf")) as fp:
        got = fp.read().rstrip("\n")
    with open(os.path.join(EXPECTED_MAF_DIR, "04-singleblock.maf")) as fp:
        expected = fp.read().rstrip("\n")
    assert got == expected


def test_scaffold_subdirs_rejected_with_single_output(tmp_path):
    cmd = [
        sys.executable, "-m", "mafutils", "fetch",
        MAF_FILE, BED_FILE,
        "--index", INDEX_FILE,
        "--scaffold-subdirs",
        "--single-output",
        "-o", os.path.join(str(tmp_path), "out.maf"),
    ]
    result = subprocess.run(cmd, check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert result.returncode != 0
    assert "--scaffold-subdirs cannot be used with --single-output" in result.stderr
