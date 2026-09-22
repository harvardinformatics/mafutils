#############################################################################
# Print a quick overview of a MAF file from its index.
#
# Gregg Thomas + Claude, September 2026
#############################################################################

"""
mafutils info

Prints a fast overview of a MAF: size, compression, index freshness, and the
block/scaffold counts recorded in the index header.

Deliberately cheap. Three tiers, in order of cost, and the MAF itself is never
scanned for these numbers:

  1. Index header carries the aggregate keys (built by mafutils >= 0.7):
     header read only, instant regardless of MAF size.
  2. Keys absent (older index): fall back to streaming the *index*, which is
     orders of magnitude smaller than the MAF -- ~52 MB/s, so ~9s for a 445MB
     index but ~12 min for a 35GB one. Warns up front with the index size and
     recommends rebuilding, so a long wait is never a surprise.
  3. No index: report only what needs no index and point at `mafutils index`.

Species are SAMPLED from the first N blocks (default 1000) and always labelled
as such: species names are not recorded in the index, and collecting them
during indexing would parse every non-reference `s` line, measured at 2.06x
slower indexing -- undoing most of the speedup in v0.6.0. For an exhaustive
per-species picture use `mafutils stats`, which writes <prefix>.species.tsv.
"""

import logging
import os
import sys
from typing import Annotated, Optional

import typer

from mafutils.lib import common as COMMON

DEFAULT_SAMPLE_BLOCKS = 1000
SPECIES_PREVIEW = 20

#############################################################################


def humanBytes(n):
    """
    Decimal units, so "GB" means 10^9 bytes and the printed number matches the
    `size=` field in the index header rather than being 7% smaller.
    """
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if n < 1000 or unit == "TB":
            return f"{n:.2f} {unit}" if unit != "B" else f"{n} B"
        n /= 1000


def readAggregates(header):
    """
    Pulls the aggregate counts out of a parsed index header, or None if this
    index predates them. All-or-nothing on purpose: a partial set would mean
    reporting some counts and silently omitting others.
    """
    if not header:
        return None
    if not all(key in header for key in COMMON.INDEX_AGGREGATE_KEYS):
        return None
    try:
        return {key: int(header[key]) for key in COMMON.INDEX_AGGREGATE_KEYS}
    except ValueError:
        return None


def scanIndexAggregates(index_file):
    """
    Recomputes the aggregates by streaming the block index. Only used for
    indexes built before the header recorded them.

    Streams and accumulates -- never materializes the rows. A whole-genome
    index is 35GB and ~233M rows, so holding them would repeat the OOM this
    project already fixed in gc/stats.
    """
    blocks = ref_bases = aln_cols = seq_lines = max_seqs = 0
    scaffolds = set()
    with open(index_file, "r", encoding="utf-8") as fp:
        for line in fp:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                continue
            try:
                ref_len = int(fields[2])
                cols = int(fields[3])
                nseq = int(fields[5])
            except ValueError:
                continue
            blocks += 1
            ref_bases += ref_len
            aln_cols += cols
            seq_lines += nseq
            if nseq > max_seqs:
                max_seqs = nseq
            scaffolds.add(fields[0])
    return {
        "blocks": blocks,
        "scaffolds": len(scaffolds),
        "ref_bases": ref_bases,
        "aln_cols": aln_cols,
        "seq_lines": seq_lines,
        "max_seqs": max_seqs,
    }


def sampleSpecies(maf_file, compression, max_blocks):
    """
    Collects species names from the first `max_blocks` blocks. Returns
    (sorted names, blocks actually seen).

    Sequence lines are detected with split()[0] == "s", never
    startswith("s ") -- real MAFs are tab-delimited and the space form matches
    nothing on them, a false negative this project has already been bitten by.
    """
    species = set()
    blocks = 0
    with COMMON.openMaf(maf_file, compression, "rt") as fp:
        for line in fp:
            if not line.strip():
                continue
            fields = line.split()
            if fields[0] == "a":
                blocks += 1
                if blocks > max_blocks:
                    blocks -= 1
                    break
            elif fields[0] == "s" and len(fields) > 1:
                # src is "species.scaffold"; the species is everything before
                # the first dot (mirrors fetch.speciesFromSrc, which lives in
                # fetch.py rather than common).
                species.add(fields[1].split(".", 1)[0])
    return sorted(species), blocks


#############################################################################


def run_info(args, LOG):
    maf_compression = COMMON.detectCompression(args.maf_file)
    if maf_compression not in ("none", "gz", "bgzip"):
        LOG.error(f"Unsupported MAF compression: {maf_compression}")
        sys.exit(1)

    index_file = args.index_file or COMMON.deriveBlockIndexPath(args.maf_file)
    scaffold_index_file = args.scaffold_index_file or COMMON.deriveScaffoldIndexPath(args.maf_file)
    have_index = os.path.isfile(index_file)

    header = COMMON.readIndexHeader(index_file) if have_index else None
    aggregates = readAggregates(header)

    if have_index and aggregates is None:
        index_size = os.path.getsize(index_file)
        LOG.warning(
            "Index predates the aggregate header fields; computing block stats by "
            "scanning the index (%s). Rebuild with `mafutils index` to make this instant.",
            humanBytes(index_size),
        )
        aggregates = scanIndexAggregates(index_file)

    # Scaffold names come from the scaffold index, which is tiny (a couple of
    # KB even for a whole genome) -- worth reading whenever it exists.
    scaffold_names = []
    if os.path.isfile(scaffold_index_file):
        seen = []
        with open(scaffold_index_file, "r", encoding="utf-8") as fp:
            for line in fp:
                if line.startswith("#") or not line.strip():
                    continue
                name = line.split("\t", 1)[0]
                if name not in seen:
                    seen.append(name)
        scaffold_names = seen

    species, sampled_blocks = ([], 0)
    if args.sample_blocks > 0:
        species, sampled_blocks = sampleSpecies(args.maf_file, maf_compression, args.sample_blocks)

    # ---- output -------------------------------------------------------
    size = os.path.getsize(args.maf_file)
    out = []
    out.append(f"{'MAF':<15}: {os.path.basename(args.maf_file)}")
    out.append(f"{'size':<15}: {humanBytes(size)} ({maf_compression})")

    if not have_index:
        out.append(f"{'index':<15}: not found at {index_file}")
    else:
        comparison = COMMON.compareIndexHeader(header, args.maf_file, maf_compression)
        problems = [k for k, v in comparison.items() if v is not None and not v["match"]]
        if header is None:
            state = "present, but has no mafutils header (built by an older version)"
        elif problems:
            state = "STALE -- " + ", ".join(f"{k} differs" for k in problems)
        else:
            state = "matches MAF (size + mtime)"
        out.append(f"{'index':<15}: {state}")

    out.append("")
    if scaffold_names:
        out.append(f"{'ref scaffolds':<15}: {len(scaffold_names):,}")
    elif aggregates:
        out.append(f"{'ref scaffolds':<15}: {aggregates['scaffolds']:,}")

    if aggregates:
        out.append(f"{'blocks':<15}: {aggregates['blocks']:,}")
        out.append(f"{'ref bases':<15}: {aggregates['ref_bases']:,}")
        out.append(f"{'aln columns':<15}: {aggregates['aln_cols']:,}")
        out.append(f"{'seq lines':<15}: {aggregates['seq_lines']:,}   (all species, including the reference)")
        out.append(f"{'max seqs/block':<15}: {aggregates['max_seqs']:,}")
    else:
        out.append(f"{'blocks':<15}: unavailable without an index -- run `mafutils index` first")

    if args.sample_blocks > 0:
        out.append("")
        out.append(
            f"{'species':<15}: {len(species):,}  "
            f"(sampled from the first {sampled_blocks:,} block(s) -- possibly incomplete)"
        )
        shown = species[:SPECIES_PREVIEW]
        line = "  "
        for name in shown:
            if len(line) + len(name) + 2 > 78:
                out.append(line.rstrip())
                line = "  "
            line += name + ", "
        out.append(line.rstrip().rstrip(","))
        if len(species) > SPECIES_PREVIEW:
            out.append(f"  ... and {len(species) - SPECIES_PREVIEW:,} more (--all-species to list them)")
        out.append("")
        out.append("For exact per-species presence across the whole file, run `mafutils stats`.")

    print("\n".join(out))


#############################################################################


def info_command(
    maf_file: Annotated[str, typer.Argument(help="Input MAF file (.maf, .maf.gz, or bgzip-compressed .maf)")],
    index_file: Annotated[Optional[str], typer.Option("--index", "-i", help="Block index (default: <MAF_FILE>.block.idx)")] = None,
    scaffold_index_file: Annotated[Optional[str], typer.Option("--scaffold-index", help="Scaffold index (default: <MAF_FILE>.scaffold.idx)")] = None,
    sample_blocks: Annotated[int, typer.Option("--sample-blocks", help="Blocks to sample for species names; 0 disables the sample and reads no MAF data at all.")] = DEFAULT_SAMPLE_BLOCKS,
    all_species: Annotated[bool, typer.Option("--all-species", help="List every sampled species name instead of the first 20.")] = False,
) -> None:
    if sample_blocks < 0:
        print("--sample-blocks must be >= 0", file=sys.stderr)
        raise typer.Exit(1)

    global SPECIES_PREVIEW
    if all_species:
        SPECIES_PREVIEW = sys.maxsize

    logging.basicConfig(format="[ %(levelname)s ] %(message)s", level=logging.INFO)
    LOG = logging.getLogger("maf_info_logger")

    from types import SimpleNamespace
    args = SimpleNamespace(
        maf_file=maf_file,
        index_file=index_file,
        scaffold_index_file=scaffold_index_file,
        sample_blocks=sample_blocks,
    )
    run_info(args, LOG)
