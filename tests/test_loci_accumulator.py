"""
ElementLociAccumulator must compute exactly what the v0.8.1 per-element summary
computed from every chunk held in memory. The oracle below is that code,
copied verbatim from mafutils v0.8.1 (renamed oracle_*), so this test pins the
new streaming version to the old results value for value.
"""
import random
import statistics
from collections import defaultdict

import pytest

from mafutils.fetch import ELEMENT_LOCI_HEADERS, ElementLociAccumulator, summarizeElementLoci


# ---------------------------------------------------------------------------
# Oracle: verbatim from mafutils v0.8.1 fetch.py
# ---------------------------------------------------------------------------

def oracle_classifyChunks(chunks):
    """
    Classifies what a species' contribution to one region actually is, from its
    per-block chunks. Returns (class, max_gap).

    A MAF's blocks are contiguous in the REFERENCE only. A non-reference
    species' chunks may be separated, inverted, or on different source
    sequences, and stitching concatenates them without padding -- so a stitched
    row is not necessarily one locus in that genome. Measured on real data
    (2,000 regions, 15 species): 77.8% single, 10.7% contiguous, 10.6% split,
    0.05% multi-strand, 0.8% multi-scaffold.

    max_gap is the largest-magnitude gap between consecutive chunks (sorted by
    start), reported SIGNED -- negative means chunks overlap in the source,
    which indicates something different again (paralogy, tandem repeat). It is
    0 for single/contiguous and None for multi_scaffold, where a gap between
    positions on different sequences is meaningless.

    A split gap is NOT an indel the alignment absorbed. The aligner threads
    query insertions inline as reference-row gap columns, but only up to ~37bp
    on real data (156,928 such runs observed, max 37, none >=40); past that it
    breaks the block instead. So a gap between chunks is query sequence the
    aligner *declined to align* -- real bases absent from the emitted FASTA
    while the header's bounding-box span still covers them. Consequence,
    verified with zero exceptions across 27,599 records: single/contiguous have
    header span == emitted bases, and split/multi_* never do.

    No gap threshold is applied on purpose: "how large a gap disqualifies a
    locus" is a judgment for the analysis (real median split 45bp, max 335Mb),
    so the number is reported and the caller decides.

    Chunks contributing zero aligned bases are excluded from the gap/class
    arithmetic: a species can have an s-line in a block whose aligned sequence
    falls entirely outside the requested region once trimmed, leaving an
    all-gap row. If NO chunk contributes bases the class is "no_bases" rather
    than "single" -- 5.46% of real rows are zero-size, and 1,648 of them used
    to report as a clean "single" with a zero-width header span.
    """
    if not chunks:
        return "none", None

    contributing = [c for c in chunks if c['size'] > 0]
    if not contributing:
        return "no_bases", 0
    if len(contributing) == 1:
        return "single", 0
    if len({c['scaffold'] for c in contributing}) > 1:
        return "multi_scaffold", None
    if len({c['strand'] for c in contributing}) > 1:
        return "multi_strand", None

    ordered = sorted(contributing, key=lambda c: c['start'])
    gaps = [
        ordered[i + 1]['start'] - (ordered[i]['start'] + ordered[i]['size'])
        for i in range(len(ordered) - 1)
    ]
    if all(g == 0 for g in gaps):
        return "contiguous", 0
    return "split", max(gaps, key=abs)


def oracle_chunkGaps(chunks):
    """
    Signed gaps between a species' consecutive contributing chunks, in
    ascending-start order. Shared by classifyChunks' caller and the per-element
    summary so both see the same definition.
    """
    contributing = sorted((c for c in chunks if c['size'] > 0), key=lambda c: c['start'])
    if len(contributing) < 2:
        return []
    if len({c['scaffold'] for c in contributing}) > 1 or len({c['strand'] for c in contributing}) > 1:
        return []
    return [
        contributing[i + 1]['start'] - (contributing[i]['start'] + contributing[i]['size'])
        for i in range(len(contributing) - 1)
    ]


def oracle_summarizeElementLoci(per_species_chunks):
    """
    Per-ELEMENT summary of how each species' contribution looks, as the values
    for ELEMENT_LOCI_HEADERS.

    Classification is irreducibly per-species -- one real element was
    contiguous for 5 species, split for 6 and empty for 3 -- so this reports
    counts per class plus the extent and the concentration of the splits,
    leaving every threshold to the caller.

    Extent needs two numbers, not one. A real element had a median gap of 63bp
    across 13 split species but one species off by 335,617,958bp: the max alone
    condemns an element where 12 species are nearly fine, the median alone
    waves through a catastrophic misalignment.
      - split.bases.max  : max over species of total unaligned bases, i.e.
                           exactly (header span - aligned bases) for the worst
                           species -- how much sequence its row is missing.
      - split.gap.median : median of every individual gap in the element.

    Concentration answers "is one shared event responsible?". Splits occur at
    block boundaries, which are single reference positions, so species that
    jump at the SAME boundary are jumping at the same place. Comparing
    split.max.boundary.n.species to n.split distinguishes one ancestral indel
    (equal, with a small gap spread -- a real 33/34/33/33/33 case gives spread
    1) from independent lineage-specific events (much lower).

    A species absent from a block inside its own span (1.6% of real multi-chunk
    cases) has a gap spanning more than one boundary; it is attributed to the
    first boundary it crosses, and the per-chunk loci table carries
    block_index/gap_to_next so that case stays inspectable.
    """
    counts = {"single": 0, "contiguous": 0, "split": 0,
              "multi_scaffold": 0, "multi_strand": 0, "no_bases": 0}
    all_gaps = []
    per_species_total = []
    # boundary -> {species: gap}, keyed by the block ordinal the gap starts at
    by_boundary = defaultdict(dict)

    for species, chunks in per_species_chunks.items():
        locus_class, _max_gap = oracle_classifyChunks(chunks)
        if locus_class in counts:
            counts[locus_class] += 1

        gaps = oracle_chunkGaps(chunks)
        nonzero = [g for g in gaps if g != 0]
        if nonzero:
            all_gaps.extend(nonzero)
            per_species_total.append(sum(nonzero))
            contributing = sorted((c for c in chunks if c['size'] > 0), key=lambda c: c['start'])
            for i, gap in enumerate(gaps):
                if gap != 0:
                    by_boundary[contributing[i].get('block_index', i + 1)][species] = gap

    if by_boundary:
        busiest = max(by_boundary.values(), key=len)
        boundary_n = len(busiest)
        spread = max(busiest.values()) - min(busiest.values())
    else:
        boundary_n = 0
        spread = 0

    return [
        len(per_species_chunks),
        counts["single"], counts["contiguous"], counts["split"],
        counts["multi_scaffold"], counts["multi_strand"], counts["no_bases"],
        max(per_species_total) if per_species_total else 0,
        int(statistics.median(all_gaps)) if all_gaps else 0,
        len(by_boundary),
        boundary_n,
        spread,
    ]


# ---------------------------------------------------------------------------
# Random elements that hit every branch
# ---------------------------------------------------------------------------


def random_element(rng):
    """
    {species: [chunk, ...]} with start ties, out-of-order and overlapping
    chunks, zero sizes, several scaffolds and strands (in either order),
    repeated and missing block_index, and species with no chunks at all.
    """
    element = {}
    for s in range(rng.randint(1, 8)):
        chunks = []
        kind = rng.random()
        scaffolds = ["s1"] if kind < 0.6 else ["s1", "s2"]
        strands = ["+"] if kind < 0.4 or kind > 0.8 else ["+", "-"]
        pos = rng.randint(0, 50)
        for _ in range(rng.randint(0, 10)):
            if rng.random() < 0.5:                       # mostly collinear...
                start = pos + rng.choice([0, 0, 0, rng.randint(-5, 40)])
            else:                                        # ...but not always
                start = rng.randint(0, 200)
            size = 0 if rng.random() < 0.15 else rng.randint(1, 30)
            pos = start + size
            c = {"scaffold": rng.choice(scaffolds), "start": start, "size": size,
                 "end": start + size, "strand": rng.choice(strands), "src_size": 1000}
            if rng.random() > 0.05:
                c["block_index"] = rng.randint(1, 6)
            chunks.append(c)
        element[f"sp{s}"] = chunks
    return element


@pytest.mark.parametrize("seed", range(20))
def test_matches_v081_on_random_elements(seed):
    rng = random.Random(seed)
    for _ in range(500):
        element = random_element(rng)
        expected = oracle_summarizeElementLoci(element)
        assert summarizeElementLoci(element) == expected, element


def test_matches_v081_when_fed_chunk_by_chunk_in_block_order():
    """fetch adds chunks block by block, interleaving species -- not per species."""
    rng = random.Random(1234)
    for _ in range(2000):
        element = random_element(rng)
        events = [(c.get("block_index", 0), sp, i, c)
                  for sp, chunks in element.items() for i, c in enumerate(chunks)]
        events.sort(key=lambda e: (e[0], e[2]))
        # the order fetch would see them in, and the per-species dicts in that order
        acc = ElementLociAccumulator()
        reordered = defaultdict(list)
        for _block, sp, _i, c in events:
            reordered[sp].append(c)
            if c["size"] > 0:
                acc.add(sp, c["scaffold"], c["strand"], c["start"], c["size"], c.get("block_index"))
            else:
                acc.add(sp, None, None, None, c["size"], None)
        assert acc.summary() == oracle_summarizeElementLoci(dict(reordered)), dict(reordered)


def test_multi_species_drop_their_chunk_arrays():
    """The memory fix: no per-chunk storage once gaps can no longer matter."""
    acc = ElementLociAccumulator()
    for i in range(1000):
        acc.add("collinear", "s1", "+", i * 10, 10, i)
        acc.add("rearranged", "s1" if i < 500 else "s2", "+", i * 10, 10, i)
        acc.add("inverted", "s1", "+" if i < 500 else "-", i * 10, 10, i)
    assert len(acc.species["collinear"][6]) == 1000
    assert acc.species["rearranged"][6] is None
    assert acc.species["inverted"][6] is None
    vals = dict(zip(ELEMENT_LOCI_HEADERS, acc.summary()))
    assert (vals["n.contiguous"], vals["n.multi.scaffold"], vals["n.multi.strand"]) == (1, 1, 1)


def test_scaffold_change_after_strand_change_is_multi_scaffold():
    """multi_scaffold outranks multi_strand even after the arrays are gone."""
    acc = ElementLociAccumulator()
    acc.add("sp", "s1", "+", 0, 10, 1)
    acc.add("sp", "s1", "-", 20, 10, 2)        # multi_strand: arrays dropped
    acc.add("sp", "s2", "-", 40, 10, 3)        # ...then a second scaffold
    element = {"sp": [
        {"scaffold": "s1", "strand": "+", "start": 0, "size": 10, "block_index": 1},
        {"scaffold": "s1", "strand": "-", "start": 20, "size": 10, "block_index": 2},
        {"scaffold": "s2", "strand": "-", "start": 40, "size": 10, "block_index": 3},
    ]}
    assert acc.summary() == oracle_summarizeElementLoci(element)
    assert dict(zip(ELEMENT_LOCI_HEADERS, acc.summary()))["n.multi.scaffold"] == 1
