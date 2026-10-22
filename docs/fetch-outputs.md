# `mafutils fetch` outputs: background and recipes

The [README](../README.md#mafutils-fetch) describes each `fetch` output and
every column. This page has the background behind the locus classes and the
`split.*` columns, and worked `awk` recipes for filtering elements and building
a per-species BED.

## A FASTA header is not always one locus

FASTA headers name the source scaffold, but for a **non-reference** species a
header's `start-end` is a *bounding box* over every block that contributed, not
necessarily one locus. MAF blocks are contiguous in the reference only: a
species' chunks can be separated, inverted, or on different scaffolds, and
stitching concatenates them without padding. On real data (2,000 regions, 15
species) 88.5% of region/species pairs are unambiguously one locus, but ~11%
are split and ~0.8% span more than one scaffold.

## Coordinates and gaps

**Coordinates are reported exactly as the MAF states them.** For `strand = -`
that means they count in the reverse-complemented source, which is MAF's own
convention. `src_size` is included so forward-genomic coordinates are one
subtraction away — `src_size - (chunk_start + chunk_size)` to
`src_size - chunk_start` — rather than mafutils choosing a frame for you. No
gap threshold is applied either. A `split` gap is **not** an indel the
alignment absorbed: the aligner threads query insertions inline as
reference-row gap columns only up to ~37 bp, and breaks the block past that, so
a gap between chunks is query sequence the aligner *declined to align* — real
bases absent from the FASTA while the header's span still covers them. The
consequence holds with no exceptions: `single`/`contiguous` records have
`header span == emitted bases`, `split`/`multi_*` never do. Real gaps range
from a median of 45 bp to 335 Mb, so `max_gap` is reported and the judgment is
yours.

## The `ref.` columns describe the reference only

This distinction matters: **the reference is contiguous by construction, so
these columns say nothing about whether a species' sequence is contiguous.** A
species can jump kilobases at a boundary where the reference does not. The old
`interblock.distances` column was removed for exactly this reason — it
reported reference-side gaps, which are all zeros on a normal MAF (measured:
19,999/19,999 consecutive block pairs adjacent), and reading it as a
contiguity check was misleading. Per-species contiguity lives in the `n.*` /
`split.*` summary columns and in the loci table's `gap_to_next`.

## Reading the `split.*` columns

Both extent columns matter: one real element had a median gap of 63 bp across
13 split species but one species off by 335 Mb. The max alone would condemn an
element where 12 species are nearly fine; the median alone would wave through a
catastrophic misalignment.

The concentration columns answer *"is one shared event responsible?"*. Splits
occur at block boundaries, which are single reference positions, so species
jumping at the same boundary are interrupted at the same place. Compare
`split.max.boundary.n.species` to `n.split`: equal, with a small
`gap.spread`, means one ancestral indel (a real case had six species at
33/34/33/33/33/33 — spread 1); much lower means independent lineage-specific
events. Where nearly every species jumps at one boundary, the most
parsimonious reading is a reference-specific deletion, which makes the
element's own reference definition the questionable thing.

## Filtering

**Filter at the granularity your product actually needs** — these are very
different yields, and using the strict one for a single-species product
throws away most of your data.

These examples look columns up by **name**, so they keep working if the
schema changes:

*All species* — "safe to claim a locus in every species", e.g. for a
concatenated alignment you want uniformly clean:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; print; next}
  $c["n.split"]==0 && $c["n.multi.scaffold"]==0 &&
  $c["n.multi.strand"]==0 && $c["n.no.bases"]==0' maf_fetch_summary.tsv
```

*One species* — for a per-species product such as a BED in a single
non-reference genome, filter that species' own rows in the loci table
instead; an element unusable for one species is usually fine for the rest.

Two things to get right when aggregating those rows:

1. **The table is one row per contributing block**, so a `contiguous` species
   with two chunks yields two rows for one element. Since `contiguous` means
   the chunks are adjacent by definition, merging them (min start, max end)
   is lossless and gives one interval per element. On one real dataset a
   species had 1,703 clean chunk rows across 1,462 elements — counting rows
   would overstate the element count by 16%.
2. **Skip `chunk_size == 0` rows.** They are emitted so nothing is hidden, but
   they contribute no sequence and are excluded from `class` — so a `single`
   record can still have a zero-size row, *possibly on a different scaffold*.
   Aggregating without filtering them produced a nonsense 5.9 Mb interval in
   testing.

A complete per-species BED, merged and converted to forward-strand
coordinates (verified: 1,462 intervals, median width 109 bp against a
reference element median of 129 bp):

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  $c["species"]=="YOUR_SPECIES" && $c["chunk_size"]+0>0 &&
  ($c["class"]=="single" || $c["class"]=="contiguous") {
    k=$c["region_scaffold"]":"$c["region_start"]
    s=$c["chunk_start"]+0; e=s+$c["chunk_size"]
    if(!(k in lo) || s<lo[k]) lo[k]=s
    if(!(k in hi) || e>hi[k]) hi[k]=e
    scaf[k]=$c["src_scaffold"]; str[k]=$c["chunk_strand"]
    sz[k]=$c["src_size"]+0; id[k]=$c["region_id"]
  }
  END{ for(k in lo){
      if(str[k]=="-"){ a=sz[k]-hi[k]; b=sz[k]-lo[k] } else { a=lo[k]; b=hi[k] }
      print scaf[k]"\t"a"\t"b"\t"id[k]"\t0\t"str[k]
  }}' maf_fetch_loci.tsv > species.bed
```

The `str[k]=="-"` branch is the forward-strand conversion described above —
`src_size - end` .. `src_size - start`, since mafutils reports MAF-frame
coordinates and leaves the choice of frame to you.

On one real 15-species dataset the all-species filter kept **65.8%** of
elements while the per-species filter kept **82-86%** (mean 83.7%) depending
on the species; on a more fragmented 45-way alignment the per-species figure
was reported at **98.6%**. The gap is dataset-specific, so measure yours
rather than assuming either number.

*One-locus rows only* — `class` is constant within a (region, species)
group, so keeping only the one-locus rows of the loci table is a one-liner:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  $c["class"]=="single" || $c["class"]=="contiguous"' out/maf_fetch_loci.tsv
```
