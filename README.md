# mafutils

`mafutils` is a command-line toolkit for indexing, extracting, and summarizing
MAF (Multiple Alignment Format) files.

It currently provides six commands:

- `mafutils index`
- `mafutils info`
- `mafutils fetch`
- `mafutils stats`
- `mafutils gc`
- `mafutils validate`

## Disclaimer

This project was developed with significant assistance from large language models (GPT-5 / Codex, Claude Sonnet).

## Third-Party Code

`mafutils/lib/bgzf.py` is vendored from [Biopython](https://biopython.org/)
(BSD 3-Clause License option) rather than pulled in as a dependency — see the
notice at the top of that file and [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)
for full attribution and license text.

## Installation

```bash
pip install mafutils
```

This installs the `mafutils` command.

## Quick Start

Show top-level help:

```bash
mafutils --help
```

Show the installed version:

```bash
mafutils --version    # also accepts -version, -v, -V
```

Create block and scaffold indexes for a MAF (defaults to
`input.maf.block.idx` / `input.maf.scaffold.idx` if output paths are omitted):

```bash
mafutils index input.maf
```

Print a quick overview of an indexed MAF (sub-second on any size):

```bash
mafutils info input.maf
```

Fetch trimmed MAF regions from a BED file (index defaults to
`input.maf.block.idx`):

```bash
mafutils fetch input.maf regions.bed -o outdir
```

Fetch FASTA output instead of MAF:

```bash
mafutils fetch input.maf regions.bed -o outdir -f -fh species-coords-id
```

Extract full scaffolds using a scaffold index:

```bash
mafutils fetch input.maf scaffolds.bed -m scaffold -o outdir
```

Summarize an indexed MAF:

```bash
mafutils stats input.maf -o summary/example
```

Calculate per-species GC content:

```bash
mafutils gc input.maf -o summary/example
```

Check whether an index is still trustworthy (see Index Integrity below):

```bash
mafutils validate input.maf
```

## Compression

`mafutils` supports three MAF compression types, detected automatically:

| Compression | Random access | Notes |
|---|---|---|
| none (plain `.maf`) | Real, parallelizable | Baseline behavior |
| bgzip (BGZF, e.g. produced by `bgzip`) | Real, parallelizable | Detected distinctly from plain gzip via its `BC` extra-field marker; uses virtual offsets for true random access |
| gzip (plain `.gz`) | Not available | Python's `gzip` module can't jump to an arbitrary offset without decompressing everything before it, so `fetch`/`stats`/`gc` all fall back to a single sequential pass regardless of `--processes`, with a warning |

If you plan to use `--processes > 1` against compressed input, compress with
`bgzip` rather than plain `gzip` to get real parallel speedup.

**No extra dependency for bgzip:** BGZF support is vendored in
`mafutils/lib/bgzf.py` rather than depending on the `biopython` package —
see `DEVELOPMENT.md` for why, and [`THIRD_PARTY_LICENSES.md`](THIRD_PARTY_LICENSES.md)
for attribution.

## Index Integrity

`mafutils index` writes a one-line `#`-prefixed header at the top of every
index recording the MAF filename, compression, file size, modification
time, and a content hash (MD5) of the exact file it was built from
(backward compatible — parsers already skip `#` lines; computed for free,
in the same read already used to find block offsets). Index paths default
to `<MAF_FILE>.block.idx` / `<MAF_FILE>.scaffold.idx` next to the MAF file,
but can always be overridden explicitly.

An index only ever means anything paired with the exact file it was built
from — if that file gets regenerated (e.g. a rerun alignment pipeline) with
the same name and compression, an index built against the old version would
otherwise look valid while pointing at stale offsets. `fetch`/`stats`/`gc`
check for this every run, at two tiers:

- **Default** (no flag, cheap — no file read): compression mismatch is an
  **error**; file size mismatch is an **error** (a strong signal — files
  differing in size are definitely not the same content); modification
  time mismatch is only a **warning** (a weak signal — files get
  touched/copied without content changing). A completely missing header
  (an index built by an older `mafutils`) only warns, since there's nothing
  to check.
- **`--verify-hash`** (opt-in, on `fetch`/`stats`/`gc`): compares against
  the stored content hash instead — authoritative, but requires reading
  the whole file to recompute it, so it's not the default.

For a "check once, trust thereafter" workflow instead of passing
`--verify-hash` on every invocation, use `mafutils validate` (see below).

## Progress Reporting

`index`, `stats`, `gc` and `fetch` log a progress line while they run, every
60 seconds by default. Each takes `--progress-interval SECONDS` to change it:

- `60` or more sets the interval;
- `0` turns progress off entirely: no lines, and none of the bookkeeping
  behind them runs;
- anything from 1 to 59, or negative, is rejected with a usage error rather
  than adjusted.

The option belongs to each command, so it goes after the command name:
`mafutils stats input.maf --progress-interval 300`. What the lines count
(blocks, bytes, regions) differs by command; see each command below.

## Commands

### `mafutils index`

Create block and scaffold indexes for a MAF file.

```bash
mafutils index MAF_FILE [BLOCK_INDEX] [SCAFFOLD_INDEX]
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf`, `.maf.gz`, or bgzip-compressed `.maf`) |
| `BLOCK_INDEX` | Output block index path (default: `<MAF_FILE>.block.idx`) |
| `SCAFFOLD_INDEX` | Output scaffold index path (default: `<MAF_FILE>.scaffold.idx`) |

Options:

| Option | Description |
|---|---|
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns progress off. See Progress Reporting above |

Every 60 seconds by default it logs how far through the file it is, e.g. `Progress:
12.30 GB / 42.70 GB (28.8%), 123.00 MB/s, 1m 40s elapsed, ~4m 07s left`.
Progress is in bytes because the block count is what indexing produces. For
gzip and bgzip input both numbers are compressed bytes, compared against the
compressed file size. Messages go to the screen only; no log file is written.

### `mafutils info`

Print a quick overview of a MAF without reading it.

```bash
mafutils info [OPTIONS] MAF_FILE
```

```
MAF            : cricetid-15spec.Mus_musculus.nodupes.maf
size           : 42.70 GB (none)
index          : matches MAF (size + mtime)

ref scaffolds  : 61
blocks         : 8,094,203
ref bases      : 2,728,222,451
aln columns    : 3,031,569,326
seq lines      : 103,950,448   (all species, including the reference)
max seqs/block : 15

species        : 15  (sampled from the first 1,000 block(s) -- possibly incomplete)
  Arvicola_amphibius, Cricetulus_griseus, Cricetus_cricetus, ...

For exact per-species presence across the whole file, run `mafutils stats`.
```

Measured on that 42.70 GB file: **1.3 s**, or **0.7 s** with `--sample-blocks 0`.

| Option | Description |
|---|---|
| `--index`, `-i` | Block index (default: `<MAF_FILE>.block.idx`) |
| `--scaffold-index` | Scaffold index (default: `<MAF_FILE>.scaffold.idx`) |
| `--sample-blocks` | Blocks to sample for species names (default 1000); `0` disables the sample and reads no MAF data at all |
| `--all-species` | List every sampled species instead of the first 20 |

**Where the numbers come from.** `mafutils index` records the block counts in
the index header, so `info` normally reads nothing but that one line. Three
tiers, and the MAF is never scanned for these counts:

1. Index built by mafutils 0.7 or later — instant, header only.
2. Older index — falls back to streaming the *index* (~52 MB/s, so ~9 s for a
   445 MB index but ~12 min for a 35 GB one), with a warning naming the index
   size and recommending a rebuild.
3. No index — reports what needs none and points at `mafutils index`.

**`species` is a sample, not a count.** Species names are not recorded in the
index, and collecting them during indexing would mean parsing every
non-reference `s` line — measured at 2.06x slower indexing. Sampling the first
1,000 blocks costs milliseconds instead. On a whole-genome MAF those blocks all
sit at the start of the first scaffold, so a species absent from that region is
missed. Use `mafutils stats` for an exhaustive per-species breakdown,
including each species' source scaffolds and their lengths.

**Column meanings.** `ref scaffolds` and `ref bases` are measured on the
*reference* sequence. `seq lines` counts every `s` line including the
reference. `aln columns` is the block width, gap columns included, which all
species share.

### `mafutils fetch`

Fetch regions or scaffolds from a MAF using an existing index.

```bash
mafutils fetch [OPTIONS] MAF_FILE BED_FILE
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf`, `.maf.gz`, or bgzip-compressed `.maf`) |
| `BED_FILE` | BED file with regions or scaffold names |

Options:

| Option | Description |
|---|---|
| `--index`, `-i` | Index file, block-level or scaffold-level matching `--mode`. **Required** — if omitted, looked up at `<MAF_FILE>.block.idx` or `<MAF_FILE>.scaffold.idx` by default; errors if not found there. |
| `--basename`, `-b` | Output basename strategy: `id` (the BED's 4th column; every region needs a unique one -- see *BED pre-check* below), `coords`, or `count` |
| `--output`, `-o` | Output directory or output filename in single-output mode |
| `--fasta`, `-f` | Write FASTA instead of MAF |
| `--fasta-header`, `-fh` | FASTA header format: `species-coords-id`, `species-coords`, or `species-only`. The coords-bearing modes always include the source scaffold (`species.scaffold`); `species-only` is deliberately bare so headers match tree tip labels |
| `--expected-species` | Comma-separated expected species list for FASTA filling |
| `--expected-species-file` | File with one expected species name per line |
| `--fasta-dedupe` | FASTA duplicate handling: `none` or `most-seq` |
| `--processes`, `-p` | Number of worker processes (see Compression above — plain gzip always runs single-process) |
| `--mode`, `-m` | Fetch mode: `block` or `scaffold` |
| `--scaffold-subdirs` | Group output files into subfolders named by reference scaffold (`<output>/<scaffold>/<basename>`) instead of one flat directory |
| `--loci-table` | Also write `<output>/maf_fetch_loci.tsv`: one row per region/species/contributing-block with each species' source scaffold, coordinates, `srcSize`, and a locus classification (see *Per-species coordinates* below). Requires `--fasta` |
| `--verbose` | Emit warning lines from each completed batch |
| `--profile` | Log internal timing breakdowns |
| `--verify-hash` | Verify the index's stored content hash against the MAF file (see Index Integrity above) |
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns progress off. See Progress Reporting above |

**BED pre-check.** In block mode, before doing any work, `fetch` checks the
whole BED and stops with a single error listing every problem it finds:

- lines it cannot use: fewer than 3 columns, non-integer start or end, a
  negative start, or end ≤ start (each reported with its line number);
- with `-b id`, regions that have no name (4th column);
- regions that would write the same output file, where the later one would
  silently replace the earlier: names used more than once with `-b id`
  (allowed on different scaffolds with `--scaffold-subdirs`), or the same
  region listed twice with `-b coords`. `-b count` numbers every region, so it
  never collides.

Blank lines, `#` comments, and UCSC `track`/`browser` header lines are
skipped. Scaffolds that are not in the MAF are not an error; they are reported
as a warning at the end, as before. The check takes under a second on a
979,442-region BED.

**Progress.** `fetch` logs one progress line every 60 seconds by default
(`--progress-interval`), combining all workers:

- **Block mode** counts regions, e.g. `Progress: 412,000 / 979,442 regions
  (42.1%), 1,850 regions/s, 3m 42s elapsed, ~5m 06s left`. Workers finish
  their shares at different times (the "Completed batch N" lines show which),
  so the combined rate drops near the end and the last estimates run low.
- **Scaffold mode** counts bytes copied against the total size of the requested
  scaffolds, and logs `Finished <scaffold> (k of M scaffolds)` as each one
  completes. For bgzip input it counts finished scaffolds instead of bytes:
  its copy loop runs once per line and is kept free of progress checks.

Each worker records its count in its own slot of shared memory, and the main
process reads them only when a line is due. No messages pass between processes
and no worker ever waits on another, so this adds no measurable time at any
scale.

### Reference-side columns in `maf_fetch_summary.tsv`

Three of the original columns describe the **reference**, not any species, and
the `ref.` prefix marks them as such:

| column | meaning |
|---|---|
| `ref.n.overlapping.blocks` | how many alignment blocks the region overlaps |
| `ref.block.bases` | reference bases in the region that alignment blocks cover |

`ref.block.bases` normally equals `end - start`. It falls short exactly when
the region contains reference positions no block covers, and in that case
`fetch` emits a `ref-coverage-gap` warning naming the distances involved.

This distinction matters: **the reference is contiguous by construction, so
these columns say nothing about whether a species' sequence is contiguous.** A
species can jump kilobases at a boundary where the reference does not. The old
`interblock.distances` column was removed for exactly this reason — it
reported reference-side gaps, which are all zeros on a normal MAF (measured:
19,999/19,999 consecutive block pairs adjacent), and reading it as a
contiguity check was misleading. Per-species contiguity lives in the `n.*` /
`split.*` columns below and in the loci table's `gap_to_next`.

### Per-element locus classification (always written)

`maf_fetch_summary.tsv` carries one row per BED element and reports how each
species' contribution to that element looks. These columns are always written
(no flag needed) and work with MAF or FASTA output:

| column | meaning |
|---|---|
| `n.species` | species contributing to this element |
| `n.single`, `n.contiguous` | species whose contribution **is** one locus — for these, the header span equals the bases emitted |
| `n.split` | species interrupted by sequence the aligner declined to align |
| `n.multi.scaffold`, `n.multi.strand` | species stitched from different source sequences or strands |
| `n.no.bases` | species present in the alignment but contributing zero bases to this element (an all-gap row) |
| `split.bases.max` | worst species' total unaligned bases, i.e. how much sequence its row is missing |
| `split.gap.median` | median of every individual gap in the element |
| `split.boundaries` | distinct block boundaries with any split |
| `split.max.boundary.n.species` | most species splitting at a single boundary |
| `split.max.boundary.gap.spread` | `max - min` of those species' gaps at that same boundary |

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

### Per-species coordinates (`--loci-table`)

FASTA headers name the source scaffold, but for a **non-reference** species a
header's `start-end` is a *bounding box* over every block that contributed, not
necessarily one locus. MAF blocks are contiguous in the reference only: a
species' chunks can be separated, inverted, or on different scaffolds, and
stitching concatenates them without padding. On real data (2,000 regions, 15
species) 88.5% of region/species pairs are unambiguously one locus, but ~11%
are split and ~0.8% span more than one scaffold.

`--loci-table` writes the per-chunk truth so you can filter before treating a
row as a locus:

```bash
mafutils fetch input.maf regions.bed -f -fh species-coords-id --loci-table -o out/
```

`out/maf_fetch_loci.tsv` has one row per (region, species, contributing block):

| column | meaning |
|---|---|
| `region_scaffold`, `region_start`, `region_end`, `region_id` | the requested reference region |
| `species`, `src_scaffold`, `src_size` | that species' source sequence for this chunk |
| `chunk_start`, `chunk_size`, `chunk_strand` | the chunk, exactly as the MAF states it |
| `chunk_index`, `n_chunks` | position among that species' chunks for this region |
| `block_index` | which block (ordinal within the element) this chunk came from |
| `gap_to_next` | signed gap to that species' next chunk; `.` for the last |
| `class` | `single`, `contiguous`, `split`, `multi_strand`, `multi_scaffold`, or `no_bases` |
| `max_gap` | largest signed gap between consecutive chunks; `0` for single/contiguous, `.` where undefined |

`class` is constant within a (region, species) group, so filtering is a one-liner:

```bash
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)c[$i]=i; next}
  $c["class"]=="single" || $c["class"]=="contiguous"' out/maf_fetch_loci.tsv
```

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

The table is roughly (regions x species) rows, so it is off by default — a
979k-region BED at 15 species produces ~13.7M rows.

On plain-gzip input, `fetch` precomputes the set of blocks needed across
*all* regions, decodes each exactly once in strictly-ascending file order,
then serves every region from that cache — this avoids the backward seeks
that would otherwise come from processing regions in arbitrary (e.g. BED
file) order, and decodes shared blocks only once even when regions overlap.

### `mafutils stats`

Summarize an indexed MAF at overall, species, and block levels.

```bash
mafutils stats [OPTIONS] MAF_FILE [INDEX_FILE]
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf`, `.maf.gz`, or bgzip-compressed `.maf`) |
| `INDEX_FILE` | Block index produced by `mafutils index`. **Required** — if omitted, looked up at `<MAF_FILE>.block.idx` by default; errors if not found there. |

Options:

| Option | Description |
|---|---|
| `--output-prefix`, `-o` | Output prefix/path |
| `--processes`, `-p` | Number of worker processes |
| `--chunk-size` | Blocks per worker task |
| `--no-block-table` | Skip writing the per-block table |
| `--expected-species` | Comma-separated species names for exact missing lists |
| `--expected-species-file` | File with one species name per line |
| `--verify-hash` | Verify the index's stored content hash against the MAF file (see Index Integrity above) |
| `--log-level` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns progress off. See Progress Reporting above |
| `--html-dashboard` | Write an HTML summary dashboard |
| `--dashboard-top-species` | Number of species shown in dashboard bar plots |
| `--dashboard-max-block-points` | Maximum block rows sampled for dashboard plots |
| `--dashboard-top-scaffolds` | Largest scaffolds listed per species in the dashboard (default 10). The TSV is never capped. |

Outputs:

| File | Contents |
|---|---|
| `<prefix>.overall.tsv` | One metric per line |
| `<prefix>.species.tsv` | One row per species |
| `<prefix>.species_scaffolds.tsv` | One row per (species, source scaffold) present in the MAF — see below |
| `<prefix>.block.tsv` | One row per block (skip with `--no-block-table`) |
| `<prefix>.dashboard.html` | Self-contained summary dashboard (with `--html-dashboard`) |

**Progress.** While it runs, `stats` logs a progress line every 60 seconds by
default (`--progress-interval`), e.g. `Progress: 1,215,000 / 8,094,203 blocks (15.0%), 2,850 blocks/s, 7m 06s
elapsed, ~40m 13s left`. The total comes from the index header, so an index
built before mafutils 0.7 gives counts and rate only (rebuild it to get the
percentage). Block sizes vary across a genome, so the time estimate is rough.
If the number of blocks processed differs from the header's count, `stats`
warns that the index rows and header disagree and the results may not cover
the whole MAF. `gc` reports progress the same way; `index` measures it in
bytes instead (see `mafutils index`).

### Per-species source scaffolds

`stats` reads every `s` line of every block, so it also records which source
scaffolds each species contributes and how long they are. Nothing else reports
this for non-reference species: the index never parses non-reference lines, and
`info` only samples species names.

**Scope: scaffolds present in the MAF, not the whole assembly.** A scaffold
appears only if at least one of its bases aligned to the reference. Unaligned
scaffolds — unplaced contigs, sequence outside the aligned regions — never
appear in the MAF, so `stats` cannot see them. Each row's length is exact (it is
the MAF's `srcSize`, which the aligner takes from the assembly), but the *set*
of scaffolds is limited to those in the MAF. The per-species totals therefore
err in opposite directions depending on the comparison:

| Compared with | Direction | Why |
|---|---|---|
| the true assembly size | **under**counts | unaligned scaffolds contribute nothing |
| aligned bases (`nongaps_total`) | **over**counts | a barely-aligned 200 Mb scaffold contributes all 200 Mb |

The gap to the true assembly cannot be measured from the MAF alone — it needs
the assembly's own scaffold list. It is small for chromosome-level assemblies
and can be large for fragmented ones. Do not read `scaffold_bases_in_maf` as an
assembly-size estimate.

`<prefix>.species_scaffolds.tsv`, sorted by species, then largest scaffold
first:

| Column | Meaning |
|---|---|
| `species` | Species (the part of the MAF `src` before the first `.`) |
| `scaffold` | Source scaffold (the part after it); `.` if the `src` has no `.` |
| `src_size` | Full length of that scaffold in that species' assembly (MAF `srcSize`); `.` if it did not parse as an integer |
| `blocks` | Blocks in which this scaffold appears |
| `aligned_bases` | Non-gap bases from this scaffold across all blocks |
| `aligned_fraction` | `aligned_bases / src_size`. **Not capped at 1.0**: a region aligned into several blocks (duplicated or paralogous alignments) can push it past 1, and that is worth seeing |

Rows are keyed by (species, scaffold), never by scaffold name alone —
`human.chr4` and `chimp.chr4` are different sequences that happen to share a
name.

Two columns are appended to the **end** of `<prefix>.species.tsv`, so no
existing column moves:

| Column | Meaning |
|---|---|
| `scaffolds_in_maf` | Distinct source scaffolds of this species present in the MAF |
| `scaffold_bases_in_maf` | Sum of those scaffolds' full lengths |

Per species, the `aligned_bases` of its scaffolds sum exactly to its
`nongaps_total`, and their `blocks` sum to at least its `blocks_present`
(more when a block holds several copies of the species).

The dashboard adds both columns to its species table and a collapsed section
per species listing the largest scaffolds. Each section's summary line reads
"showing N of M scaffolds", so a truncated list is visible before expanding it.

**`srcSize` integrity.** If one species.scaffold reports two different
`srcSize` values, the MAF is internally inconsistent and `stats` exits with an
error naming the scaffold and both values — picking one would hide the problem.
An `s` line whose `srcSize` is not an integer still counts toward every other
statistic. Its scaffold takes its length from any other line of that scaffold
that does parse, and is written with `src_size` `.` only if none does. The
number of such lines is logged as a warning and written to `overall.tsv` as
`scaffold_size_parse_errors`.

### `mafutils gc`

Calculate per-species GC content from a MAF file. Writes `<prefix>.gc.csv`
(one `species,gc` row per species, GC as a 0-1 fraction) and
`<prefix>.gc.mean.txt` (a single float: the unweighted mean of the `gc`
column across species).

GC is computed as `(G+C) / (A+C+G+T)`, case-insensitive; gaps, `N`, and other
ambiguity codes are excluded from both the numerator and denominator.

```bash
mafutils gc MAF_FILE [INDEX_FILE] [OPTIONS]
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf`, `.maf.gz`, or bgzip-compressed `.maf`) |
| `INDEX_FILE` | Block index produced by `mafutils index`. **Not required** — `gc` works fine without one. If given, or found at `<MAF_FILE>.block.idx` by default, it can speed things up via `--processes > 1` on uncompressed/bgzip MAFs. |

Options:

| Option | Description |
|---|---|
| `--output-prefix`, `-o` | Output prefix/path |
| `--processes`, `-p` | Number of worker processes (see note below) |
| `--chunk-size` | Blocks per worker task in parallel mode |
| `--verify-hash` | Verify the index's stored content hash against the MAF file (see Index Integrity above; only checked when an index is actually used) |
| `--log-level` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR` |
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns progress off. See Progress Reporting above |

**Sequential vs. parallel processing:** GC counting is cheap per base, so the
bottleneck is I/O/decompression rather than CPU. There are two processing
paths:

- **Sequential** (default): streams through the MAF once, top to bottom, with
  no index required. Used whenever no index is found, the MAF is plain-gzip
  compressed, or `--processes` is `1`.
- **Parallel**: requires both an index and `--processes > 1`, and is used on
  **uncompressed or bgzip-compressed** MAFs. It chunks index entries across
  worker processes that each seek to their own byte ranges (real random
  access for both of these types), the same architecture as `mafutils stats`.

Parallel mode is skipped for plain-gzip input even if requested: the index's
byte offsets are positions in the *decompressed* stream, so each worker would
have to redundantly re-decompress everything up to its own start point,
making things slower rather than faster. A warning is logged and the tool
falls back to the sequential path instead of silently ignoring the request —
use a bgzip-compressed MAF instead if you need real parallel speedup on
compressed input.

### `mafutils validate`

Checks a MAF file against its block index's stored header (size, mtime, and
content hash) and reports whether the index can still be trusted — a
"check once, trust thereafter" alternative to passing `--verify-hash` on
every `fetch`/`stats`/`gc` invocation. Also looks up the scaffold index at
its default location (`<MAF_FILE>.scaffold.idx`) and cross-checks its
header against the block index's — `mafutils index` always writes both
from the same size/mtime/hash in the same run, so any difference means the
two indexes weren't built together (e.g. only one was rebuilt) and
shouldn't be trusted as a matched pair.

```bash
mafutils validate MAF_FILE [INDEX_FILE]
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf`, `.maf.gz`, or bgzip-compressed `.maf`) |
| `INDEX_FILE` | Block index produced by `mafutils index`. **Required** — if omitted, looked up at `<MAF_FILE>.block.idx` by default; errors if not found there. The scaffold index is always looked up at `<MAF_FILE>.scaffold.idx` (no separate argument for it). |

Prints one line per checked field (compression, size, hash, mtime, and the
block/scaffold header pair) with its match/mismatch status, then an overall
verdict. Exit codes:

| Exit code | Verdict | Meaning |
|---|---|---|
| `0` | VERIFIED | The index's stored hash matches this file exactly, and the block/scaffold index headers agree. |
| `1` | MISMATCH | A conclusive difference was found (compression, size, or hash against the MAF file; or the block and scaffold index headers disagree) — rebuild the index(es). |
| `2` | UNVERIFIABLE | Nothing contradicts, but there's no stored hash to be fully sure (the index predates this feature), or the scaffold index is missing/headerless so the pair couldn't be cross-checked. |

## Notes

For development setup (installing from a source checkout), running the test
suite, and internal implementation notes, see [`DEVELOPMENT.md`](DEVELOPMENT.md).
