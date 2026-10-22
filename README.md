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
notice at the top of that file and [`THIRD_PARTY_LICENSES.md`](https://github.com/harvardinformatics/mafutils/blob/main/THIRD_PARTY_LICENSES.md)
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
see `DEVELOPMENT.md` for why, and [`THIRD_PARTY_LICENSES.md`](https://github.com/harvardinformatics/mafutils/blob/main/THIRD_PARTY_LICENSES.md)
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
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns them off |

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
| `--basename`, `-b` | Output basename strategy: `id` (the BED's 4th column, which must be present and unique), `coords`, or `count` |
| `--output`, `-o` | Output directory or output filename in single-output mode |
| `--fasta`, `-f` | Write FASTA instead of MAF |
| `--fasta-header`, `-fh` | FASTA header format: `species-coords-id`, `species-coords`, or `species-only`. The coords-bearing modes always include the source scaffold (`species.scaffold`); `species-only` is deliberately bare so headers match tree tip labels |
| `--expected-species` | Comma-separated expected species list for FASTA filling |
| `--expected-species-file` | File with one expected species name per line |
| `--fasta-dedupe` | FASTA duplicate handling: `none` or `most-seq` |
| `--processes`, `-p` | Number of worker processes (see Compression above — plain gzip always runs single-process) |
| `--mode`, `-m` | Fetch mode: `block` or `scaffold` |
| `--scaffold-subdirs` | Group output files into subfolders named by reference scaffold (`<output>/<scaffold>/<basename>`) instead of one flat directory |
| `--loci-table` | Also write `<output>/maf_fetch_loci.tsv` (see below). Requires `--fasta` |
| `--verbose` | Emit warning lines from each completed batch |
| `--profile` | Log internal timing breakdowns |
| `--verify-hash` | Verify the index's stored content hash against the MAF file (see Index Integrity above) |
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns them off |

**BED pre-check.** In block mode, `fetch` checks the whole BED before doing
any work and stops with one error listing every problem: unusable lines
(malformed, or end ≤ start), regions without a name under `-b id`, and regions
that would overwrite each other's output file (a repeated name under `-b id`,
or a repeated region under `-b coords`). `#` comments and `track`/`browser`
lines are skipped.

**Outputs:**

- One file per region (`.maf`, or `.fa` with `--fasta`), named by `--basename`,
  in `--output` (or in per-scaffold subfolders with `--scaffold-subdirs`).
  Scaffold mode writes one `<scaffold>.maf` per scaffold instead.
- `maf_fetch_summary.tsv`: one row per BED region, in genomic order.
- `maf_fetch_loci.tsv` (with `--loci-table`): one row per region, species and
  contributing block.

Both tables describe how each species' sequence for a region is put together.
A region's alignment blocks are contiguous in the reference, but a species'
pieces of those blocks ("chunks") need not be contiguous in its own genome, so
each species is given a class:

| class | meaning |
|---|---|
| `single` | one chunk |
| `contiguous` | several chunks, adjacent in the species' own sequence |
| `split` | chunks separated by sequence the aligner did not align |
| `multi_strand` | chunks on both strands |
| `multi_scaffold` | chunks from more than one source scaffold |
| `no_bases` | present in the alignment, but all gaps in this region |

Only `single` and `contiguous` are one locus: for those, the FASTA header's
span equals the bases written. For the others it never does; the span is a
bounding box over chunks that are separated, on both strands, or on different
scaffolds.

For background on these classes, filtering recipes, and an example of building
a per-species BED, see [docs/fetch-outputs.md](https://github.com/harvardinformatics/mafutils/blob/main/docs/fetch-outputs.md).

### `maf_fetch_summary.tsv`

| column | meaning |
|---|---|
| `scaffold`, `start`, `end`, `basename` | the BED region and its output file name |
| `ref.n.overlapping.blocks` | alignment blocks the region overlaps |
| `ref.block.bases` | reference bases in the region covered by blocks; less than `end - start` when part of the region is unaligned (logged as a `ref-coverage-gap` warning) |
| `n.sequences` | sequences in the region's FASTA file, including `--expected-species` filled with Ns; `NA` for MAF output |
| `n.species` | species contributing to the region |
| `n.single`, `n.contiguous`, `n.split`, `n.multi.strand`, `n.multi.scaffold`, `n.no.bases` | species in each class |
| `split.bases.max` | the most unaligned bases in any one split species |
| `split.gap.median` | median of every individual gap between chunks |
| `split.boundaries` | block boundaries where at least one species splits |
| `split.max.boundary.n.species` | the most species splitting at a single boundary |
| `split.max.boundary.gap.spread` | `max - min` of those species' gaps at that boundary |

The `ref.` columns describe the reference only. They say nothing about whether
a species' sequence is contiguous; the `n.*` and `split.*` columns do.

### `maf_fetch_loci.tsv` (`--loci-table`)

```bash
mafutils fetch input.maf regions.bed -f -fh species-coords-id --loci-table -o out/
```

| column | meaning |
|---|---|
| `region_scaffold`, `region_start`, `region_end`, `region_id` | the requested reference region |
| `species`, `src_scaffold`, `src_size` | that species' source sequence for this chunk |
| `chunk_start`, `chunk_size`, `chunk_strand` | the chunk, exactly as the MAF states it |
| `chunk_index`, `n_chunks` | position among that species' chunks for this region |
| `block_index` | which block (ordinal within the region) this chunk came from |
| `gap_to_next` | signed gap to that species' next chunk; `.` for the last |
| `class` | the species' class for this region (table above); the same on all its rows |
| `max_gap` | largest signed gap between consecutive chunks; `0` for single/contiguous, `.` where undefined |

Coordinates are exactly as in the MAF, so on the `-` strand they count from the
end of the source sequence. `src_size` is included so forward coordinates are
one subtraction away: `src_size - (chunk_start + chunk_size)` to
`src_size - chunk_start`. The table has roughly (regions × species) rows, so it
is off by default.

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
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns them off |
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
| `--progress-interval` | Seconds between progress lines (default and minimum 60); `0` turns them off |

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
suite, and internal implementation notes, see [`DEVELOPMENT.md`](https://github.com/harvardinformatics/mafutils/blob/main/DEVELOPMENT.md).
