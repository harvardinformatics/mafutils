# mafutils

`mafutils` is a command-line toolkit for indexing, extracting, and summarizing
MAF (Multiple Alignment Format) files.

It currently provides five commands:

- `mafutils index`
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

Create block and scaffold indexes for a MAF (defaults to
`input.maf.block.idx` / `input.maf.scaffold.idx` if output paths are omitted):

```bash
mafutils index input.maf
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
| `--basename`, `-b` | Output basename strategy: `id`, `coords`, or `count` |
| `--output`, `-o` | Output directory or output filename in single-output mode |
| `--fasta`, `-f` | Write FASTA instead of MAF |
| `--fasta-header`, `-fh` | FASTA header format: `species-coords-id`, `species-coords`, or `species-only` |
| `--expected-species` | Comma-separated expected species list for FASTA filling |
| `--expected-species-file` | File with one expected species name per line |
| `--fasta-dedupe` | FASTA duplicate handling: `none` or `most-seq` |
| `--processes`, `-p` | Number of worker processes (see Compression above — plain gzip always runs single-process) |
| `--mode`, `-m` | Fetch mode: `block` or `scaffold` |
| `--verbose` | Emit warning lines from each completed batch |
| `--profile` | Log internal timing breakdowns |
| `--verify-hash` | Verify the index's stored content hash against the MAF file (see Index Integrity above) |

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
| `--html-dashboard` | Write an HTML summary dashboard |
| `--dashboard-top-species` | Number of species shown in dashboard bar plots |
| `--dashboard-max-block-points` | Maximum block rows sampled for dashboard plots |

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
