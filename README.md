# mafutils

`mafutils` is a command-line toolkit for indexing, extracting, and summarizing
MAF (Multiple Alignment Format) files.

It currently provides four commands:

- `mafutils index`
- `mafutils fetch`
- `mafutils stats`
- `mafutils gc`

## Disclaimer

This project was developed with significant assistance from a large language model (GPT-5 / Codex).

## Installation

From this `mafutils/` directory:

```bash
pip install -e .
```

If you are using the current shared workflow environment instead of a separate
package environment, you can also run it directly with:

```bash
python -m mafutils --help
```

## Quick Start

Show top-level help:

```bash
python -m mafutils --help
```

Create block and scaffold indexes for a MAF:

```bash
python -m mafutils index input.maf output.block.idx output.scaffold.idx
```

Fetch trimmed MAF regions from a BED file:

```bash
python -m mafutils fetch input.maf input.block.idx regions.bed -o outdir
```

Fetch FASTA output instead of MAF:

```bash
python -m mafutils fetch input.maf input.block.idx regions.bed -o outdir -f -fh species-coords-id
```

Extract full scaffolds using a scaffold index:

```bash
python -m mafutils fetch input.maf input.scaffold.idx scaffolds.bed -m scaffold -o outdir
```

Summarize an indexed MAF:

```bash
python -m mafutils stats input.maf input.block.idx -o summary/example
```

Calculate per-species GC content:

```bash
python -m mafutils gc input.maf -o summary/example
```

## Commands

### `mafutils index`

Create block and scaffold indexes for a MAF file.

```bash
python -m mafutils index MAF_FILE BLOCK_INDEX SCAFFOLD_INDEX
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf` or `.maf.gz`) |
| `BLOCK_INDEX` | Output block index path |
| `SCAFFOLD_INDEX` | Output scaffold index path |

### `mafutils fetch`

Fetch regions or scaffolds from a MAF using an existing index.

```bash
python -m mafutils fetch [OPTIONS] MAF_FILE INDEX_FILE BED_FILE
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf` or `.maf.gz`) |
| `INDEX_FILE` | Block-level or scaffold-level index |
| `BED_FILE` | BED file with regions or scaffold names |

Options:

| Option | Description |
|---|---|
| `--basename`, `-b` | Output basename strategy: `id`, `coords`, or `count` |
| `--output`, `-o` | Output directory or output filename in single-output mode |
| `--fasta`, `-f` | Write FASTA instead of MAF |
| `--fasta-header`, `-fh` | FASTA header format: `species-coords-id`, `species-coords`, or `species-only` |
| `--expected-species` | Comma-separated expected species list for FASTA filling |
| `--expected-species-file` | File with one expected species name per line |
| `--fasta-dedupe` | FASTA duplicate handling: `none` or `most-seq` |
| `--processes`, `-p` | Number of worker processes |
| `--mode`, `-m` | Fetch mode: `block` or `scaffold` |
| `--verbose` | Emit warning lines from each completed batch |
| `--profile` | Log internal timing breakdowns |

### `mafutils stats`

Summarize an indexed MAF at overall, species, and block levels.

```bash
python -m mafutils stats [OPTIONS] MAF_FILE INDEX_FILE
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf` or `.maf.gz`) |
| `INDEX_FILE` | Block index produced by `mafutils index` |

Options:

| Option | Description |
|---|---|
| `--output-prefix`, `-o` | Output prefix/path |
| `--processes`, `-p` | Number of worker processes |
| `--chunk-size` | Blocks per worker task |
| `--no-block-table` | Skip writing the per-block table |
| `--expected-species` | Comma-separated species names for exact missing lists |
| `--expected-species-file` | File with one species name per line |
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
python -m mafutils gc MAF_FILE [INDEX_FILE] [OPTIONS]
```

Arguments:

| Argument | Description |
|---|---|
| `MAF_FILE` | Input MAF file (`.maf` or `.maf.gz`) |
| `INDEX_FILE` | (Optional) block index produced by `mafutils index` — only used to enable `--processes > 1` |

Options:

| Option | Description |
|---|---|
| `--output-prefix`, `-o` | Output prefix/path |
| `--processes`, `-p` | Number of worker processes (see note below) |
| `--chunk-size` | Blocks per worker task in parallel mode |
| `--log-level` | Logging level: `DEBUG`, `INFO`, `WARNING`, `ERROR` |

**Sequential vs. parallel processing:** GC counting is cheap per base, so the
bottleneck is I/O/decompression rather than CPU. There are two processing
paths:

- **Sequential** (default): streams through the MAF once, top to bottom, with
  no index required. Used whenever no index is given, or the MAF is
  gzip-compressed, or `--processes` is `1`.
- **Parallel**: requires both an `INDEX_FILE` and `--processes > 1`, and is
  only used on an **uncompressed** `.maf`. It chunks index entries across
  worker processes that each seek to their own byte ranges, the same
  architecture as `mafutils stats`.

Parallel mode is skipped for `.maf.gz` input even if requested: the index's
byte offsets are positions in the *decompressed* stream, so each worker would
have to redundantly re-decompress everything up to its own start point,
making things slower rather than faster. A warning is logged and the tool
falls back to the sequential path instead of silently ignoring the request.

## Testing

From inside this `mafutils/` directory:

```bash
pytest tests/test_fetch.py
```

## Notes

- `python -m mafutils ...` is the supported invocation from a source checkout.
- The installed console entrypoint is `mafutils ...` once the package is
  installed.
