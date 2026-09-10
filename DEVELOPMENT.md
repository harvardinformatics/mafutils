# mafutils Development Notes

Internal implementation notes, gotchas, and maintainer workflows (testing,
releasing) for developing `mafutils` itself — see [`README.md`](README.md)
for user-facing installation/usage docs.

- `tests/example.maf.scaffold.idx` is a **current-format** index, consistent
  with `example.maf.block.idx`. It used to be a deliberately preserved
  headerless/off-by-one fixture, but `fetch` block mode now reads it to order
  regions (see the streaming note below), and a stale pair silently drops
  regions — an off-by-one run boundary is enough to skip every run's first
  block. Tests that need a *headerless* index now build one in `tmp_path` by
  stripping the header (`writeHeaderlessIndex` in `tests/test_validate.py`),
  so no committed file's staleness is load-bearing any more. Keep these two
  index files in sync: regenerate both with `mafutils index` together.
- `tests/example.maf.scaffold.regenerated.idx` is the scaffold index produced
  by the current `mafutils index` implementation.
- The real production scaffold indexes checked so far match the regenerated
  convention, not the preserved older fixture.
- **`mafutils validate` cross-checks the block and scaffold index headers
  against each other**, not just the block index against the MAF file.
  `mafutils index` always writes both from the identical
  size/mtime/hash/compression values in the same run
  (`index.run_index`), so if they differ, the two indexes weren't built
  together (e.g. only one was rebuilt) and shouldn't be trusted as a
  matched pair — this downgrades what would otherwise be a VERIFIED verdict
  to MISMATCH (if headers differ) or UNVERIFIABLE (if the scaffold index is
  missing or headerless), even when the block index matches the MAF file
  perfectly.

## Compression / index internals

- **`mafutils/lib/bgzf.py` is vendored from Biopython, not a dependency.**
  BGZF support needs real random-access `seek()`/`tell()` via virtual
  offsets; the only maintained pure-Python implementation is Biopython's
  `Bio.bgzf`, but depending on `biopython` pulls in `numpy` (a hard
  dependency of biopython itself) for no other reason, adding ~170MB to
  every install. `Bio/bgzf.py` is self-contained (stdlib-only: `io`,
  `struct`, `sys`, `zlib`) and permissively (dual-)licensed, so it's copied
  verbatim into `mafutils/lib/bgzf.py` instead — see the notice at the top
  of that file and `THIRD_PARTY_LICENSES.md` at the repo root. **Tradeoff:**
  this means upstream bug fixes to `Bio/bgzf.py` are not picked up
  automatically (it changes ~1-3 times/year, mostly minor, but not always —
  e.g. a real byte-mode-detection bug fix landed in `BgzfWriter` in May
  2024). If you hit a bgzip-related bug, check
  [upstream's commit history](https://github.com/biopython/biopython/commits/master/Bio/bgzf.py)
  before assuming it's mafutils-specific. Keep this file otherwise
  unmodified so future diffs against upstream stay clean.
- Index files (`.block.idx`/`.scaffold.idx`) carry a `# mafutils-index
  format=2 maf=<name> compression=<mode> size=<bytes> mtime=<epoch>
  hash=<algo>:<hex>` header line (helpers in `mafutils/lib/common.py`:
  `writeIndexHeader`/`readIndexHeader`/`compareIndexHeader`/
  `validateIndexHeader`). Existing parsers already skip `#` lines, so this
  is backward compatible. Format-1 indexes (pre-dating size/mtime/hash)
  simply lack those keys in the parsed header dict — every consumer treats
  their absence as "can't check this dimension," not an error, same as a
  fully missing header. Tests build a headerless index on the fly rather than
  relying on a committed one (see above).
- **Computing the hash costs nothing extra at `mafutils index` time, for
  all three compression types** — this was the whole point of always
  storing it. `common._HashingRawIO` (an `io.RawIOBase` subclass) sits
  *underneath* the existing decompression/decoding layer during indexing
  (`common.openMafHashing`), hashing every byte physically read from disk
  before it's handed onward — for `gz`/`bgzip` this means hashing the raw
  *compressed* on-disk bytes (correct, since the index's offsets are tied
  to that exact representation), not the decompressed content. One
  non-obvious gotcha: `bgzf.BgzfReader` issues `seek()` calls even during
  purely sequential reading (e.g. re-asserting the current position when
  loading each new block) — `_HashingRawIO.seek()` allows only true no-op
  seeks (target == current position) and raises otherwise, since a real
  jump would hash bytes out of order and silently corrupt the digest. This
  was verified empirically against a 27-block/18k-alignment-block synthetic
  bgzip fixture (many internal seek calls), not just small test files.
- `computeFileHash()` is the separate, dedicated-read version (used by
  `mafutils validate` and `--verify-hash`, which only have a file path, not
  a mid-stream position) — there's no read to piggyback on at validation
  time, unlike at build time, so this one really does cost a full read.
- Default (non-`--verify-hash`) validation treats a size mismatch as an
  **error** but an mtime-only mismatch as a **warning** — deliberately
  asymmetric: size differing is a strong, nearly-conclusive signal (content
  is almost certainly different), mtime differing alone is weak and often
  benign (files get touched/copied without content changing). See
  `common.MTIME_TOLERANCE_SECONDS` for the comparison tolerance.
- bgzip is detected distinctly from plain gzip by checking for BGZF's `BC`
  extra-field subfield (magic `\x1f\x8b\x08\x04`, `XLEN=6`, subfield id
  `BC`) after the shared gzip-family magic-byte match — see
  `common._isBgzip`. Verified against real `bgzip`-CLI output, not just
  vendored `bgzf.BgzfWriter`'s own output, since encoder details could in
  principle differ.
- **BGZF virtual offsets cannot be subtracted to get a byte count.** They
  pack a compressed block offset and an in-block offset into one int
  (`bgzf.make_virtual_offset`/`split_virtual_offset`); `offset_end -
  offset_start` is meaningless for them (unlike plain byte offsets or
  gzip's decompressed-stream positions, where it gives a valid read
  length). `common.readMafBlockBytes` handles this by reading forward
  line-by-line until the handle's position (comparable, not subtractable)
  reaches `offset_end` for bgzip, vs. a direct `read(offset_end -
  offset_start)` for the other two. Every block-read call site in
  `fetch.py`/`stats.py`/`gc.py` goes through this helper — don't
  reintroduce a raw `seek()`+subtract read against a bgzip handle.
- **`common.iterMafBlocks` must never compute a block boundary by arithmetic
  on a `tell()` value** (e.g. `tell() - len(line)`) — a real bug of exactly
  this shape shipped in the original bgzip-support work and was only caught
  by testing against a real ~42GB production file. Whenever an `a`-line
  happened to straddle a genuine BGZF block boundary, subtracting a plain
  byte count from the packed virtual offset "borrowed" across the packing
  and produced a byte position that was never a valid block start (observed:
  computed compressed-block-offset 235 vs. true 236, confirmed independently
  via `bgzf.BgzfBlocks()` as ground truth). The fix: capture `stream.tell()`
  immediately before each `readline()` call and use only those captured
  values as boundaries, never derived by subtraction — correct for bgzip's
  packed offsets and equally correct (just simpler) for `none`/`gz`'s plain
  byte/decompressed-stream offsets. `tests/test_compression.py::test_itermafblocks_handles_bgzf_block_boundary_straddling_a_line`
  engineers this exact scenario deterministically (using the vendored
  `bgzf.BgzfWriter`'s 65536-byte block flushing) so it can't regress
  silently. **If you have a bgzip index built before this fix, rebuild it**
  — any block whose header line happened to straddle a BGZF block boundary
  would have a silently corrupted offset.
- **`gc`'s memory usage has a sharp cliff between `--processes 1` and `--processes
  2`, not a smooth scaling curve** — confirmed via real `benchmarks/` data on the
  ~42GB production file: `-p 1` uses ~48-53MB `max_rss`, while `-p 2/4/8` all sit
  at ~1.3-1.5GB regardless of worker count. This is `gc.py`'s own branching
  (`use_parallel = ... and args.processes > 1`, `gc.py:292`): `-p 1` skips
  `ProcessPoolExecutor` entirely and runs a lean single-process loop
  (`runSequentialGC`), while `-p 2+` pays a large, mostly-fixed pool-creation
  cost that doesn't grow further with more workers. `stats` did *not* show this
  cliff (flat ~2.8-3.0GB across every process count including `-p 1`) because it
  used to always construct a `ProcessPoolExecutor` regardless of worker count —
  fixed to match `gc`'s pattern (see below).
- **`stats` now skips `ProcessPoolExecutor` at `--processes 1`, mirroring `gc`'s
  parallel/sequential split above.** `workerTask` is a plain function with no
  pool-initializer dependency (it opens its own MAF handle internally), so
  calling it directly in a loop at `-p 1` is exactly equivalent to submitting it
  to a 1-worker pool, just without the pool's fixed memory overhead. This was a
  real, unnecessary cost for the common case: most invocations don't pass
  `--processes` at all, defaulting to 1.
- **`fetch` on `gz` uses *more* peak memory than `none`/`bgzip` for the identical
  full-BED workload** (~36.2GB vs ~26.3-26.4GB, observed) despite being the
  single-process path, not the 8-worker one. This follows directly from the
  decode-once cache design above (`prefetchBlockCache`): gz holds the *entire*
  needed block cache in one process, while none/bgzip spread the same total
  memory demand across 8 separate worker processes' `max_rss` (which Snakemake's
  `benchmark:` reports as a sum across the process tree). So gz costs more time
  *and* more peak memory than the parallel-eligible compressions for the same
  `fetch` workload — worth knowing before running it on a memory-constrained node.
- **`mafutils index`'s hot loop: `TextIOWrapper.tell()` was the bottleneck, and
  the fix is binary mode + a byte counter, NOT a chunked scanner.** Profiling a
  real 1.5GB MAF slice showed `tell()` at ~60% of runtime — it was called once
  per line (`iterMafBlocks` needs a position before every line) and cost ~4.6x
  more than the `readline()` it accompanied, because on a `TextIOWrapper`
  `tell()` must snapshot the incremental UTF-8 decoder state to build a
  seekable cookie (visible as ~1.1M `codecs.getstate`/`setstate` calls). Since
  `none`/`gz` have plain additive stream positions, `iterMafBlocks(...,
  binary=True)` now reads a binary stream and tracks position with `pos +=
  len(line)`. Measured on 1.5GB of real data, scan-only: **11.89s -> 1.38s
  (8.6x)**; end-to-end `mafutils index`: **23.82s -> 11.63s (2.05x)**, output
  byte-identical.
  Two approaches that sound faster but measured **slower**, so don't retry them
  without new evidence: a chunked `bytes.find(b"\n")` scanner (2.17s — loses to
  C-implemented `BufferedReader.readline()` while adding Python-level buffer
  slicing), and a block-level scanner that finds `b"\na"` boundaries to avoid
  per-line iteration entirely (2.71s — still has to walk every line to count
  non-blank lines for `num_seqs`, plus buffer-concatenation cost). Remaining
  headroom is real but modest: end-to-end sits at ~129 MB/s against a ~640 MB/s
  ceiling (raw read + MD5 with no parsing at all), and closing it would require
  assuming things about intra-block blank/comment lines that risk silently
  corrupting an index.
  **bgzip deliberately stays on the text+`tell()` path** — its virtual offsets
  are not additive, so a byte counter cannot reproduce them (same root cause as
  the straddling-`a`-line bug above).
- **Whole-scaffold extraction must stream; `readMafBlockBytes` is only for
  block-sized ranges.** `fetch -m scaffold` used to do
  `readMafBlockBytes(...)` for a whole scaffold, which materializes the range
  as one bytes object. On a real 7.4TB 241-species MAF, chr1 spans bytes
  16 -> 597,495,640,164, so that became a single **~597GB `read()`** and died
  with a `MemoryError`. Two traps worth knowing: `str(MemoryError())` is the
  **empty string**, so the failure logged as a completely blank
  `[ERROR] chr1:` message; and with memory overcommit the doomed read
  actually starts, so at `-p 1` it looks like a hang rather than an error.
  `common.copyMafRangeToStream()` now copies in bounded chunks (O(chunk_size)
  memory regardless of range size) straight into the output handle. Note its
  bgzip branch also had to change: the old line-by-line loop accumulated
  lines into a list to `join()` at the end, so it was O(range) too --
  fixing only the uncompressed branch would have left bgzip broken at scale.
- **Scaffold extraction failures now fail the run.** Both `fetchByScaffold`
  and `fetchScaffoldsSequential` used to catch any exception, log it, and
  then report `"Wrote <scaffold>"` regardless, so `fetch` exited **0** with a
  0-byte output file -- a downstream Snakemake pipeline consumed that as
  success and only broke several steps later. They now return a status dict,
  the caller `sys.exit(1)`s if any scaffold failed *or* extracted 0 bytes,
  "Wrote" is logged only for genuinely non-empty output, and errors are
  logged with `repr(e)` plus a debug traceback (precisely because the
  exception that actually occurred had an empty `str()`). This was also a
  direct violation of `AGENTS.md`'s "never insert silent fallbacks or error
  handling that could hide errors" rule -- worth re-reading the other
  `except Exception` sites against that standard.
- **`gc` streams the block index rather than loading it.** `gc.parseIndex`
  built a list of every block's `(offset_start, offset_end)`, then
  `list(chunker(...))` **sliced that list**, making a second full copy.
  Measured peak, scaling perfectly linearly: 34MB @ 250k blocks, 136MB @ 1M,
  563MB @ 4M -- extrapolating to **~32GB at the 233M blocks** of a real
  whole-genome index (observed as a 62GB-RSS OOM kill before any GC was
  computed). Replaced by `iterIndexEntries()` (a generator) plus an
  `islice`-based `chunker()` that never materializes its source, with a
  bounded number of in-flight futures and results merged as they arrive via
  `mergeCountsInto()` -- submitting all ~47k chunks upfront, or collecting
  every worker's count dict before merging, would each have reintroduced the
  same scaling. Peak is now **flat at ~1.5MB** regardless of index size.
- **`stats` streams its index the same way**, with three wrinkles `gc` didn't
  have. (1) `block_id` is assigned during index parsing, and the original
  incremented it *before* the `int()` conversions, so a line with enough
  fields that fails to parse still consumes an id; `iterIndexEntries`
  reproduces that exactly, since `block_id` is written into the per-block
  output table. (2) The gz path used to set `chunk_size = len(entries)` to
  force one single task, keeping that one worker's seeks forward-only --
  which requires knowing the entry count upfront, the one thing streaming
  can't do. `workerTask` now takes an optional `maf_fp`, and the
  single-process path opens one handle and reuses it across chunks, which
  preserves forward-only seeking without materializing anything (measured:
  no gz slowdown, 41.3s -> 40.4s). (3) Results are merged on arrival via
  `mergeStatsInto`/`newStatsAccumulator` rather than collected -- ~47k chunk
  results, each carrying a per-species dict, would otherwise scale with
  index size just like the entries list did.
- **`stats`' per-block output rows were being concatenated in lexicographic
  filename order** -- a pre-existing bug, unrelated to streaming, found while
  testing the above. Each chunk writes `maf_stats.blocks.task<N>.tsv` and
  `writeBlockTable` did `sorted(block_tmp_paths)`, so `task10` sorted before
  `task2` and `block.tsv` came out badly out of genomic order for any run
  with >=10 chunks (a real 269k-block run has 54). Results are now tracked
  as `(task_id, path)` and sorted numerically; verified 0 out-of-order
  transitions across 269,308 rows. Note this *changes* `block.tsv` row order
  versus older mafutils for multi-chunk runs -- it was wrong before.
- **`fetch` block mode streams the index too, by sorting regions into index
  order instead of doing random access.** It used to build a 5-key **dict per
  block** (~288 bytes/block, ~67GB at 233M blocks) so it could `bisect` by
  coordinate. Now regions are sorted into index order and each worker walks its
  own span of the index forward, matching regions as it passes them -- regions
  are the cheap side (979k regions vs 233M blocks, ~240x fewer). Measured on
  20k real regions against the 8.09M-block hamster index: **3.30GB -> 47MB peak
  RSS and 8.8x faster** (68.7s -> 7.7s; the old path spent most of its time
  building that structure and pickling it to every worker), with all 20,000
  per-region outputs byte-identical for `none`/`gz`/`bgzip`.
  Four things that make this correct, each of which broke a first attempt:
  1. **Scaffold order is MAF file order, not alphabetical.** Real data runs
     `CM000994.3, GL456210.1, ..., CM000995.3`. The ordering comes from
     `<maf>.scaffold.idx`, whose line order is exactly the block index's.
     Sorting by scaffold *name* would seek backwards through the index.
  2. **A scaffold can occupy several disjoint runs.** MAFs may interleave
     scaffolds -- `tests/example.maf` does (`chr4 ... chrX ... chr4`), so
     `chr4` has two runs. `readScaffoldRuns` returns one entry per *run*, not
     per scaffold; assuming one contiguous run per scaffold silently loses the
     later blocks.
  3. **Runs are bounded by their MAF byte range, not by index byte offsets.**
     `findIndexOffsetForMafByte` binary-searches the index *by byte position*
     (seek to a midpoint, discard the partial line, compare the next complete
     record -- ~35 probes for 35GB, the `look(1)` technique) using
     `offset_start`, the only column that increases monotonically through the
     whole file. Its result is a valid place to *start* scanning, deliberately
     allowed to sit early, so it must never be treated as an exact boundary.
  4. **Offsets passed between helpers are exact record boundaries.**
     `iterIndexRecordsFrom` does NOT skip a leading partial line, because
     resume offsets point at real record starts; skipping one dropped the first
     record of every region after the first.
  Entries reach `fetchByRegion` as an **iterator**, never a list: one BED region
  can span a whole scaffold, so materializing its blocks would reintroduce an
  O(blocks) blowup for a single region. (FASTA output is still inherently
  O(region length x species), since `fasta_seqs` accumulates before writing.)
  Block mode now also requires the block and scaffold indexes to come from the
  same `mafutils index` run, and errors if their headers disagree -- it locates
  blocks via byte ranges recorded in the scaffold index, so a mismatched pair
  silently drops regions.
- **`prefetchBlockCache` is gone.** It existed to give gzip "read every needed
  block once, in ascending file order", but it materialized the decoded text of
  every needed block (a real run logged "Prefetching 1031475 distinct
  block(s)") -- its own O(needed blocks) memory blowup. Streaming the index in
  order gives that property to *all three* compression types structurally, and
  reuse between neighbouring regions is covered by the bounded
  `WORKER_BLOCK_CACHE` LRU, since sorting makes such regions adjacent.
- **`fetch`'s no-output accounting read the wrong column, and it mattered a
  lot.** The per-region summary line is `scaffold, start, end, basename,
  n.overlapping.blocks, ...`; the code checked field **1** (`start`) instead of
  field **4** (`n.overlapping.blocks`), so it was really asking "does this
  region start at coordinate 0?". Two opposite symptoms: a region starting at 0
  that worked fine was counted as "no output" (a one-region BED at coordinate 0
  wrote correct output and then **failed the run**), while a region elsewhere
  that found zero blocks was counted as written. On `tests/example.bed` the two
  errors cancelled to exactly 1/10 = 0.10 -- precisely the threshold, which
  fires on `>` -- so every test passed and the thresholds were never really
  exercised.
- **Non-overlap outcomes are advisory, not fatal.** With the count fixed, the
  leftover hardcoded thresholds from the abandoned `--max-no-overlap-*` feature
  started failing legitimate runs (`example.bed` genuinely has 2/10 = 20%
  non-overlapping). `fetch` is normally pointed at many intervals, so
  `ADVISORY_NO_OVERLAP_REGIONS`/`ADVISORY_NO_OVERLAP_FRACTION`, "all regions
  produced nothing", and "BED names a scaffold absent from the index" all now
  **warn and continue**. Real misuse is still caught loudly and earlier:
  index/MAF size+hash mismatch, a mismatched block/scaffold index pair, and
  per-scaffold extraction failures.
- The vendored `bgzf.BgzfReader` has no `.name` attribute (unlike `gzip.GzipFile`),
  so `fetch.py` gets the MAF's display filename from the known `maf_file`
  path (`WORKER_MAF_FILE`), not from the open file handle.
- Plain gzip can't get real random access no matter how the index is built:
  Python's `gzip.GzipFile.seek()` implements a forward seek as "decompress
  and discard" (roughly free) and a backward seek as a full restart from
  byte 0 (expensive). `stats`/`gc` exploit this by forcing a single
  process and processing all blocks in the index's already-ascending file
  order. `fetch` is trickier since it only needs a sparse subset of
  blocks (whichever BED regions overlap) — it precomputes the global set
  of needed blocks across *all* regions, reads each exactly once in
  ascending order, and serves every region from that cache
  (`prefetchBlockCache`/`WORKER_BLOCK_CACHE`), which sidesteps
  backward-seek risk entirely rather than relying on any per-region
  processing order.

## Development setup

From a source checkout, an editable install picks up code changes without
reinstalling:

```bash
pip install -e .
```

Or run directly against the checkout without installing at all:

```bash
python -m mafutils --help
```

## Testing

From inside this `mafutils/` directory:

```bash
pytest tests/
```

`tests/test_compression.py` covers cross-cutting compression behavior
(bgzip/gzip detection, index headers, index auto-derivation, and the
per-command parallel/sequential/fallback logic) across all three compression
types. `tests/test_validate.py` covers index integrity (size/mtime/hash
header fields, `mafutils validate`'s three-way verdict, and `--verify-hash`
on `fetch`/`stats`/`gc`). `tests/test_stats.py` checks `mafutils stats`'s
computed `overall.tsv`/`species.tsv` values against numbers hand-derived
directly from `tests/example.maf`'s raw content (not just "doesn't crash"),
including a regression guard that the sequential (`-p 1`, no
`ProcessPoolExecutor`) and real multi-worker (`-p 2`+) paths produce
identical output.

`tests/test_real_data.py` runs `stats`/`gc`/`fetch` against
`tests/real-excerpt.maf` -- 8 real, complete alignment blocks extracted
read-only from gwct's actual ~42GB production MAF
(`data/hamsters/uncompressed/...`), not hand-crafted. Unlike
`example.maf`-based tests (which verify exact hand-computed values),
these check plausibility/robustness on genuinely real data -- real
species-naming conventions, real gap patterns, real block-size
distribution -- since a hand-crafted fixture wouldn't think to include
whatever real data actually looks like. This is deliberately still a
tiny, committed fixture, not the real dataset itself: `data/hamsters/`
stays untouched and out of git (see the `.gitignore` `data/` entry).

## Releasing to PyPI

- Publishing is triggered by pushing a Git tag that matches `v*` (for example
  `v0.1.1`).
- Ensure the GitHub repository secret `PYPI_API_TOKEN` is set before releasing.
- Commit release-related changes first, then create and push the tag:

```bash
git add .github/workflows/publish-pypi.yml pyproject.toml README.md DEVELOPMENT.md
git commit -m "release: prepare v0.1.1"
git tag -a v0.1.1 -m "v0.1.1"
git push origin main v0.1.1
```
