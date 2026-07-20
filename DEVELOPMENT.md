# mafutils Development Notes

- `tests/example.maf.scaffold.idx` is preserved as an older, headerless
  scaffold-index fixture for comparison — it deliberately sits at
  `example.maf`'s *default* scaffold-index path
  (`common.deriveScaffoldIndexPath`), so `mafutils validate` against the
  checked-in `example.maf`/`example.maf.block.idx` pair reports
  UNVERIFIABLE rather than VERIFIED (it can't cross-check a headerless
  scaffold index against the block index — see below). Tests exercising
  the "clean, fully verified match" case build a fresh pair with
  `mafutils index` in a tmp dir instead of using this fixture pair directly.
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
  fully missing header. `tests/example.maf.scaffold.idx` is deliberately
  kept in format-1 shape (see above) specifically to exercise this path.
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
