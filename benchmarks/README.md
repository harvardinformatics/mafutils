# mafutils benchmarks

Snakemake pipeline that times `index`/`validate`/`gc`/`stats`/`fetch` across
compression types and process counts against gwct's real ~42GB production
data (`data/hamsters/`), plus a notebook that turns the results into plots
and tables.

`baseline-manual-timings.txt` is the original hand-run (`time -p mafutils
...`) log this pipeline formalizes -- kept as a sanity-check baseline.

## Setup

None of the existing conda envs in this project have the full set needed.
Install via conda/mamba, not pip: Snakemake's own docs recommend the
conda/bioconda install (it pulls in several compiled/scheduling
dependencies -- e.g. `pulp` -- that conda resolves more reliably than pip),
and mixing pip and conda installs in one env risks the conda solver not
knowing about pip-installed constraints on a later `conda update`/install.
Since `mafutils-dev` already has `mafutils` itself, add everything else
there in one shot with the same tool:

```bash
conda activate mafutils-dev
mamba install -c conda-forge -c bioconda snakemake snakemake-executor-plugin-slurm pandas matplotlib pyyaml jupyter nbconvert ipykernel
# or, if mamba isn't installed: conda install ...
```

(`snakemake-executor-plugin-slurm` is only needed for cluster execution,
see below -- harmless to install either way.)

## Running

Locally:

```bash
cd benchmarks
snakemake --cores 8 --resources io_none=1 io_gz=1 io_bgzip=1 io_heavy=2
```

The `--resources` flags enforce the per-compression file-read locks and the
global heavy-I/O cap (see below) -- without them, Snakemake has no global
cap to serialize against and would happily run several jobs against the
same MAF, or several heavy jobs against different MAFs, concurrently.

Via SLURM (recommended -- see below for why):

```bash
cd benchmarks
snakemake --workflow-profile profiles/slurm
```

`profiles/slurm/config.yaml` has placeholder `slurm_account`/`slurm_partition`
values that need filling in for your cluster before this will submit
anything. Per-rule `mem_mb`/`runtime` in the Snakefile have been calibrated
from a real completed run's `results/*.benchmark.tsv` values (comments on
each rule cite the observed peak), not just guessed from
`baseline-manual-timings.txt` -- re-tighten them further as more runs
complete, especially after any code change that could plausibly shift a
command's memory profile. `threads` is set to match each rule's
`--processes` value (or 1 for single-process commands) so SLURM reserves
the right number of cores for `mafutils`'s internal multiprocessing pool.

If you've just changed `mafutils` itself (not this Snakefile), note that
Snakemake has no visibility into that -- it only tracks its own declared
rule inputs/outputs/shell text, not the installed package's source. A
rerun will skip every rule whose declared inputs/outputs haven't changed,
even if the underlying `mafutils` command would now behave differently.
Force a re-run of the affected rules explicitly, e.g.:
```bash
snakemake --workflow-profile profiles/slurm --forcerun stats_headline stats_sweep
```

This never touches the real MAFs or their indexes under `data/hamsters/` --
it only reads them. All pipeline output (benchmark TSVs, rebuilt scratch
indexes, command outputs, the rendered notebook) goes under
`benchmarks/results/` and `benchmarks/summary.html`, both gitignored;
`benchmarks/history.csv` is the one exception -- see below.

Given some cells take upwards of an hour on the full dataset (see
`baseline-manual-timings.txt`), running this via SLURM rather than
locally/interactively is recommended.

## Design decisions

- **`index` writes to a throwaway path** (`results/{compression}/rebuilt.*.idx`,
  `temp()`-marked) rather than the real `<MAF>.block.idx`/`.scaffold.idx`
  next to each MAF in `data/hamsters/` -- those are gwct's own
  already-verified production indexes and must never be rebuilt/overwritten
  by this pipeline. `validate`/`gc`/`stats`/`fetch` all reuse those existing
  indexes directly (read-only) rather than depending on the rebuilt one, so
  the (multi-minute, for `gz`) index cost isn't paid again for every
  downstream rule.
- **`gz` has no `--processes` sweep.** `gc`/`stats`/`fetch` on gzip-compressed
  input always fall back to a single sequential pass regardless of
  `--processes` (see `../DEVELOPMENT.md`), so sweeping process counts would
  just re-measure the same sequential run several times over at real cost
  (hours, on this dataset) for no new information. `gz` still gets one
  headline data point (`--processes 8`, matching
  `baseline-manual-timings.txt`), plotted as a flat reference line in the
  scaling charts.
- **`fetch` writes real one-file-per-region output** (`mafutils fetch`'s
  default), not `--single-output` -- `--single-output` is a materially
  cheaper operation (one open handle, sequential writes) than the default
  per-region open/close + directory-metadata path, so it would benchmark a
  different, lighter operation than what `fetch` normally does. Output
  directories are `temp()`-wrapped so Snakemake cleans them up once no
  longer needed; gwct's own manual runs already produced this same output
  once under `data/hamsters/*/combined-ucsc-peaks*-fa/`, so this is a
  second, throwaway copy purely for timing. Consequently the headline
  `fetch` runs are genuinely expensive: real multi-file output across the
  full BED means ~2.9M temp file creates and deletes total across the 3
  compressions, and roughly the same wall-clock as
  `baseline-manual-timings.txt`'s manual runs (up to ~70 minutes for the
  slowest single compression).
- **`fetch`'s process-count sweep uses a small BED subset**
  (`fetch_sweep_regions` in `config.yaml`, default 5000 regions) rather than
  the full ~979k-region BED, so the sweep matrix stays fast; the headline
  `fetch` run still uses the full BED to stay comparable with
  `baseline-manual-timings.txt`. The subset is a reproducible (fixed-seed)
  uniform random sample of the full BED, not its first N lines -- the BED
  is scaffold-ordered, so `head -n N` alone would under-represent the real
  diversity of scaffolds/region sizes and could skew sweep timings relative
  to the full-BED headline run.
- **`fetch_sweep_gz` runs `gz` on that same 5000-region sweep BED**, giving
  it a real, directly comparable reference point for the process-count
  scaling plot -- reusing `gz`'s full-BED *headline* number there (as
  `gc`/`stats` do) isn't apples-to-apples, since none/bgzip's *sweep* runs
  are on a ~196x smaller BED, and tried normalizing by region count doesn't
  fix that either (the smaller sweep's per-region time is dominated by fixed
  per-run overhead in a way the full run isn't). `--processes` is ignored
  for `gz` regardless, so this is one job, not a real sweep, matching
  `gc_headline`/`stats_headline`'s existing precedent for `gz`'s single data
  point.
- **`stats`'s sweep cells write the full block table**, same as the headline
  run, varying only `--processes` -- writing it is a real (largely
  non-parallelized) part of what `stats` costs, so skipping it would
  benchmark a different, lighter operation than what `--processes` actually
  sweeps over, and would hide the Amdahl's-law-style flattening that fixed
  serial write cost should produce as process count increases.
- Heavy non-benchmark outputs (stats tables, gc CSVs, fetch output) are all
  `temp()`-marked so Snakemake cleans them up once no longer needed --
  only the small benchmark TSVs, `history.csv`, and the final `summary.html`
  persist.
- **`history.csv` is appended to, not regenerated,** by the `record_history`
  rule -- a deliberate deviation from normal Snakemake idempotency norms.
  Every full pipeline run gets tagged with the current git commit (plus
  `-dirty` if the working tree has uncommitted changes) and a UTC timestamp,
  and appended as new rows, so performance can be tracked across code
  changes over time rather than only compared within a single run. Unlike
  `results/`, it's intentionally **not** gitignored -- it's meant to be
  committed and grow. `notebook.ipynb`'s final section reads it and plots
  wall-clock time across all recorded runs (skipped gracefully if fewer
  than 2 runs are recorded yet).
- **Every rule that reads a MAF declares a per-compression file lock**
  (`IO_LOCKS` in the Snakefile: `io_none`/`io_gz`/`io_bgzip`, each capped at
  1 globally -- see `profiles/slurm/config.yaml`'s `resources:` key, or
  `--resources io_none=1 io_gz=1 io_bgzip=1` for local runs). At most one
  job reads a given MAF at a time, so concurrent readers of the *same* file
  can't inflate each other's measured time (an early SLURM run showed
  ~7-18% slower results than the solo `baseline-manual-timings.txt` numbers,
  plausibly from exactly this). Jobs on *different* compressions are
  unaffected and still run fully concurrently with each other -- this
  serializes per-file access, not the whole pipeline.
- **A global `io_heavy` cap (`IO_HEAVY` in the Snakefile, default 2 --
  `config.yaml`'s `io_heavy_concurrency`) additionally limits how many
  heavy MAF-reading jobs run at once ACROSS ALL compressions combined.**
  `io_locks` alone only stops two jobs reading the *same* file at once; a
  real run showed an `index` build and an unrelated `gc` sequential pass
  (different compressions, different files) still both running slower
  while overlapping, from contending for the same underlying shared Lustre
  bandwidth. Lower `io_heavy` to 1 for maximally clean/isolated timings, at
  the cost of much longer total wall-clock.
- **No `--keep-going`.** A failure anywhere stops all further job
  submission, on purpose -- a partial, silently-incomplete-looking result
  set is worse than an obvious hard stop, even though it means independent,
  unrelated branches of the DAG sit idle until the failure is dealt with.
- **`stats` skips `ProcessPoolExecutor` entirely at `--processes 1`**, matching
  `gc`'s existing parallel/sequential split (`mafutils/stats.py`, alongside
  `mafutils/gc.py`'s `use_parallel` pattern) -- pool creation carries a large,
  mostly-fixed memory cost regardless of true worker count (`gc`'s own data
  showed ~48MB with no pool vs ~1.3-1.5GB with one, at any worker count),
  which `stats` used to pay even when asked to run on a single process. This
  is a `mafutils` code fix discovered *from* this benchmark data, not a
  benchmark-pipeline change -- see `../DEVELOPMENT.md`.
