#############################################################################
# Calculate per-species GC content from an indexed MAF file.
#
# Gregg Thomas + Codex, July 2026
#############################################################################

"""
mafutils gc

Computes overall GC content per species across every sequence in a MAF file,
and the unweighted mean GC content across species.

Outputs:
    <prefix>.gc.csv      : species,gc (fraction 0-1)
    <prefix>.gc.mean.txt : single float, the mean of the gc column above

GC is computed as (G+C) / (A+C+G+T), case-insensitive; gaps, N, and other
ambiguity codes are excluded from both numerator and denominator.

Two processing paths:
    Sequential (default, required for gzip input):
        - Streams through the MAF once, top to bottom, no seeking.
        - Does not require a block index.

    Parallel (opt-in speedup, uncompressed input only):
        - Requires a block index (from `mafutils index`) and --processes > 1.
        - Chunks index entries across worker processes that seek to their
          own byte ranges independently, same architecture as `mafutils stats`.
        - Not used for gzip input: the index's byte offsets are positions in
          the decompressed stream, so parallel seeking would force each
          worker to redundantly re-decompress everything before its start
          point. A warning is logged and the sequential path is used instead.
"""

import gzip
import logging
import os
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from enum import Enum
from types import SimpleNamespace
from typing import Annotated, Optional

import typer

from mafutils.lib import common as COMMON
from mafutils.lib import loginit as LOGINIT

#############################################################################


class LogLevel(str, Enum):
    DEBUG = "DEBUG"
    INFO = "INFO"
    WARNING = "WARNING"
    ERROR = "ERROR"


#############################################################################


def speciesFromSrc(src):
    if "." in src:
        return src.split(".", 1)[0]
    return src


def chromFromSrc(src):
    if "." in src:
        return src.split(".", 1)[1]
    return src


def newBaseCounts():
    # Named (not lambda) so defaultdicts using it can be pickled across process boundaries.
    return {"A": 0, "C": 0, "G": 0, "T": 0}


def updateCountsFromSeqLine(line, counts):
    """
    Parses one MAF 's' line and adds its A/C/G/T tallies (case-insensitive;
    gaps, N, and other ambiguity codes ignored) into counts[(species, chrom)].
    """
    fields = line.split()
    if not fields or fields[0] != "s" or len(fields) < 7:
        return

    src = fields[1]
    species = speciesFromSrc(src)
    chrom = chromFromSrc(src)
    seq_counts = Counter(fields[6].upper())

    rec = counts[(species, chrom)]
    rec["A"] += seq_counts["A"]
    rec["C"] += seq_counts["C"]
    rec["G"] += seq_counts["G"]
    rec["T"] += seq_counts["T"]


#############################################################################


def parseIndex(index_file, LOG):
    """
    Parses a block-level index file (same 8-column format as mafutils stats).
    Only offset_start/offset_end are used here.
    """
    entries = []
    with open(index_file, "r", encoding="utf-8") as fp:
        for line in fp:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 8:
                LOG.warning(f"Skipping malformed index line: {line.strip()}")
                continue
            try:
                offset_start = int(fields[6])
                offset_end = int(fields[7])
            except ValueError:
                LOG.warning(f"Skipping malformed index line: {line.strip()}")
                continue
            entries.append((offset_start, offset_end))
    return entries


def chunker(seq, size):
    for i in range(0, len(seq), size):
        yield seq[i : i + size]


#############################################################################


def runSequentialGC(maf_file, maf_compression, LOG):
    counts = defaultdict(newBaseCounts)

    opener = gzip.open if maf_compression == "gz" else open
    with opener(maf_file, "rt") as maf_fp:
        for line in maf_fp:
            if line.startswith("s"):
                updateCountsFromSeqLine(line, counts)

    return counts


#############################################################################


def workerTask(task_id, maf_file, maf_compression, entries):
    counts = defaultdict(newBaseCounts)

    opener = gzip.open if maf_compression == "gz" else open
    with opener(maf_file, "rb") as maf_fp:
        for offset_start, offset_end in entries:
            maf_fp.seek(offset_start)
            block_bytes = maf_fp.read(offset_end - offset_start)
            block_text = block_bytes.decode("utf-8", errors="replace")
            for line in block_text.splitlines():
                if line.startswith("s"):
                    updateCountsFromSeqLine(line, counts)

    return counts


def mergeCounts(results):
    merged = defaultdict(newBaseCounts)
    for counts in results:
        for key, rec in counts.items():
            m = merged[key]
            m["A"] += rec["A"]
            m["C"] += rec["C"]
            m["G"] += rec["G"]
            m["T"] += rec["T"]
    return merged


def runParallelGC(maf_file, maf_compression, index_file, processes, chunk_size, LOG):
    LOG.info(f"Parsing index file: {index_file}")
    entries = parseIndex(index_file, LOG)
    if not entries:
        LOG.error("No valid index entries found.")
        sys.exit(1)
    LOG.info(f"Loaded {len(entries)} indexed blocks.")

    chunks = list(chunker(entries, chunk_size))
    LOG.info(f"Running {len(chunks)} tasks with chunk size {chunk_size} across {processes} process(es).")

    results = []
    with ProcessPoolExecutor(max_workers=processes) as executor:
        futures = [
            executor.submit(workerTask, task_id, maf_file, maf_compression, chunk)
            for task_id, chunk in enumerate(chunks, start=1)
        ]
        for future in futures:
            results.append(future.result())

    return mergeCounts(results)


#############################################################################


def computeSpeciesGC(counts, LOG):
    """
    Rolls up (species, chrom) counts to per-species totals and computes
    gc = (G+C) / (A+C+G+T) for each species.
    """
    species_totals = defaultdict(newBaseCounts)
    for (species, _chrom), rec in counts.items():
        s = species_totals[species]
        s["A"] += rec["A"]
        s["C"] += rec["C"]
        s["G"] += rec["G"]
        s["T"] += rec["T"]

    species_gc = {}
    for species, rec in species_totals.items():
        total = rec["A"] + rec["C"] + rec["G"] + rec["T"]
        if total == 0:
            LOG.warning(f"Species '{species}' has zero countable (A/C/G/T) bases; reporting gc=0.0.")
            species_gc[species] = 0.0
        else:
            species_gc[species] = (rec["G"] + rec["C"]) / total

    return species_gc


def writeSpeciesCSV(out_path, species_gc):
    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write("species,gc\n")
        for species in sorted(species_gc):
            fp.write(f"{species},{species_gc[species]:.8f}\n")


def writeMeanFile(out_path, species_gc):
    values = list(species_gc.values())
    mean_gc = sum(values) / len(values) if values else 0.0
    with open(out_path, "w", encoding="utf-8") as fp:
        fp.write(f"{mean_gc:.8f}\n")
    return mean_gc


#############################################################################


def run_gc(args, cmdline="mafutils gc"):
    out_prefix = args.output_prefix
    out_dir = os.path.dirname(out_prefix) if os.path.dirname(out_prefix) else "."
    os.makedirs(out_dir, exist_ok=True)

    log_file = out_prefix + ".log"
    logger_name = LOGINIT.configureLogging(
        log_level=args.log_level,
        log_verbosity="BOTH",
        log_filename=log_file,
        logger_name="maf_gc_logger",
        overwrite_log_file=True,
    )
    LOG = logging.getLogger(logger_name)

    LOG.info(f"mafutils gc called as: {cmdline}")
    LOG.info("-" * 40)
    for arg, value in vars(args).items():
        LOG.info(f"{arg:16} : {value}")
    LOG.info("-" * 40)

    if args.processes < 1:
        LOG.error("--processes must be >= 1")
        sys.exit(1)
    if args.chunk_size < 1:
        LOG.error("--chunk-size must be >= 1")
        sys.exit(1)

    maf_compression = COMMON.detectCompression(args.maf_file)
    if maf_compression not in ("none", "gz"):
        LOG.error(f"Unsupported MAF compression: {maf_compression}")
        sys.exit(1)

    use_parallel = args.index_file is not None and maf_compression == "none" and args.processes > 1

    if args.index_file is not None and maf_compression != "none" and args.processes > 1:
        LOG.warning(
            "Parallel processing requires random access into the MAF, which is inefficient on "
            "gzip-compressed files: the index's byte offsets are positions in the decompressed "
            "stream, so each worker would have to redundantly re-decompress everything before "
            "its own start point. Falling back to a single sequential pass; --processes will be "
            "ignored. (bgzip support that would allow safe parallel access is planned separately.)"
        )
    elif args.index_file is None and maf_compression == "none" and args.processes > 1:
        LOG.warning(
            "An index (from `mafutils index`) is required for parallel processing. "
            "Falling back to a single sequential pass; --processes will be ignored."
        )

    if use_parallel:
        LOG.info("Running in PARALLEL mode (uncompressed MAF + index + processes > 1)")
        counts = runParallelGC(
            args.maf_file, maf_compression, args.index_file, args.processes, args.chunk_size, LOG
        )
    else:
        LOG.info("Running in SEQUENTIAL mode")
        counts = runSequentialGC(args.maf_file, maf_compression, LOG)

    species_gc = computeSpeciesGC(counts, LOG)
    if not species_gc:
        LOG.error("No species found in MAF file.")
        sys.exit(1)

    csv_file = out_prefix + ".gc.csv"
    mean_file = out_prefix + ".gc.mean.txt"

    LOG.info(f"Writing per-species GC: {csv_file}")
    writeSpeciesCSV(csv_file, species_gc)

    LOG.info(f"Writing mean GC: {mean_file}")
    mean_gc = writeMeanFile(mean_file, species_gc)

    LOG.info(f"Computed GC for {len(species_gc)} species; mean GC = {mean_gc:.8f}")
    LOG.info("Done.")


def gc_command(
    maf_file: Annotated[str, typer.Argument(help="Input MAF file (.maf or .maf.gz)")],
    index_file: Annotated[
        Optional[str],
        typer.Argument(
            help="Block index from mafutils index (optional; only used to enable --processes > 1 on uncompressed MAFs)"
        ),
    ] = None,
    output_prefix: Annotated[str, typer.Option("--output-prefix", "-o", help="Output prefix/path (default: maf_gc)")] = "maf_gc",
    processes: Annotated[int, typer.Option("--processes", "-p", help="Number of worker processes; only effective with an index on an uncompressed MAF (default: 1)")] = 1,
    chunk_size: Annotated[int, typer.Option("--chunk-size", help="Blocks per worker task in parallel mode (default: 5000)")] = 5000,
    log_level: Annotated[LogLevel, typer.Option("--log-level", help="Logging level (default: INFO)")] = LogLevel.INFO,
) -> None:
    args = SimpleNamespace(
        maf_file=maf_file,
        index_file=index_file,
        output_prefix=output_prefix,
        processes=processes,
        chunk_size=chunk_size,
        log_level=log_level.value,
    )
    cmdline = "mafutils gc"
    if len(sys.argv) > 2:
        cmdline += " " + " ".join(sys.argv[2:])
    run_gc(args, cmdline=cmdline)
