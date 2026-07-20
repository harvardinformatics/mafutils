#############################################################################
# Check whether a MAF file's index is still trustworthy.
#
# Gregg Thomas + Codex, July 2026
#############################################################################

"""
mafutils validate

Checks a MAF file against its block index's stored header (size, mtime, and
content hash, from `mafutils index`) and reports whether the index can
still be trusted -- a "check once, trust thereafter" alternative to passing
--verify-hash on every fetch/stats/gc invocation. Also looks up the
scaffold index at its default location (<MAF_FILE>.scaffold.idx) and
cross-checks its header against the block index's -- `mafutils index`
always writes both from the same size/mtime/hash in the same run, so any
difference means the two indexes were not built together (e.g. one was
rebuilt and the other wasn't) and should not be trusted as a matched pair.

Exit codes:
    0  VERIFIED      -- the index's stored hash matches this file exactly,
                         and the block/scaffold index headers agree.
    1  MISMATCH      -- a conclusive difference was found (compression,
                         size, or hash against the MAF file; or the block
                         and scaffold index headers disagree); rebuild the
                         index(es) with `mafutils index`.
    2  UNVERIFIABLE  -- nothing contradicts, but there's no stored hash to
                         be fully sure (the index predates this feature),
                         or the scaffold index is missing/headerless so the
                         pair couldn't be cross-checked.
"""

import logging
import os
import sys
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


def run_validate(args, cmdline="mafutils validate"):
    logger_name = LOGINIT.configureLogging(
        log_level=args.log_level,
        log_verbosity="SCREEN",
        logger_name="maf_validate_logger",
    )
    LOG = logging.getLogger(logger_name)

    LOG.info(f"mafutils validate called as: {cmdline}")

    index_file = args.index_file
    if index_file is None:
        index_file = COMMON.deriveBlockIndexPath(args.maf_file)
        if not os.path.isfile(index_file):
            LOG.error(
                f"No index given and none found at default location ({index_file})."
            )
            sys.exit(1)
        LOG.info(f"Using index at default location: {index_file}")

    if not os.path.isfile(index_file):
        LOG.error(f"Index file not found: {index_file}")
        sys.exit(1)

    detected_compression = COMMON.detectCompression(args.maf_file)
    header = COMMON.readIndexHeader(index_file)

    if header is None:
        LOG.warning(
            f"{index_file} has no mafutils header (built by an older mafutils "
            f"version) -- there is nothing to check it against."
        )
        LOG.warning(
            "VERDICT: UNVERIFIABLE -- rebuild the index with `mafutils index` "
            "to enable verification."
        )
        sys.exit(2)

    result = COMMON.compareIndexHeader(header, args.maf_file, detected_compression, check_hash=True)

    conclusive_mismatch = False
    for field in ("compression", "size", "hash"):
        check = result[field]
        if check is None:
            LOG.info(f"{field:11}: not recorded in index header")
            continue
        if check["match"]:
            LOG.info(f"{field:11}: match (index={check['expected']})")
        else:
            LOG.error(f"{field:11}: MISMATCH (index={check['expected']}, actual={check['actual']})")
            conclusive_mismatch = True

    mtime_check = result["mtime"]
    if mtime_check is None:
        LOG.info(f"{'mtime':11}: not recorded in index header")
    elif mtime_check["match"]:
        LOG.info(f"{'mtime':11}: match (index={mtime_check['expected']})")
    else:
        LOG.warning(
            f"{'mtime':11}: differs (index={mtime_check['expected']}, "
            f"actual={mtime_check['actual']}) -- not necessarily a problem; "
            f"the hash check above is the authoritative signal."
        )

    scaffold_index_file = COMMON.deriveScaffoldIndexPath(args.maf_file)
    scaffold_unverifiable = False

    if not os.path.isfile(scaffold_index_file):
        LOG.warning(
            f"Scaffold index not found at default location ({scaffold_index_file}); "
            f"cannot cross-check it against the block index."
        )
        scaffold_unverifiable = True
    else:
        scaffold_header = COMMON.readIndexHeader(scaffold_index_file)
        if scaffold_header is None:
            LOG.warning(
                f"{scaffold_index_file} has no mafutils header (built by an older "
                f"mafutils version); cannot cross-check it against the block index."
            )
            scaffold_unverifiable = True
        elif scaffold_header == header:
            LOG.info(f"{'pair':11}: block and scaffold index headers match")
        else:
            LOG.error(
                f"{'pair':11}: MISMATCH -- block index ({index_file}) and scaffold "
                f"index ({scaffold_index_file}) headers differ, so they were not "
                f"built together (block={header}, scaffold={scaffold_header}). "
                f"Rebuild both with `mafutils index`."
            )
            conclusive_mismatch = True

    if conclusive_mismatch:
        LOG.error(
            f"VERDICT: MISMATCH -- index for {args.maf_file} does not match this "
            f"file. Rebuild it with `mafutils index`."
        )
        sys.exit(1)

    if result["hash"] is not None and not scaffold_unverifiable:
        LOG.info(
            f"VERDICT: VERIFIED -- {args.maf_file} matches its index exactly, "
            f"and the block/scaffold index headers agree."
        )
        sys.exit(0)
    else:
        LOG.warning(
            "VERDICT: UNVERIFIABLE -- nothing contradicts, but this index has no "
            "stored hash to be fully sure (built by an older mafutils version), "
            "or the scaffold index couldn't be cross-checked. Rebuild the "
            "index(es) with `mafutils index` to enable full verification."
        )
        sys.exit(2)


def validate_command(
    maf_file: Annotated[str, typer.Argument(help="Input MAF file (.maf, .maf.gz, or bgzip-compressed .maf)")],
    index_file: Annotated[
        Optional[str],
        typer.Argument(help="Block index from mafutils index. Required -- if omitted, looked up at <MAF_FILE>.block.idx by default; errors if not found there. The scaffold index is always looked up at <MAF_FILE>.scaffold.idx and its header is cross-checked against this one."),
    ] = None,
    log_level: Annotated[LogLevel, typer.Option("--log-level", help="Logging level (default: INFO)")] = LogLevel.INFO,
) -> None:
    args = SimpleNamespace(
        maf_file=maf_file,
        index_file=index_file,
        log_level=log_level.value,
    )
    cmdline = "mafutils validate"
    if len(sys.argv) > 2:
        cmdline += " " + " ".join(sys.argv[2:])
    run_validate(args, cmdline=cmdline)
