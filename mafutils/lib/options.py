#############################################################################
# Command-line options shared by several mafutils subcommands, defined once so
# their name, default, help and validation cannot drift apart between them.
#############################################################################

from typing import Annotated

import typer

from mafutils.lib import common as COMMON


def validateProgressInterval(value: int) -> int:
    """
    0 (off) or at least the default. Anything else is a usage error rather than
    being clamped -- a silently changed value would hide the mistake. The
    minimum applies at the command line only; the run_* functions accept any
    interval, so tests can use short ones.
    """
    if value == 0 or value >= COMMON.PROGRESS_INTERVAL_SECONDS:
        return value
    raise typer.BadParameter(
        f"must be 0 (off) or at least {COMMON.PROGRESS_INTERVAL_SECONDS} seconds, not {value}."
    )


# Per command (index, stats, gc, fetch) rather than one global option: a global
# option must come before the subcommand (`mafutils --progress-interval 300
# stats ...`), so the natural `mafutils stats ... --progress-interval 300` would
# fail, and it would not appear in each command's --help.
ProgressIntervalOption = Annotated[
    int,
    typer.Option(
        "--progress-interval",
        metavar="SECONDS",
        callback=validateProgressInterval,
        help=(
            f"Seconds between progress lines (default and minimum: {COMMON.PROGRESS_INTERVAL_SECONDS}). "
            f"0 turns progress reporting off."
        ),
    ),
]
