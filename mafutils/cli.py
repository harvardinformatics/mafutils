import typer

from . import fetch as fetch_mod
from . import gc as gc_mod
from . import index as index_mod
from . import stats as stats_mod
from . import validate as validate_mod


app = typer.Typer(
    help="MAF toolkit for indexing, fetching, and summarizing MAF alignments.",
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)
app.command("fetch", help="Fetch alignment blocks from a MAF using a BED file and an existing index.")(fetch_mod.fetch_command)
app.command("gc", help="Calculate per-species GC content from a MAF file.")(gc_mod.gc_command)
app.command("index", help="Create block and scaffold indexes for a MAF file.")(index_mod.index_command)
app.command("stats", help="Summarize an indexed MAF at overall, species, and block levels.")(stats_mod.stats_command)
app.command("validate", help="Check whether a MAF file's index is still trustworthy.")(validate_mod.validate_command)


def main() -> None:
    app()
