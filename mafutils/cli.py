import typer

from . import fetch as fetch_mod
from . import index as index_mod
from . import stats as stats_mod


app = typer.Typer(
    help="MAF toolkit for indexing, fetching, and summarizing MAF alignments.",
    no_args_is_help=True,
    pretty_exceptions_show_locals=False,
    context_settings={"help_option_names": ["-h", "--help"]},
)
app.command("fetch", help="Fetch alignment blocks from a MAF using a BED file and an existing index.")(fetch_mod.fetch_command)
app.command("index", help="Create block and scaffold indexes for a MAF file.")(index_mod.index_command)
app.command("stats", help="Summarize an indexed MAF at overall, species, and block levels.")(stats_mod.stats_command)


def main() -> None:
    app()
