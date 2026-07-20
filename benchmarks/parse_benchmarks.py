"""
Shared parsing for results/{compression}/{command}[.headline|.sweep.p{N}].benchmark.tsv
files -- used by both the Snakefile's record_history rule and notebook.ipynb,
so the two can't drift out of sync with each other.
"""

import re
from pathlib import Path

import pandas as pd

FILENAME_RE = re.compile(
    r"^(?P<command>[a-z]+)(?:\.(?P<kind>headline|sweep)(?:\.p(?P<processes>\d+))?)?\.benchmark\.tsv$"
)


def load_benchmarks(results_dir, headline_processes):
    """Load every results/{compression}/*.benchmark.tsv under results_dir into one long DataFrame."""
    rows = []
    for tsv_path in sorted(Path(results_dir).glob("*/*.benchmark.tsv")):
        compression = tsv_path.parent.name
        m = FILENAME_RE.match(tsv_path.name)
        if not m:
            continue
        fields = m.groupdict()
        df = pd.read_csv(tsv_path, sep="\t")
        df["compression"] = compression
        df["command"] = fields["command"]
        df["kind"] = fields["kind"] or "single"
        df["processes"] = int(fields["processes"]) if fields["processes"] else headline_processes
        rows.append(df)

    if not rows:
        return pd.DataFrame(
            columns=["s", "max_rss", "mean_load", "compression", "command", "kind", "processes"]
        )

    bench = pd.concat(rows, ignore_index=True)
    if "mean_load" not in bench.columns or bench["mean_load"].isna().all():
        bench["mean_load"] = 100 * bench["cpu_time"] / bench["s"]
    return bench
