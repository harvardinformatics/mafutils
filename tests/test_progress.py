import hashlib
import io
import logging
import os
import subprocess
import sys

import pytest

from mafutils.lib import common as COMMON

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))
MAF_FILE = os.path.join(TEST_DIR, "example.maf")


class FakeClock:
    def __init__(self):
        self.now = 1000.0

    def __call__(self):
        return self.now


class ListHandler(logging.Handler):
    def __init__(self):
        super().__init__()
        self.messages = []

    def emit(self, record):
        self.messages.append(record.getMessage())


@pytest.fixture
def log():
    logger = logging.getLogger("test_progress")
    logger.setLevel(logging.INFO)
    handler = ListHandler()
    logger.addHandler(handler)
    yield logger, handler.messages
    logger.removeHandler(handler)


# ---------------------------------------------------------------------------
# BlockProgress
# ---------------------------------------------------------------------------


def test_progress_is_silent_until_the_interval_passes(log):
    logger, messages = log
    clock = FakeClock()
    progress = COMMON.BlockProgress(logger, 1000, interval=60, clock=clock)
    for _ in range(5):
        clock.now += 10
        progress.update(100)
    assert messages == []
    clock.now += 10
    progress.update(100)
    assert messages == ["Progress: 600 / 1,000 blocks (60.0%), 10 blocks/s, 1m 00s elapsed, ~40s left"]


def test_progress_logs_at_most_once_per_interval(log):
    logger, messages = log
    clock = FakeClock()
    progress = COMMON.BlockProgress(logger, 10_000, interval=60, clock=clock)
    for _ in range(100):
        clock.now += 3          # 300s of updates
        progress.update(10)
    assert len(messages) == 5
    assert progress.done == 1000    # every block counted, logged or not


def test_progress_without_a_total_reports_no_percentage_or_estimate(log):
    logger, messages = log
    clock = FakeClock()
    progress = COMMON.BlockProgress(logger, None, interval=60, clock=clock)
    clock.now += 60
    progress.update(1200)
    assert messages == ["Progress: 1,200 blocks, 20 blocks/s, 1m 00s elapsed"]


def test_progress_omits_estimate_once_total_is_reached(log):
    logger, messages = log
    clock = FakeClock()
    progress = COMMON.BlockProgress(logger, 500, interval=60, clock=clock)
    clock.now += 60
    progress.update(500)
    assert messages == ["Progress: 500 / 500 blocks (100.0%), 8 blocks/s, 1m 00s elapsed"]


def test_byte_progress_reports_sizes_and_throughput(log):
    logger, messages = log
    clock = FakeClock()
    progress = COMMON.ByteProgress(logger, 42_700_000_000, interval=60, clock=clock)
    clock.now += 100
    progress.updateTo(12_300_000_000)
    assert messages == [
        "Progress: 12.30 GB / 42.70 GB (28.8%), 123.00 MB/s, 1m 40s elapsed, ~4m 07s left"
    ]


def test_interval_defaults_to_the_module_constant_at_construction(log, monkeypatch):
    logger, _messages = log
    monkeypatch.setattr(COMMON, "PROGRESS_INTERVAL_SECONDS", 5.0)
    assert COMMON.BlockProgress(logger, None).interval == 5.0


@pytest.mark.parametrize("header,expected", [
    ({"blocks": "8094203"}, 8094203),
    ({"size": "777"}, None),          # index built before the aggregates
    (None, None),                     # headerless index
    ({"blocks": "lots"}, None),       # corrupt value: warned about
])
def test_header_block_count(header, expected, log):
    logger, messages = log
    assert COMMON.headerBlockCount(header, logger) == expected
    if expected is None:
        assert len(messages) == 1     # always says why there is no percentage


@pytest.mark.parametrize("fixture,compression", [
    ("example.maf", "none"), ("example.maf.gz", "gz"), ("example.maf.bgz", "bgzip"),
])
def test_disk_position_is_compressed_bytes_read(fixture, compression):
    """
    Progress for index and sequential gc is the position in the file ON DISK,
    so for gz/bgzip it must count compressed bytes (comparable with the file
    size), from both openers -- openMafHashing's wrapper has to pass fileno()
    through for this to work at all.
    """
    path = os.path.join(TEST_DIR, fixture)
    size = os.path.getsize(path)
    with COMMON.openMaf(path, compression, "rt") as fp:
        # bgzip loads its first block on open, so this need not be 0 -- it
        # runs ahead of parsing by one read buffer, never past the file.
        assert 0 <= COMMON.diskPosition(fp) <= size
        while fp.read(65536):       # bgzf refuses an unbounded read()
            pass
        assert COMMON.diskPosition(fp) == size
    with COMMON.openMafHashing(path, compression, hashlib.md5()) as fp:
        while fp.read(65536):
            pass
        assert COMMON.diskPosition(fp) == size


@pytest.mark.parametrize("seconds,expected", [
    (0, "0s"), (59.4, "59s"), (60, "1m 00s"), (125, "2m 05s"),
    (3600, "1h 00m"), (3725, "1h 02m"), (50 * 3600 + 180, "50h 03m"),
])
def test_format_duration(seconds, expected):
    assert COMMON.formatDuration(seconds) == expected


# ---------------------------------------------------------------------------
# Parameter block in the run log: colons line up whatever the longest name
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("command,longest", [
    ("stats", "dashboard_max_block_points"),
    ("gc", "output_prefix"),
])
def test_parameter_log_colons_are_aligned(command, longest, tmp_path):
    prefix = os.path.join(str(tmp_path), "out")
    r = subprocess.run([sys.executable, "-m", "mafutils", command, MAF_FILE, "-o", prefix],
                       check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    log_text = open(prefix + ".log").read()
    block = log_text.split("-" * 40)[1]     # between the two separator lines
    lines = [ln for ln in block.splitlines() if " : " in ln]
    assert any(longest in ln for ln in lines)
    assert len({ln.index(" : ") for ln in lines}) == 1, "\n".join(lines)


# ---------------------------------------------------------------------------
# index now reports what it did
# ---------------------------------------------------------------------------


def test_index_reports_size_and_result(tmp_path):
    out = os.path.join(str(tmp_path), "x")
    r = subprocess.run([sys.executable, "-m", "mafutils", "index", MAF_FILE, out + ".block.idx", out + ".scaffold.idx"],
                       check=False, cwd=REPO_ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert "(777 B, none). Progress is logged every 60s." in r.stdout
    assert "Indexed 10 blocks on 5 reference scaffold(s)" in r.stdout
    assert not os.path.exists(os.path.join(REPO_ROOT, "maf_index.log"))   # screen only, no log file


# ---------------------------------------------------------------------------
# fetch: progress units, the wait() timeout, and the 8 MB copy hook
# ---------------------------------------------------------------------------


def test_progress_unit_names_what_is_counted(log):
    logger, messages = log
    clock = FakeClock()
    progress = COMMON.BlockProgress(logger, 1000, interval=60, clock=clock, unit="regions")
    clock.now += 60
    progress.updateTo(600)
    assert messages == ["Progress: 600 / 1,000 regions (60.0%), 10 regions/s, 1m 00s elapsed, ~40s left"]


def test_due_in_counts_down_to_the_next_line(log):
    """fetch passes dueIn() as its wait() timeout, so it sleeps until a line is due."""
    logger, _messages = log
    clock = FakeClock()
    progress = COMMON.BlockProgress(logger, None, interval=60, clock=clock)
    assert progress.dueIn() == 60
    clock.now += 45
    assert progress.dueIn() == 15
    clock.now += 30
    assert progress.dueIn() == 0            # overdue: wait() returns at once
    progress.updateTo(1)                    # logs, which restarts the countdown
    assert progress.dueIn() == 60


def _whole_file_range(fixture):
    with open(os.path.join(TEST_DIR, fixture + ".block.idx")) as fp:
        rows = [ln.split("\t") for ln in fp if ln.strip() and not ln.startswith("#")]
    return int(rows[0][6]), int(rows[-1][7])


@pytest.mark.parametrize("fixture,compression", [("example.maf", "none"), ("example.maf.gz", "gz")])
def test_copy_progress_hook_reports_cumulative_bytes(fixture, compression):
    start, end = _whole_file_range(fixture)
    path = os.path.join(TEST_DIR, fixture)

    plain = io.BytesIO()
    with COMMON.openMaf(path, compression, "rb") as handle:
        COMMON.copyMafRangeToStream(handle, compression, start, end, plain, chunk_size=64)

    hooked, calls = io.BytesIO(), []
    with COMMON.openMaf(path, compression, "rb") as handle:
        written = COMMON.copyMafRangeToStream(handle, compression, start, end, hooked, chunk_size=64,
                                              on_progress=calls.append)

    assert hooked.getvalue() == plain.getvalue()
    assert len(calls) == -(-written // 64)                  # once per chunk, no more
    assert calls == sorted(set(calls)) and calls[-1] == written


def test_copy_progress_hook_is_never_called_on_the_bgzip_line_loop():
    start, end = _whole_file_range("example.maf.bgz")
    calls = []
    with COMMON.openMaf(os.path.join(TEST_DIR, "example.maf.bgz"), "bgzip", "rb") as handle:
        written = COMMON.copyMafRangeToStream(handle, "bgzip", start, end, io.BytesIO(), on_progress=calls.append)
    assert written > 0 and calls == []


def _fetch(args):
    return subprocess.run([sys.executable, "-m", "mafutils", "fetch"] + args,
                          check=False, cwd=REPO_ROOT, capture_output=True, text=True)


def test_fetch_block_mode_reports_progress_setup(tmp_path):
    r = _fetch([MAF_FILE, os.path.join(TEST_DIR, "example.bed"), "-o", str(tmp_path), "-p", "2"])
    assert r.returncode == 0, r.stderr
    assert "Progress is logged every 60s." in r.stdout


@pytest.mark.parametrize("fixture", ["example.maf", "example.maf.gz", "example.maf.bgz"])
def test_fetch_scaffold_mode_reports_each_finished_scaffold(fixture, tmp_path):
    bed = os.path.join(str(tmp_path), "s.bed")
    with open(bed, "w") as fp:
        fp.write("chr1\t0\t50\tchr1\nchr4\t0\t50\tchr4\n")
    maf = os.path.join(TEST_DIR, fixture)
    r = _fetch([maf, bed, "-i", maf + ".scaffold.idx", "-m", "scaffold", "-p", "2",
                "-o", os.path.join(str(tmp_path), "out")])
    assert r.returncode == 0, r.stdout + r.stderr
    assert "of 2 scaffolds)" in r.stdout
    assert "Finished chr1" in r.stdout and "Finished chr4" in r.stdout
    assert ("counts finished scaffolds" in r.stdout) == (fixture.endswith(".bgz"))


# ---------------------------------------------------------------------------
# --progress-interval: default and minimum 60, 0 = off, anything else rejected
# ---------------------------------------------------------------------------


class CountingClock(FakeClock):
    def __init__(self):
        super().__init__()
        self.calls = 0

    def __call__(self):
        self.calls += 1
        return self.now


def test_interval_zero_turns_progress_off_without_reading_the_clock(log):
    logger, messages = log
    clock = CountingClock()
    progress = COMMON.BlockProgress(logger, 1000, interval=0, clock=clock)
    calls_at_construction = clock.calls
    for _ in range(100):
        clock.now += 3600
        progress.update(10)
        progress.updateTo(500)
    assert not progress.enabled
    assert clock.calls == calls_at_construction       # nothing ran
    assert messages == [] and progress.done == 0
    assert progress.dueIn() is None                   # wait(timeout=None): plain blocking wait


def test_progress_startup_line():
    assert COMMON.progressStartupLine(60) == "Progress is logged every 60s."
    assert COMMON.progressStartupLine(300) == "Progress is logged every 300s."
    assert COMMON.progressStartupLine(0) == "Progress reporting is off (--progress-interval 0)."


def test_header_block_count_is_quiet_about_progress_when_it_is_off(log):
    logger, messages = log
    assert COMMON.headerBlockCount({"size": "777"}, logger, progress=False) is None
    assert messages == []
    # a corrupt count is still worth a warning, just without the progress clause
    assert COMMON.headerBlockCount({"blocks": "lots"}, logger, progress=False) is None
    assert messages == ["Index header has a non-integer block count ('lots')."]


def _command_args(command, tmp_path):
    out = str(tmp_path)
    return {
        "index": [MAF_FILE, os.path.join(out, "x.block.idx"), os.path.join(out, "x.scaffold.idx")],
        "stats": [MAF_FILE, "-o", os.path.join(out, "s")],
        "gc": [MAF_FILE, "-o", os.path.join(out, "g")],
        "fetch": [MAF_FILE, os.path.join(TEST_DIR, "example.bed"), "-b", "id", "-o", os.path.join(out, "f")],
    }[command]


def _mafutils(args):
    return subprocess.run([sys.executable, "-m", "mafutils"] + args,
                          check=False, cwd=REPO_ROOT, capture_output=True, text=True)


COMMANDS = ["index", "stats", "gc", "fetch"]


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize("extra,expected", [
    ([], "Progress is logged every 60s."),
    (["--progress-interval", "120"], "Progress is logged every 120s."),
    (["--progress-interval", "0"], "Progress reporting is off (--progress-interval 0)."),
])
def test_progress_interval_values_accepted(command, extra, expected, tmp_path):
    r = _mafutils([command] + _command_args(command, tmp_path) + extra)
    assert r.returncode == 0, r.stdout + r.stderr
    assert expected in r.stdout
    assert "Progress:" not in r.stdout                 # runs this short never reach a line


@pytest.mark.parametrize("command", COMMANDS)
@pytest.mark.parametrize("value", ["30", "59", "-1"])
def test_progress_interval_below_minimum_is_a_usage_error(command, value, tmp_path):
    """Rejected, never clamped: a silently changed value would hide the mistake."""
    r = _mafutils([command] + _command_args(command, tmp_path) + ["--progress-interval", value])
    assert r.returncode == 2
    # typer draws the error in a box and wraps it; compare with the box and
    # line breaks collapsed
    flat = " ".join((r.stdout + r.stderr).replace("\u2502", " ").split())
    assert f"must be 0 (off) or at least 60 seconds, not {value}." in flat, flat


@pytest.mark.parametrize("command", COMMANDS)
def test_progress_interval_is_in_help(command):
    r = _mafutils([command, "--help"])
    assert r.returncode == 0 and "--progress-interval" in r.stdout


@pytest.mark.parametrize("command", ["info", "validate"])
def test_commands_without_progress_do_not_take_the_option(command):
    r = _mafutils([command, "--help"])
    assert r.returncode == 0 and "--progress-interval" not in r.stdout


@pytest.mark.parametrize("mode_args", [
    ["-b", "id", "-p", "2"],
    ["-m", "scaffold", "-i", MAF_FILE + ".scaffold.idx", "-p", "2"],
])
def test_fetch_outputs_identical_with_progress_off(mode_args, tmp_path):
    """Off only removes bookkeeping; the extracted alignments must not change."""
    bed = os.path.join(TEST_DIR, "example.bed")
    if "scaffold" in mode_args:
        bed = os.path.join(str(tmp_path), "s.bed")
        with open(bed, "w") as fp:
            fp.write("chr1\t0\t50\tchr1\nchr4\t0\t50\tchr4\n")
    trees = {}
    for interval in ("60", "0"):
        out = os.path.join(str(tmp_path), f"out{interval}")
        r = _fetch([MAF_FILE, bed, "-o", out, "--progress-interval", interval] + mode_args)
        assert r.returncode == 0, r.stdout + r.stderr
        trees[interval] = {
            name: open(os.path.join(out, name), "rb").read()
            for name in sorted(os.listdir(out)) if name != "maf_fetch.log"
        }
    assert trees["60"] == trees["0"] and trees["0"]


@pytest.mark.parametrize("parallel", [False, True])
def test_gc_outputs_identical_with_progress_off(parallel, tmp_path):
    results = {}
    for interval in ("60", "0"):
        prefix = os.path.join(str(tmp_path), f"g{interval}")
        args = ["gc", MAF_FILE, "-o", prefix, "--progress-interval", interval]
        if parallel:
            args[2:2] = [MAF_FILE + ".block.idx"]
            args += ["-p", "2", "--chunk-size", "3"]
        r = _mafutils(args)
        assert r.returncode == 0, r.stdout + r.stderr
        assert ("PARALLEL" in r.stdout) == parallel
        results[interval] = open(prefix + ".gc.csv").read()
    assert results["60"] == results["0"]
