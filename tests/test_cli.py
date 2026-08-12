import os
import re
import subprocess
import sys

import pytest

TEST_DIR = os.path.dirname(__file__)
REPO_ROOT = os.path.abspath(os.path.join(TEST_DIR, ".."))


def run(args):
    return subprocess.run(
        [sys.executable, "-m", "mafutils"] + args,
        check=False, cwd=REPO_ROOT, capture_output=True, text=True,
    )


@pytest.mark.parametrize("flag", ["--version", "-version", "-v", "-V"])
def test_version_flag_prints_version_and_exits_cleanly(flag):
    result = run([flag])
    assert result.returncode == 0, result.stderr
    # e.g. "mafutils 0.3.0" / "mafutils 0.3.1.dev3+g1234abc" / "mafutils unknown"
    assert re.match(r"^mafutils \S+$", result.stdout.strip()), result.stdout


def test_version_flag_works_without_a_subcommand():
    """
    --version is is_eager so it short-circuits before subcommand parsing --
    without that, `mafutils --version` would fail asking for a COMMAND
    (the app is built with no_args_is_help=True).
    """
    result = run(["--version"])
    assert result.returncode == 0
    assert "Usage:" not in result.stdout


def test_single_dash_version_is_not_parsed_as_grouped_short_flags():
    """
    `-version` is declared alongside the short flag `-v`, so a parser could
    plausibly read it as grouped short options (-v -e -r -s -i -o -n) and
    either error or silently do the wrong thing. Assert it's treated as the
    version alias it's declared as.
    """
    result = run(["-version"])
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip().startswith("mafutils ")
    assert "no such option" not in (result.stdout + result.stderr).lower()


def test_top_level_callback_does_not_break_subcommands():
    """
    Guard against the @app.callback() added for --version interfering with
    normal subcommand dispatch.
    """
    result = run(["index", "--help"])
    assert result.returncode == 0, result.stderr
    assert "block index" in result.stdout.lower()


def test_package_exposes_dunder_version():
    import mafutils

    assert isinstance(mafutils.__version__, str)
    assert mafutils.__version__
