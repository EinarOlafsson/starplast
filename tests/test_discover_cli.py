#!/usr/bin/env python3
"""The search from a terminal.

Two properties matter more than the argument parsing. The first: nothing in this path imports Qt,
because a headless option that pulls in a GUI toolkit fails on the machine it exists for. The
second: a batch survives one task failing, and survives being killed -- the failure this command
was written for is a long run dying two thirds of the way through.
"""
from __future__ import annotations

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import discover as D  # noqa: E402


@pytest.fixture
def table(tmp_path):
    """A small node table on disk, with enough structure to find something in."""
    rng = np.random.default_rng(0)
    n = 240
    rows = []
    for i in range(n):
        group = i % 4
        rows.append({"gene_id": f"G{i:04d}", "product": "hypothetical protein",
                     "n_publications": 0,
                     "compartment": None if i % 9 == 0 else ["IMC", "ER", "nucleus", "rhoptry"][group],
                     "fitness": float(rng.normal(-3 if group == 0 else 0, 0.4)),
                     **{f"expr_{j}": float(rng.normal(group * 4, 0.6)) for j in range(5)}})
    nodes = pd.DataFrame(rows)
    where = str(tmp_path / "nodes.parquet")
    nodes.to_parquet(where)
    return where, str(tmp_path / "searches")


# --------------------------------------------------------------------------- the arguments
def test_a_task_is_mode_layer_and_optionally_what_it_disagrees_with():
    assert D.parse_task("guilt:compartment") == {"mode": "guilt", "layer": "compartment",
                                                 "against": ""}
    assert D.parse_task("disagreement:a:b")["against"] == "b"


@pytest.mark.parametrize("bad", ["nonsense", "wrongmode:x", ":x", "guilt:", ""])
def test_a_task_that_is_not_a_task_says_what_was_wrong(bad):
    with pytest.raises(ValueError) as caught:
        D.parse_task(bad)
    assert "mode:layer" in str(caught.value) or "unknown mode" in str(caught.value)


def test_nothing_to_do_is_an_error_rather_than_a_silent_success(capsys):
    with pytest.raises(SystemExit):
        D.main([])
    assert "nothing to do" in capsys.readouterr().err


def test_a_bad_task_is_refused_before_anything_expensive_starts(capsys):
    assert D.main(["--task", "telepathy:compartment"]) == 2
    assert "unknown mode" in capsys.readouterr().err


# --------------------------------------------------------------------------- running
def test_a_search_runs_and_is_saved_where_the_window_will_find_it(table, capsys):
    nodes, out = table
    code = D.main(["--task", "guilt:compartment", "--budget", "6", "--restarts", "1",
                   "--nodes", nodes, "--out", out, "--name", "mine", "--quiet"])
    assert code == 0
    printed = capsys.readouterr().out
    assert "mine:" in printed and "configs" in printed
    from starplast.searches import SearchStore
    saved = SearchStore(out).load("mine")
    assert not saved.configs.empty
    assert saved.manifest["mode"] == "guilt" and saved.manifest["layer"] == "compartment"
    assert saved.matches(pd.read_parquet(nodes))


def test_a_task_already_saved_is_skipped_so_a_killed_batch_can_be_rerun(table, capsys):
    nodes, out = table
    args = ["--task", "guilt:compartment", "--budget", "4", "--restarts", "1",
            "--nodes", nodes, "--out", out, "--name", "once"]
    D.main(args)
    capsys.readouterr()
    D.main(args)
    assert "already saved" in capsys.readouterr().out


def test_one_failing_task_does_not_take_the_batch_with_it(table, capsys):
    """A batch that dies on its third task and loses the other seven is the failure this exists to
    avoid."""
    nodes, out = table
    code = D.main(["--task", "guilt:not_a_column", "--task", "guilt:compartment",
                   "--budget", "4", "--restarts", "1", "--nodes", nodes, "--out", out,
                   "--prefix", "batch_", "--quiet"])
    seen = capsys.readouterr()
    assert code == 1, "a failed task did not show in the exit code"
    assert "FAILED" in seen.err
    assert "batch_01_guilt_compartment" in seen.out, "the second task did not run"


def test_a_search_over_no_blocks_at_all_is_refused(table):
    nodes, out = table
    with pytest.raises(ValueError):
        D.run_task(pd.read_parquet(nodes), D.parse_task("guilt:compartment"),
                   blocks=None, exclude=["expression_summary", "literature"],
                   store=D.store_for(out), log=lambda *_a: None)


# --------------------------------------------------------------------------- reading back
def test_saved_searches_can_be_listed_and_read_without_a_window(table, capsys):
    nodes, out = table
    D.main(["--task", "guilt:compartment", "--budget", "6", "--restarts", "1",
            "--nodes", nodes, "--out", out, "--name", "readable", "--quiet"])
    capsys.readouterr()
    assert D.main(["--list", "--out", out]) == 0
    assert "readable" in capsys.readouterr().out
    assert D.main(["--read", "readable", "--out", out, "--nodes", nodes]) == 0
    text = capsys.readouterr().out
    assert "claims survive correction" in text or "Nothing to report" in text


def test_reading_a_search_that_is_not_there_says_so(table, capsys):
    _nodes, out = table
    assert D.main(["--read", "never_ran", "--out", out]) == 1
    assert "no search called" in capsys.readouterr().err


def test_an_empty_store_lists_nothing_rather_than_failing(tmp_path, capsys):
    assert D.main(["--list", "--out", str(tmp_path)]) == 0
    assert "no saved searches" in capsys.readouterr().out


def test_a_missing_node_table_names_the_command_that_builds_it(tmp_path):
    with pytest.raises(SystemExit) as caught:
        D.load_nodes(str(tmp_path / "absent.parquet"))
    assert "build_graph" in str(caught.value)


def test_the_summary_line_survives_a_search_that_found_nothing():
    from starplast.searches import Search
    assert D._summary(None) == "nothing evaluated"
    assert D._summary(Search()) == "nothing evaluated"


# --------------------------------------------------------------------------- and no window
def test_the_command_does_not_import_a_gui_toolkit(table):
    """The whole point. Checked in a fresh interpreter, because this test session has Qt loaded
    already and asking `sys.modules` here would prove nothing."""
    nodes, out = table
    code = ("import sys; from starplast import discover; "
            "discover.main(['--list', '--out', %r]); "
            "print('QT' if any(m.startswith('PyQt') for m in sys.modules) else 'CLEAN')" % out)
    done = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True,
                          cwd=os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    assert done.returncode == 0, done.stderr
    assert done.stdout.strip().endswith("CLEAN"), done.stdout
