#!/usr/bin/env python3
"""Choosing a CUDA wheel set for the machine this is running on.

The rule under test, and the reason the command exists at all: **the newest set the driver could run
is not the right answer.** NVIDIA drivers are backward compatible, so a driver serving CUDA 13 runs
CUDA 12 wheels perfectly, while cu13 wheels need a CUDA 13 driver and are the newer build. Getting
this wrong costs two gigabytes and a broken environment, which is why nothing is downloaded before
the plan is printed.
"""
from __future__ import annotations

import os
import sys

import pytest

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from starplast import install_gpu as IG  # noqa: E402

SMI = """Sat Aug 13 09:00:00 2026
+---------------------------------+
| NVIDIA-SMI 580.17  Driver Version: 580.173.02  CUDA Version: 13.0 |
+---------------------------------+
"""


def _driver(monkeypatch, text, exe="/usr/bin/nvidia-smi"):
    import subprocess
    monkeypatch.setattr(IG.shutil, "which", lambda name: exe)
    monkeypatch.setattr(IG.subprocess, "run",
                        lambda *a, **k: subprocess.CompletedProcess(a, 0, stdout=text, stderr=""))


def test_the_driver_decides_not_the_toolkit(monkeypatch):
    """nvidia-smi reports the highest CUDA the DRIVER supports; nvcc reports whichever toolkit is
    installed, which is frequently older or absent. The wheels load against the driver."""
    _driver(monkeypatch, SMI)
    assert IG.driver_cuda() == (13, "driver serves CUDA 13.0")


def test_no_card_is_reported_as_no_card(monkeypatch):
    monkeypatch.setattr(IG.shutil, "which", lambda name: None)
    major, why = IG.driver_cuda()
    assert major == 0 and "no NVIDIA driver" in why
    p = IG.plan()
    assert p["packages"] == [] and "runs without any of this" in p["note"]


def test_a_driver_that_says_nothing_useful_is_not_guessed_at(monkeypatch):
    _driver(monkeypatch, "some other output entirely")
    assert IG.driver_cuda()[0] == 0
    assert IG.plan()["packages"] == []


def test_the_set_matching_torch_wins(monkeypatch):
    """Mixing CUDA majors in one environment is how two libraries load two runtimes. This machine
    has torch built for cu12, so cuml for cu12 is what belongs beside it."""
    _driver(monkeypatch, SMI)
    monkeypatch.setattr(IG, "torch_cuda", lambda: 12)
    p = IG.plan()
    assert p["cuda"] == 12 and any("cu12" in w for w in p["packages"])
    assert "matching torch" in p["why"]


def test_without_torch_the_better_travelled_build_wins(monkeypatch):
    """Not the newest the driver could run: cu12 wheels run on a CUDA 13 driver, and cu13 wheels do
    not run on a CUDA 12 one."""
    _driver(monkeypatch, SMI)
    monkeypatch.setattr(IG, "torch_cuda", lambda: 0)
    p = IG.plan()
    assert p["cuda"] == 12
    assert "better-travelled" in p["why"]


def test_an_old_driver_gets_no_wheels_rather_than_wrong_ones(monkeypatch):
    _driver(monkeypatch, "CUDA Version: 11.4")
    monkeypatch.setattr(IG, "torch_cuda", lambda: 0)
    p = IG.plan()
    assert p["packages"] == [] and "older than any wheel set" in p["why"]


def test_the_choice_can_be_overridden(monkeypatch):
    _driver(monkeypatch, SMI)
    monkeypatch.setattr(IG, "torch_cuda", lambda: 12)
    assert IG.plan(prefer=13)["cuda"] == 13
    unknown = IG.plan(prefer=11)
    assert unknown["packages"] == [] and "no wheel set" in unknown["why"]


def test_torch_reports_the_cuda_it_was_built_for():
    """Reads torch.version.cuda where torch is installed, and 0 where it is not -- never raises,
    because this runs before anything is known about the machine."""
    assert isinstance(IG.torch_cuda(), int)


def test_nothing_is_downloaded_before_the_plan_is_shown(monkeypatch, capsys):
    """Two gigabytes is not something to start on a user's behalf and explain afterwards."""
    _driver(monkeypatch, SMI)
    monkeypatch.setattr(IG, "torch_cuda", lambda: 12)
    ran = []
    monkeypatch.setattr(IG.subprocess, "call", lambda cmd: ran.append(cmd) or 0)
    assert IG.main(["--dry-run"]) == 0
    out = capsys.readouterr().out
    assert "would install" in out and "pip install" in out
    assert ran == [], "a dry run downloaded something"


def test_it_asks_before_installing(monkeypatch, capsys):
    _driver(monkeypatch, SMI)
    monkeypatch.setattr(IG, "torch_cuda", lambda: 12)
    ran = []
    monkeypatch.setattr(IG.subprocess, "call", lambda cmd: ran.append(cmd) or 0)
    monkeypatch.setattr("builtins.input", lambda *_: "n")
    assert IG.main([]) == 0
    assert ran == [] and "nothing installed" in capsys.readouterr().out

    monkeypatch.setattr("builtins.input", lambda *_: "y")
    assert IG.main([]) == 0
    assert ran and "cuml-cu12>=24.10" in ran[0]


def test_yes_skips_the_question(monkeypatch, capsys):
    _driver(monkeypatch, SMI)
    monkeypatch.setattr(IG, "torch_cuda", lambda: 12)
    monkeypatch.setattr(IG.subprocess, "call", lambda cmd: 0)
    monkeypatch.setattr("builtins.input", lambda *_: pytest.fail("it asked anyway"))
    assert IG.main(["--yes"]) == 0
    assert "restart starplast" in capsys.readouterr().out


def test_a_failed_pip_is_reported_with_its_code(monkeypatch, capsys):
    _driver(monkeypatch, SMI)
    monkeypatch.setattr(IG, "torch_cuda", lambda: 12)
    monkeypatch.setattr(IG.subprocess, "call", lambda cmd: 1)
    assert IG.main(["--yes"]) == 1
    assert "pip exited 1" in capsys.readouterr().out


def test_with_no_card_it_exits_nonzero(monkeypatch, capsys):
    monkeypatch.setattr(IG.shutil, "which", lambda name: None)
    assert IG.main(["--yes"]) == 1
    assert "runs without any of this" in capsys.readouterr().out


def test_a_malformed_cuda_flag_is_refused(monkeypatch, capsys):
    _driver(monkeypatch, SMI)
    assert IG.main(["--cuda"]) == 2
    assert IG.main(["--cuda", "twelve"]) == 2
    assert "takes a major version" in capsys.readouterr().out


def test_the_command_is_this_interpreter_not_whichever_pip_is_on_path():
    """`pip` on PATH can belong to another environment entirely, and installing cuml into the wrong
    one is a confusing hour."""
    cmd = IG.command(["cuml-cu12"])
    assert cmd[0] == sys.executable and cmd[1:4] == ["-m", "pip", "install"]


def test_an_nvidia_smi_that_will_not_run_is_not_fatal(monkeypatch):
    """It is on PATH and it fails -- a driver mid-upgrade, a container without /dev/nvidia. This
    command exists to tell someone what they can install, so it must survive being unable to."""
    monkeypatch.setattr(IG.shutil, "which", lambda name: "/usr/bin/nvidia-smi")

    def boom(*a, **k):
        raise OSError("no such device")

    monkeypatch.setattr(IG.subprocess, "run", boom)
    major, why = IG.driver_cuda()
    assert major == 0 and "would not run" in why and "OSError" in why


def test_torch_without_cuda_reports_nothing(monkeypatch):
    """A CPU-only torch build has version.cuda of None, which is not a version."""
    import types
    fake = types.ModuleType("torch")
    fake.version = types.SimpleNamespace(cuda=None)
    monkeypatch.setitem(sys.modules, "torch", fake)
    assert IG.torch_cuda() == 0


def test_a_torch_that_cannot_answer_is_not_a_failure(monkeypatch):
    """Half-installed torch raises on attribute access rather than importing cleanly, and this runs
    before anything is known about the machine."""
    import types
    monkeypatch.setitem(sys.modules, "torch", types.ModuleType("torch"))   # no .version at all
    assert IG.torch_cuda() == 0
