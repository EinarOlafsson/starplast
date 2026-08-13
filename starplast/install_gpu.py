#!/usr/bin/env python3
"""Installing the GPU stack for the machine this is actually running on.

Why this is a run-time command and not something the installer does: a wheel's dependencies are
fixed when the wheel is BUILT. pip resolves them from static metadata, and the only conditionals
available are environment markers -- `sys_platform`, `python_version`, `platform_machine`. There is
no marker for "this machine has an NVIDIA card". Probing at install time would mean shipping an
sdist with executable setup.py, which loses the wheel for the common case, breaks `--only-binary`,
and makes one version resolve to different dependency trees on different machines, so lockfiles stop
meaning anything. In CI or a container it would probe the build machine rather than the target.

There is a second reason, which is that finding a card does not tell you which CUDA to install. The
driver does. Guess wrong and it is two gigabytes of wheels that cannot load.

So detection happens here, when the real machine is known, and the user is told what will be
installed before anything is downloaded.
"""
from __future__ import annotations

import re
import shutil
import subprocess
import sys

#: What each CUDA major version needs. RAPIDS and CuPy publish to PyPI proper, so no extra index.
WHEELS = {
    12: ["cuml-cu12>=24.10", "cupy-cuda12x>=13.0"],
    13: ["cuml-cu13", "cupy-cuda13x"],
}


def driver_cuda() -> tuple:
    """(major, detail) for the CUDA the driver serves, or (0, why not).

    Read from `nvidia-smi`, which reports the highest CUDA a driver supports -- not from `nvcc`,
    which reports whatever toolkit happens to be installed and is frequently older or absent. The
    wheels are built against a CUDA major version and load against the driver, so the driver is the
    one that decides.
    """
    exe = shutil.which("nvidia-smi")
    if exe is None:
        return 0, "nvidia-smi is not on PATH, so there is no NVIDIA driver to build against"
    try:
        out = subprocess.run([exe], capture_output=True, text=True, timeout=30).stdout
    except Exception as exc:                    # pragma: no cover - defensive
        return 0, f"nvidia-smi would not run: {type(exc).__name__}: {exc}"
    m = re.search(r"CUDA Version:\s*(\d+)\.(\d+)", out)
    if not m:
        return 0, "nvidia-smi ran but reported no CUDA version"
    return int(m.group(1)), f"driver serves CUDA {m.group(1)}.{m.group(2)}"


def torch_cuda() -> int:
    """The CUDA major version torch was built against, or 0 -- another vote on which set to install.

    Worth asking because mixing CUDA majors in one environment is how two libraries end up loading
    two runtimes: this machine has torch built for cu12, and installing cuml for cu13 beside it is
    asking for it.
    """
    try:
        import torch
        m = re.match(r"(\d+)\.", torch.version.cuda or "")
        return int(m.group(1)) if m else 0
    except Exception:
        return 0


def plan(prefer: int = 0) -> dict:
    """What would be installed, and why -- decided before anything is downloaded.

    **The newest set the driver could run is not the right answer.** NVIDIA drivers are backward
    compatible, so a driver serving CUDA 13 runs CUDA 12 wheels perfectly, while cu13 wheels need a
    CUDA 13 driver and are the newer, less-travelled build. So the rule is: take what torch in this
    environment was built against if it is installed and supported, else the lowest set the driver
    can run -- which is cu12 on everything current. `--cuda N` overrides, for a machine that wants
    the newer set on purpose.
    """
    major, detail = driver_cuda()
    if not major:
        return {"packages": [], "cuda": 0, "why": detail,
                "note": "starplast runs without any of this; the GPU switch simply has nothing to "
                        "turn on."}
    runnable = sorted(k for k in WHEELS if k <= major)
    if not runnable:
        return {"packages": [], "cuda": major,
                "why": f"{detail}, which is older than any wheel set here",
                "note": f"known sets: {', '.join(f'CUDA {k}' for k in sorted(WHEELS))}."}
    if prefer:
        if prefer not in WHEELS:
            return {"packages": [], "cuda": prefer, "why": f"no wheel set for CUDA {prefer}",
                    "note": f"known sets: {', '.join(f'CUDA {k}' for k in sorted(WHEELS))}."}
        chosen, why = prefer, f"{detail}; CUDA {prefer} asked for"
    elif torch_cuda() in runnable:
        chosen = torch_cuda()
        why = f"{detail}; matching torch, which is built for CUDA {chosen}"
    else:
        chosen = runnable[0]
        why = (f"{detail}; taking CUDA {chosen}, which that driver runs and which is the "
               f"better-travelled build")
    return {"packages": WHEELS[chosen], "cuda": chosen, "why": why,
            "note": "about 2 GB. cuml moves UMAP and HDBSCAN to the GPU; cupy moves the array work."}


def command(packages) -> list:
    """The exact pip command, so it can be shown before it is run."""
    return [sys.executable, "-m", "pip", "install", *packages]


def main(argv=None) -> int:
    """`starplast-install-gpu`: report the plan, then run it unless asked not to."""
    argv = list(sys.argv[1:] if argv is None else argv)
    dry = "--dry-run" in argv or "-n" in argv
    yes = "--yes" in argv or "-y" in argv
    prefer = 0
    if "--cuda" in argv:
        try:
            prefer = int(argv[argv.index("--cuda") + 1])
        except (IndexError, ValueError):
            print("--cuda takes a major version, e.g. --cuda 12")
            return 2
    p = plan(prefer)
    print(p["why"])
    if not p["packages"]:
        print(p["note"])
        return 1
    print(f"would install: {' '.join(p['packages'])}")
    print(p["note"])
    print("  " + " ".join(command(p["packages"])))
    if dry:
        return 0
    if not yes:
        try:
            if input("run it? [y/N] ").strip().lower() not in ("y", "yes"):
                print("nothing installed")
                return 0
        except EOFError:                        # pragma: no cover - no terminal to ask on
            print("nothing installed (no terminal to ask on; pass --yes)")
            return 0
    code = subprocess.call(command(p["packages"]))
    print("installed -- restart starplast" if code == 0 else f"pip exited {code}")
    return code


if __name__ == "__main__":                      # pragma: no cover - the console script calls main
    sys.exit(main())
