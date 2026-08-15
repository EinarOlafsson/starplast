#!/usr/bin/env python3
"""Time identical repeated clusterings on Starplast's CPU and GPU paths.

The benchmark uses an 8,140 x 3 fixture, matching the shipped proteome map, warms each backend once,
and reports GPU utilization sampled from ``nvidia-smi`` while the timed repetitions run. It never
touches saved searches, so the CPU baseline under ``~/.cache/starplast/searches`` remains intact.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import threading
import time

import numpy as np


def _utilization(stop: threading.Event, samples: list) -> None:
    """Sample device utilization until the timed section ends."""
    while not stop.is_set():
        try:
            out = subprocess.check_output(
                ["nvidia-smi", "--query-gpu=utilization.gpu", "--format=csv,noheader,nounits"],
                text=True, timeout=2)
            samples.extend(int(line.strip()) for line in out.splitlines() if line.strip())
        except (OSError, subprocess.SubprocessError, ValueError):
            return
        stop.wait(0.2)


def run(algorithm: str, repeats: int, use_gpu: bool, seed: int = 42) -> dict:
    """Run one warmed benchmark and return its timing and utilization record."""
    os.environ["STARPLAST_GPU"] = "1" if use_gpu else "0"
    from starplast import clustering, gpu
    X = np.random.default_rng(seed).normal(size=(8140, 3)).astype(np.float32)
    options = ({"n_clusters": 60} if algorithm == "kmeans"
               else {"eps": 0.25, "min_samples": 10})
    clustering.cluster(X, algorithm=algorithm, **options)
    stop, samples = threading.Event(), []
    watcher = threading.Thread(target=_utilization, args=(stop, samples), daemon=True)
    watcher.start()
    started = time.perf_counter()
    for _ in range(int(repeats)):
        clustering.cluster(X, algorithm=algorithm, **options)
    seconds = time.perf_counter() - started
    stop.set()
    watcher.join(timeout=3)
    return {"algorithm": algorithm, "backend": gpu.backend()[algorithm], "repeats": int(repeats),
            "seconds": seconds, "per_clustering": seconds / max(int(repeats), 1),
            "gpu_util_mean": float(np.mean(samples)) if samples else np.nan,
            "gpu_util_peak": int(max(samples)) if samples else None}


def main(argv=None) -> int:
    """Command-line entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("algorithm", choices=("kmeans", "dbscan"))
    parser.add_argument("--repeats", type=int, default=100)
    parser.add_argument("--gpu", action="store_true")
    args = parser.parse_args(argv)
    result = run(args.algorithm, args.repeats, args.gpu)
    print("  ".join(f"{key}={value:.4f}" if isinstance(value, float) else f"{key}={value}"
                    for key, value in result.items()))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
