#!/usr/bin/env python3
"""Run rigid-contact benchmark executables once they are built."""

import argparse
import pathlib
import subprocess
import sys


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", default="Release")
    args = parser.parse_args()

    root = pathlib.Path(__file__).resolve().parents[1]
    exe_dir = root / "build" / args.config
    benchmarks = [
        exe_dir / "run_sphere_plane.exe",
        exe_dir / "run_sphere_sphere.exe",
        exe_dir / "run_benchmarks.exe",
    ]

    missing = [str(p) for p in benchmarks if not p.exists()]
    if missing:
        print("Benchmark executables are not built yet:")
        for path in missing:
            print(f"  {path}")
        return 2

    for exe in benchmarks:
        subprocess.run([str(exe)], cwd=root, check=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())

