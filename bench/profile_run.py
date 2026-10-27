#!/usr/bin/env python3
"""Item 1: profile a run in a child; measure Python/C++ time and CPU occupancy.

The tester compares top-20 cumulative entries, attributed C++ self-time fraction,
thread count and ps CPU percent before/after each item on the same Aron inputs.
perf is a separate unprofiled replay: profiling overhead must not enter its counts.
No model files are changed; --moments/--states override only the in-memory run.
"""
from __future__ import annotations

import argparse
import cProfile
import io
import json
import os
from pathlib import Path
import pstats
import shutil
import subprocess
import sys
import time


def child(args):
    from lsquant.runfile import RunFile, execute_runfile

    run = RunFile.from_file(str(Path(args.runfile).resolve()))
    if args.moments is not None:
        run.expansion["moments"] = args.moments
        run.expansion["batch"] = None
    if args.states is not None:
        run.expansion["states"] = args.states
    profiler = cProfile.Profile()
    started = time.perf_counter()
    try:
        if not args.unprofiled:
            profiler.enable()
        execute_runfile(run, out_dir_override=str(Path(args.out) / "run"))
    finally:
        profiler.disable()
        wall = time.perf_counter() - started
        if not args.unprofiled:
            profiler.dump_stats(str(Path(args.out) / "profile.pstats"))
            stream = io.StringIO()
            stats = pstats.Stats(profiler, stream=stream).sort_stats("cumulative")
            stats.print_stats(20)
            # pybind11 exports appear as PyCapsule builtins in CPython's profiler.
            # Self time avoids double-counting nested C++/Python callback frames.
            entries = {str(key): val[2] for key, val in stats.stats.items()
                       if "_core" in str(key) or "PyCapsule" in str(key)}
            report = {"profiled_wall_seconds": wall,
                      "profile_total_seconds": stats.total_tt,
                      "core_attributed_self_seconds": sum(entries.values()),
                      "core_attributed_fraction": (sum(entries.values()) / stats.total_tt
                                                   if stats.total_tt else 0.0),
                      "core_attributed_entries": entries}
            (Path(args.out) / "profile.txt").write_text(stream.getvalue())
            (Path(args.out) / "profile.json").write_text(json.dumps(report, indent=2) + "\n")
            print(stream.getvalue(), flush=True)
            print(json.dumps(report, indent=2), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runfile")
    parser.add_argument("--out", required=True)
    parser.add_argument("--moments", type=int)
    parser.add_argument("--states", type=int)
    parser.add_argument("--interval", type=float, default=1.0)
    parser.add_argument("--perf", action="store_true", help="also replay without cProfile")
    parser.add_argument("--child", action="store_true", help=argparse.SUPPRESS)
    parser.add_argument("--unprofiled", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.interval <= 0:
        parser.error("--interval must be positive")
    out = Path(args.out).resolve()
    out.mkdir(parents=True, exist_ok=True)
    args.out = str(out)
    if args.child:
        child(args)
        return 0
    command = [sys.executable, str(Path(__file__).resolve()), args.runfile,
               "--out", str(out), "--child"]
    for name in ("moments", "states"):
        if getattr(args, name) is not None:
            command += ["--" + name, str(getattr(args, name))]
    with (out / "run.log").open("w") as log, (out / "occupancy.jsonl").open("w") as samples:
        process = subprocess.Popen(command, stdout=log, stderr=subprocess.STDOUT)
        (out / "run.pid").write_text(str(process.pid) + "\n")
        started = time.monotonic()
        try:
            while process.poll() is None:
                status_path = Path(f"/proc/{process.pid}/status")
                try:
                    status = dict(line.split(":", 1) for line in status_path.read_text().splitlines()
                                  if ":" in line)
                    cpu = subprocess.run(["ps", "-p", str(process.pid), "-o", "%cpu="],
                                         capture_output=True, text=True, check=False)
                    sample = {"elapsed_seconds": time.monotonic() - started,
                              "threads": int(status.get("Threads", "0").strip()),
                              "cpu_percent": float(cpu.stdout.strip() or "0"),
                              "allowed_cpus": status.get("Cpus_allowed_list", "").strip(),
                              "allowed_memory_nodes": status.get("Mems_allowed_list", "").strip()}
                    samples.write(json.dumps(sample) + "\n")
                    samples.flush()
                except (FileNotFoundError, ProcessLookupError):
                    pass
                time.sleep(args.interval)
        except KeyboardInterrupt:
            process.terminate()  # Only the PID this launcher owns.
        code = process.wait()
    print(f"profile exit={code}; artifacts: {out}")
    if code:
        return code
    if args.perf:
        if shutil.which("perf"):
            replay = command + ["--unprofiled"]
            with (out / "perf-run.log").open("w") as log:
                perf = subprocess.Popen(["perf", "stat", "-x", ",", "-o", str(out / "perf.csv"),
                                         "-e", "task-clock,cycles,instructions,cache-misses",
                                         "--", *replay], stdout=log, stderr=subprocess.STDOUT)
                (out / "perf.pid").write_text(str(perf.pid) + "\n")
                perf_code = perf.wait()
            (out / "perf-status.txt").write_text(f"exit={perf_code}\n")
        else:
            (out / "perf-status.txt").write_text("perf unavailable\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
