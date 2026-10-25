#!/usr/bin/env python3
"""Item 6: time a production run and engine/contraction calls with measured receipts.

Tester records three repeats at each revision, profile and affinity. Instrumentation
is local to this process. FLOP/byte rates are explicit leading-work stream models,
not hardware counter measurements; missing models stay null. No golden is blessed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import resource
import time


def rates(flops, byte_count, seconds):
    if flops is None or byte_count is None or seconds <= 0:
        return {"modeled_gflops": None, "modeled_gbs": None, "one_socket_roof_gflops": None}
    intensity = flops / byte_count if byte_count else 0.0
    # Reference ceilings from the brief, not measured STREAM/compute results.
    roof = min(48 * 35.0, 190.0 * intensity)
    achieved = flops / seconds / 1e9
    return {"modeled_gflops": achieved, "modeled_gbs": byte_count / seconds / 1e9,
            "one_socket_roof_gflops": roof,
            "fraction_of_reference_one_socket_roof": achieved / roof if roof else None}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("runfile")
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--revision", required=True, help="tester-recorded git revision")
    args = parser.parse_args()
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    (out / "run.pid").write_text(str(os.getpid()) + "\n")
    from lsquant import _core
    from lsquant.runfile import RunFile, execute_runfile
    import numpy as np

    rf = RunFile.from_file(str(Path(args.runfile).resolve()))
    phases = []
    original_engine = _core.MomentEngine
    saved_functions = {}
    input_nnz = None
    if rf.model.get("kind") == "csr":
        csr = Path(rf.model["dir"]) / f"{rf.model['label']}.HAM.CSR"
        with csr.open() as handle:
            _, input_nnz = map(int, handle.readline().split())

    class MeasuredEngine:
        def __init__(self, *values, **kwargs):
            start = time.perf_counter()
            self.engine = original_engine(*values, **kwargs)
            phases.append({"phase": "engine_setup", "wall_seconds": time.perf_counter() - start})

        def run(self, plan):
            start = time.perf_counter()
            mu = self.engine.run(plan)
            seconds = time.perf_counter() - start
            n, r = int(mu.system_size), int(mu.number_of_states)
            dims = list(mu.dims)
            flops = byte_count = None
            model = "unmodeled"
            if mu.rank == 1 and input_nnz is not None:
                m = dims[0]
                flops = 8.0 * (input_nnz + n) * m * r
                byte_count = (24.0 * input_nnz + 8.0 * (n + 1) + 80.0 * n) * m * r
                model = "input CSR actions + dot products; omits added diagonal, scalar updates and setup"
            elif mu.rank == 2 and not mu.number_of_time_steps:
                m0, m1 = dims
                flops = 8.0 * n * m0 * m1 * r
                byte_count = (16.0 * n * (m0 + m1) + 16.0 * m0 * m1) * r
                model = "dot-table leading work only; unique vector/output stream lower bound; omits H/operators/packing traffic"
            raw = np.ascontiguousarray(mu.data)
            phases.append({"phase": "moments", "wall_seconds": seconds,
                           "dimension": n, "states": r, "dims": dims, "model": model,
                           "modeled_flops": flops, "modeled_bytes": byte_count,
                           "raw_moment_sha256": hashlib.sha256(memoryview(raw).cast("B")).hexdigest(),
                           **rates(flops, byte_count, seconds)})
            return mu

    def wrap_contract(name, function):
        def measured(*values, **kwargs):
            start = time.perf_counter()
            result = function(*values, **kwargs)
            seconds = time.perf_counter() - start
            mu = values[0]
            dims = list(mu.dims) if hasattr(mu, "dims") else list(mu.shape)
            flops = byte_count = None
            model = "FFT work unmodeled; wall time only"
            if name != "FftMomentContraction":
                m0, m1 = dims
                points = values[1].shape[1]
                flops = 8.0 * m0 * m1 * points
                byte_count = 16.0 * m0 * m1 + sum(x.nbytes for x in values[1:]) + result.nbytes
                model = "complex contraction leading work; unique input/output stream lower bound"
            phases.append({"phase": name, "wall_seconds": seconds, "dims": dims,
                           "modeled_flops": flops, "modeled_bytes": byte_count,
                           "model": model, **rates(flops, byte_count, seconds)})
            return result
        return measured

    _core.MomentEngine = MeasuredEngine
    for name in ("ContractRank2Imag", "ContractRank2Real", "ContractRank2RealReal", "FftMomentContraction"):
        if hasattr(_core, name):
            saved_functions[name] = getattr(_core, name)
            setattr(_core, name, wrap_contract(name, saved_functions[name]))
    before = resource.getrusage(resource.RUSAGE_SELF)
    start = time.perf_counter()
    summary, error = None, None
    try:
        summary = execute_runfile(rf, out_dir_override=str(out / "run"))
    except BaseException as exc:
        error = repr(exc)
        raise
    finally:
        seconds = time.perf_counter() - start
        after = resource.getrusage(resource.RUSAGE_SELF)
        _core.MomentEngine = original_engine
        for name, function in saved_functions.items():
            setattr(_core, name, function)
        env_names = ("OMP_NUM_THREADS", "OMP_PROC_BIND", "OMP_PLACES", "OMP_DYNAMIC",
                     "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "LSQUANT_SPMV_MIN_NNZ",
                     "LSQUANT_KERNEL_BATCH", "LSQUANT_FFT_BACKEND", "KPM_SEED")
        report = {"revision": args.revision, "build_profile": getattr(_core, "build_profile", "golden"),
                  "runfile": str(Path(args.runfile).resolve()), "host": platform.node(),
                  "python": platform.python_version(), "numpy": np.__version__,
                  "affinity": sorted(os.sched_getaffinity(0)),
                  "environment": {key: os.environ.get(key) for key in env_names},
                  "wall_total_seconds": seconds,
                  "cpu_seconds": after.ru_utime + after.ru_stime - before.ru_utime - before.ru_stime,
                  "max_rss_kib": after.ru_maxrss, "phases": phases,
                  "summary": summary, "error": error,
                  "timing_note": "imports excluded; total includes IO/setup/reconstruction/hash recording; phase timers exclude hashing"}
        (out / "benchmark.json").write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
