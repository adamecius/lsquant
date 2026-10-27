"""CPU affinity and kernel tuning defaults applied before loading OpenMP or NumPy."""
from __future__ import annotations

import os
from pathlib import Path

DEFAULT_SPMV_MIN_NNZ = "16384"
DEFAULT_KERNEL_BATCH = "0"


def physical_core_count() -> int:
    """Count distinct physical cores in the process CPU affinity, including cpusets."""
    try:
        allowed = os.sched_getaffinity(0)
    except AttributeError:
        return os.cpu_count() or 1
    cores = set()
    for cpu in allowed:
        topology = Path(f"/sys/devices/system/cpu/cpu{cpu}/topology")
        try:
            cores.add(((topology / "physical_package_id").read_text().strip(),
                       (topology / "core_id").read_text().strip()))
        except OSError:
            cores.add(("cpu", str(cpu)))
    return max(1, len(cores))


def configure_threads() -> None:
    """Respect caller overrides; close binding keeps static work near its first touch."""
    os.environ.setdefault("OMP_PROC_BIND", "close")
    os.environ.setdefault("OMP_PLACES", "cores")
    if "OMP_NUM_THREADS" not in os.environ:
        os.environ["OMP_NUM_THREADS"] = str(physical_core_count())


def configure_kernels(core) -> None:
    """Send numeric tuning settings once; the engine never reads the environment."""
    if core.build_profile == "native":
        core.ConfigurePerformance(
            int(os.environ.get("LSQUANT_SPMV_MIN_NNZ", DEFAULT_SPMV_MIN_NNZ)),
            int(os.environ.get("LSQUANT_KERNEL_BATCH", DEFAULT_KERNEL_BATCH)),
        )
