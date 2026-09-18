"""Progress reporting and resource measurement for long runs.

A closed network run can take minutes, so the CLI shows progress on a terminal and every run
manifest records the peak resident memory of the process.
"""

from __future__ import annotations

import os
import resource
import sys
from collections.abc import Iterable, Iterator
from typing import TypeVar

T = TypeVar("T")


def track(items: Iterable[T], *, total: int | None = None, description: str = "") -> Iterator[T]:
    """Yield `items`, showing a progress bar only when stderr is an interactive terminal."""
    disable = not sys.stderr.isatty() or os.environ.get("PARTICLE_SIM_NO_PROGRESS") == "1"
    try:
        from tqdm import tqdm
    except ImportError:  # progress is optional, never a hard dependency at runtime
        yield from items
        return
    yield from tqdm(items, total=total, desc=description, disable=disable, file=sys.stderr)


def peak_memory_mb() -> float:
    """Peak resident set size of this process in MiB."""
    try:
        import psutil

        process = psutil.Process()
        info = process.memory_info()
        peak = getattr(info, "peak_wset", None)  # Windows
        if peak:
            return float(peak) / (1024 * 1024)
    except Exception:  # psutil is optional; fall back to the stdlib
        pass
    usage = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    # Linux reports KiB, macOS reports bytes
    divisor = 1024 * 1024 if sys.platform == "darwin" else 1024
    return float(usage) / divisor
