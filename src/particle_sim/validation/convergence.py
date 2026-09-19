"""Scenario-count convergence check for probability CA.

The tolerance is not fixed in code: the caller decides it from pilot results.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

CONVERGENCE_COLUMNS: tuple[str, ...] = (
    "scramble_seed",
    "source_step",
    "particle_size_um",
    "pair_id",
    "failure_mode",
    "scenario_count",
    "ca_um2",
    "delta_previous",
)


@dataclass(frozen=True)
class ConvergenceStep:
    scenario_count: int
    ca_um2: float
    delta_previous: float  # |CA_n - CA_prev| / CA_n, NaN for the first entry


def convergence_series(samples: dict[int, float]) -> list[ConvergenceStep]:
    """Relative change of CA between successive scenario counts (64 -> 128 -> 256 ...)."""
    if not samples:
        raise ValueError("수렴 검사에는 최소 한 개의 결과가 필요합니다")
    steps: list[ConvergenceStep] = []
    previous: float | None = None
    for count in sorted(samples):
        current = float(samples[count])
        if previous is None:
            delta = float("nan")
        elif current == 0.0:
            delta = 0.0 if previous == 0.0 else float("inf")
        else:
            delta = abs(current - previous) / abs(current)
        steps.append(
            ConvergenceStep(scenario_count=int(count), ca_um2=current, delta_previous=delta)
        )
        previous = current
    return steps


def is_converged(steps: list[ConvergenceStep], tolerance: float) -> bool:
    """True when the last comparison is within `tolerance` (caller-provided, never default)."""
    if tolerance < 0:
        raise ValueError(f"허용 오차는 0 이상이어야 합니다: {tolerance}")
    comparisons = [s.delta_previous for s in steps[1:]]
    if not comparisons:
        return False
    last = comparisons[-1]
    if not np.isfinite(last):
        # an infinite or undefined last comparison is never evidence of convergence
        return False
    return bool(last <= tolerance)
