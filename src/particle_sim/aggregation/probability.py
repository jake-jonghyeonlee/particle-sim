"""Probability hotspots accumulated over process-variation scenarios.

The value is a geometric contact risk conditioned on the scenario model, not an actual
failure probability.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

from particle_sim.aggregation.critical_area import critical_area_um2
from particle_sim.io.masks import GridSpec


@dataclass
class ProbabilityAccumulator:
    """Counts failures per pixel over scenarios and keeps the per-scenario CA values."""

    shape_yx: tuple[int, int]
    grid: GridSpec
    failure_count: np.ndarray = field(init=False)
    margin_sum: np.ndarray = field(init=False)
    margin_samples: int = field(init=False, default=0)
    scenario_ca_um2: list[float] = field(init=False, default_factory=list)

    def __post_init__(self) -> None:
        self.failure_count = np.zeros(self.shape_yx, dtype=np.uint32)
        self.margin_sum = np.zeros(self.shape_yx, dtype=np.float64)

    @property
    def scenario_count(self) -> int:
        return len(self.scenario_ca_um2)

    def add(self, failure: np.ndarray, margin_um: np.ndarray | None = None) -> None:
        """Accumulate one scenario's binary hotspot (and optionally its margin map)."""
        failure = np.asarray(failure)
        if failure.shape != self.shape_yx:
            raise ValueError(f"failure shape {failure.shape} != {self.shape_yx}")
        self.failure_count += failure.astype(np.uint32)
        self.scenario_ca_um2.append(critical_area_um2(failure.astype(bool), self.grid))
        if margin_um is not None:
            margin = np.asarray(margin_um, dtype=np.float64)
            if margin.shape != self.shape_yx:
                raise ValueError(f"margin shape {margin.shape} != {self.shape_yx}")
            self.margin_sum += np.where(np.isfinite(margin), margin, 0.0)
            self.margin_samples += 1

    def probability(self) -> np.ndarray:
        """Failure fraction per pixel as float32 in [0, 1]."""
        if self.scenario_count == 0:
            raise ValueError("시나리오가 하나도 누적되지 않았습니다")
        return (self.failure_count.astype(np.float32) / np.float32(self.scenario_count))

    def mean_margin_um(self) -> np.ndarray:
        """Mean contact margin over the scenarios that provided one (NaN when none)."""
        if self.margin_samples == 0:
            return np.full(self.shape_yx, np.nan, dtype=np.float32)
        return (self.margin_sum / self.margin_samples).astype(np.float32)


@dataclass(frozen=True)
class ProbabilityCA:
    ca_probability_mean_um2: float  # spatial integral of the probability map
    ca_scenario_mean_um2: float  # mean of the per-scenario CA values (identical by linearity)
    ca_p05_um2: float
    ca_p50_um2: float
    ca_p95_um2: float
    scenario_count: int


def probability_ca(accumulator: ProbabilityAccumulator) -> ProbabilityCA:
    """CA of the probability map plus percentiles of the per-scenario CA distribution."""
    samples = np.asarray(accumulator.scenario_ca_um2, dtype=np.float64)
    if samples.size == 0:
        raise ValueError("시나리오가 하나도 누적되지 않았습니다")
    probability_map = accumulator.probability()
    return ProbabilityCA(
        ca_probability_mean_um2=critical_area_um2(probability_map, accumulator.grid),
        ca_scenario_mean_um2=float(samples.mean()),
        ca_p05_um2=float(np.percentile(samples, 5)),
        ca_p50_um2=float(np.percentile(samples, 50)),
        ca_p95_um2=float(np.percentile(samples, 95)),
        scenario_count=int(samples.size),
    )
