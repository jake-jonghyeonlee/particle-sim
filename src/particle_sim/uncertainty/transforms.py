"""단위 균등 표본을 공정 산포 분포 값으로 변환하는 함수."""

from __future__ import annotations

import numpy as np
from scipy.stats import norm

from particle_sim.errors import InputContractError
from particle_sim.schema import NormalSpec, UniformSpec


def _check_unit(u: np.ndarray) -> None:
    if np.any(u < 0.0) or np.any(u > 1.0):
        raise InputContractError("단위 표본은 [0, 1] 범위에 있어야 합니다", field="unit_sample")


def normal_from_unit(u: np.ndarray | float, spec: NormalSpec) -> np.ndarray:
    """단위 균등 표본을 정규 분포 값으로 역 CDF 변환한다.

    sigma == 0이면 ppf를 호출하지 않고 평균 상수값을 반환한다.
    """
    arr = np.asarray(u, dtype=np.float64)
    _check_unit(arr)
    if spec.sigma == 0.0:
        return np.full(arr.shape, spec.mean, dtype=np.float64)
    return norm.ppf(arr, loc=spec.mean, scale=spec.sigma)


def uniform_from_unit(u: np.ndarray | float, spec: UniformSpec) -> np.ndarray:
    """단위 균등 표본을 구간 균등 분포 값으로 변환한다.

    low == high이면 고정값을 반환한다.
    """
    arr = np.asarray(u, dtype=np.float64)
    _check_unit(arr)
    if spec.low == spec.high:
        return np.full(arr.shape, spec.low, dtype=np.float64)
    return spec.low + (spec.high - spec.low) * arr


def clamp(values: np.ndarray | float, low: float | None, high: float | None) -> np.ndarray:
    """값을 [low, high] 구간으로 제한한다. None 경계는 무시한다."""
    arr = np.asarray(values, dtype=np.float64)
    out = arr
    if low is not None:
        out = np.maximum(out, low)
    if high is not None:
        out = np.minimum(out, high)
    return out
