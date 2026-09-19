"""결정적 Sobol 기반 공정 산포 시나리오 생성."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.stats import qmc

from particle_sim.errors import InputContractError
from particle_sim.schema import NormalSpec, OverlayXYSpec, ScenarioConfig
from particle_sim.uncertainty.transforms import normal_from_unit, uniform_from_unit


@dataclass(frozen=True)
class Scenario:
    """하나의 공정 산포 시나리오."""

    scenario_id: int
    thickness_delta_um: dict[str, float]
    overlay_xy_um: dict[str, tuple[float, float]]
    cd_bias_um: dict[str, float]
    particle_diameter_delta_um: float
    center_offset_xy_um: tuple[float, float]
    embedding_depth_um: float | None
    planarization_factor: dict[str, float]
    contact_margin_um: float | None
    scramble_seed: int = 0
    sample_index: int = 0


def variable_names(config: ScenarioConfig) -> list[str]:
    """설정에서 존재하는 변수를 결정적 순서로 나열한다."""
    names: list[str] = []
    for key in sorted(config.thickness_delta_um):
        names.append(f"thickness_delta_um.{key}")
    for key in sorted(config.overlay_um):
        names.append(f"overlay_um.{key}.x")
        names.append(f"overlay_um.{key}.y")
    for key in sorted(config.cd_bias_um):
        names.append(f"cd_bias_um.{key}")
    if config.particle_diameter_delta_um is not None:
        names.append("particle_diameter_delta_um")
    if config.center_offset_um is not None:
        names.append("center_offset_um.x")
        names.append("center_offset_um.y")
    if config.embedding_depth_um is not None:
        names.append("embedding_depth_um")
    for key in sorted(config.planarization_factor):
        names.append(f"planarization_factor.{key}")
    if config.contact_margin_um is not None:
        names.append("contact_margin_um")
    return names


def _overlay_axis_spec(config: ScenarioConfig, layer: str, axis: str) -> NormalSpec:
    spec = config.overlay_um[layer]
    if isinstance(spec, OverlayXYSpec):
        return spec.x if axis == "x" else spec.y
    return spec


def make_sobol_unit_samples(dimensions: int, exponent_m: int, seed: int) -> np.ndarray:
    """scrambled Sobol 단위 표본을 생성한다."""
    if dimensions <= 0 or exponent_m < 0:
        raise InputContractError("Sobol 차원과 지수는 유효해야 합니다", field="sobol")
    sobol = qmc.Sobol(d=dimensions, scramble=True, seed=seed)
    return sobol.random_base2(m=exponent_m)


def _fixed_values(config: ScenarioConfig) -> dict[str, float]:
    """산포 0 시나리오의 고정값을 변수 이름으로 매핑한다."""
    values: dict[str, float] = {}
    for key, spec in config.thickness_delta_um.items():
        values[f"thickness_delta_um.{key}"] = spec.mean
    for key in config.overlay_um:
        values[f"overlay_um.{key}.x"] = _overlay_axis_spec(config, key, "x").mean
        values[f"overlay_um.{key}.y"] = _overlay_axis_spec(config, key, "y").mean
    for key, spec in config.cd_bias_um.items():
        values[f"cd_bias_um.{key}"] = spec.mean
    if config.particle_diameter_delta_um is not None:
        values["particle_diameter_delta_um"] = config.particle_diameter_delta_um.mean
    if config.center_offset_um is not None:
        values["center_offset_um.x"] = config.center_offset_um.mean
        values["center_offset_um.y"] = config.center_offset_um.mean
    if config.embedding_depth_um is not None:
        values["embedding_depth_um"] = config.embedding_depth_um.low
    for key, spec in config.planarization_factor.items():
        values[f"planarization_factor.{key}"] = spec.low
    if config.contact_margin_um is not None:
        values["contact_margin_um"] = config.contact_margin_um.low
    return values


def _build_one(
    config: ScenarioConfig,
    scenario_id: int,
    values: dict[str, float],
    *,
    seed: int = 0,
    sample_index: int | None = None,
) -> Scenario:
    thickness = {
        key: values[f"thickness_delta_um.{key}"] for key in config.thickness_delta_um
    }
    overlay = {
        key: (values[f"overlay_um.{key}.x"], values[f"overlay_um.{key}.y"])
        for key in config.overlay_um
    }
    cd_bias = {key: values[f"cd_bias_um.{key}"] for key in config.cd_bias_um}
    particle_diameter = values.get("particle_diameter_delta_um", 0.0)
    center_offset = (
        (values["center_offset_um.x"], values["center_offset_um.y"])
        if config.center_offset_um is not None
        else (0.0, 0.0)
    )
    embedding_depth = values.get("embedding_depth_um")
    planarization = {
        key: values[f"planarization_factor.{key}"] for key in config.planarization_factor
    }
    contact_margin = values.get("contact_margin_um")
    return Scenario(
        scenario_id=scenario_id,
        thickness_delta_um=thickness,
        overlay_xy_um=overlay,
        cd_bias_um=cd_bias,
        particle_diameter_delta_um=particle_diameter,
        center_offset_xy_um=center_offset,
        embedding_depth_um=embedding_depth,
        planarization_factor=planarization,
        contact_margin_um=contact_margin,
        scramble_seed=seed,
        sample_index=scenario_id if sample_index is None else sample_index,
    )


def build_scenarios(config: ScenarioConfig, exponent_m: int, seed: int) -> list[Scenario]:
    """설정과 Sobol 매개변수로 시나리오 목록을 생성한다.

    모든 산포가 0이면 고정값을 가진 시나리오 1개만 반환한다.
    """
    if config.is_degenerate():
        return [_build_one(config, 0, _fixed_values(config), seed=seed, sample_index=0)]

    names = variable_names(config)
    unit = make_sobol_unit_samples(len(names), exponent_m, seed)
    scenarios: list[Scenario] = []
    for row, unit_row in enumerate(unit):
        values: dict[str, float] = {}
        for name, u in zip(names, unit_row, strict=True):
            if name.startswith("thickness_delta_um."):
                key = name.split(".", 1)[1]
                values[name] = float(
                    normal_from_unit(u, config.thickness_delta_um[key])
                )
            elif name.startswith("overlay_um."):
                key, axis = name.split(".", 2)[1:3]
                values[name] = float(
                    normal_from_unit(u, _overlay_axis_spec(config, key, axis))
                )
            elif name.startswith("cd_bias_um."):
                key = name.split(".", 1)[1]
                values[name] = float(normal_from_unit(u, config.cd_bias_um[key]))
            elif name == "particle_diameter_delta_um":
                assert config.particle_diameter_delta_um is not None
                values[name] = float(
                    normal_from_unit(u, config.particle_diameter_delta_um)
                )
            elif name.startswith("center_offset_um."):
                assert config.center_offset_um is not None
                values[name] = float(normal_from_unit(u, config.center_offset_um))
            elif name == "embedding_depth_um":
                assert config.embedding_depth_um is not None
                values[name] = float(uniform_from_unit(u, config.embedding_depth_um))
            elif name.startswith("planarization_factor."):
                key = name.split(".", 1)[1]
                values[name] = float(
                    uniform_from_unit(u, config.planarization_factor[key])
                )
            elif name == "contact_margin_um":
                assert config.contact_margin_um is not None
                values[name] = float(uniform_from_unit(u, config.contact_margin_um))
        scenarios.append(_build_one(config, row, values, seed=seed, sample_index=row))
    return scenarios
