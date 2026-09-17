"""uncertainty sampler 단위 테스트."""

from __future__ import annotations

import numpy as np
import pytest

from particle_sim.errors import InputContractError
from particle_sim.schema import (
    NormalSpec,
    ScenarioConfig,
    UniformSpec,
)
from particle_sim.uncertainty.sampler import (
    build_scenarios,
    make_sobol_unit_samples,
    variable_names,
)
from particle_sim.uncertainty.transforms import (
    clamp,
    normal_from_unit,
    uniform_from_unit,
)


@pytest.mark.parametrize(
    ("u", "expected"),
    [(0.5, 10.0), (0.0, 0.0), (1.0, 20.0), (0.25, 5.0), (0.75, 15.0)],
)
def test_uniform_from_unit(u: float, expected: float) -> None:
    spec = UniformSpec(low=0.0, high=20.0)
    result = uniform_from_unit(u, spec)
    assert result.shape == ()
    assert result.dtype == np.float64
    assert result == pytest.approx(expected)


def test_uniform_from_unit_constant() -> None:
    spec = UniformSpec(low=3.0, high=3.0)
    result = uniform_from_unit(0.7, spec)
    assert result == pytest.approx(3.0)


def test_uniform_from_unit_array() -> None:
    spec = UniformSpec(low=1.0, high=3.0)
    result = uniform_from_unit(np.array([0.0, 0.5, 1.0]), spec)
    assert result.shape == (3,)
    assert result.dtype == np.float64
    np.testing.assert_allclose(result, [1.0, 2.0, 3.0])


@pytest.mark.parametrize("u", [0.0, 0.5, 1.0])
def test_normal_from_unit_sigma_zero(u: float) -> None:
    spec = NormalSpec(mean=10.0, sigma=0.0)
    result = normal_from_unit(u, spec)
    assert result == pytest.approx(10.0)


def test_normal_from_unit_symmetry() -> None:
    spec = NormalSpec(mean=0.0, sigma=1.0)
    low = normal_from_unit(0.2, spec)
    high = normal_from_unit(0.8, spec)
    assert low + high == pytest.approx(0.0, abs=1e-9)
    assert low < 0.0 < high


def test_normal_from_unit_array() -> None:
    spec = NormalSpec(mean=0.0, sigma=1.0)
    result = normal_from_unit(np.array([0.5, 0.5]), spec)
    assert result.shape == (2,)
    np.testing.assert_allclose(result, [0.0, 0.0])


@pytest.mark.parametrize("u", [-0.1, 1.1])
def test_normal_from_unit_out_of_range(u: float) -> None:
    with pytest.raises(InputContractError):
        normal_from_unit(u, NormalSpec(mean=0.0, sigma=1.0))


@pytest.mark.parametrize("u", [-0.1, 1.1])
def test_uniform_from_unit_out_of_range(u: float) -> None:
    with pytest.raises(InputContractError):
        uniform_from_unit(u, UniformSpec(low=0.0, high=1.0))


def test_clamp() -> None:
    values = np.array([-1.0, 0.5, 2.0])
    result = clamp(values, 0.0, 1.0)
    np.testing.assert_allclose(result, [0.0, 0.5, 1.0])


def test_clamp_no_bounds() -> None:
    values = np.array([-1.0, 0.5, 2.0])
    result = clamp(values, None, None)
    np.testing.assert_allclose(result, values)


def test_clamp_scalar() -> None:
    result = clamp(5.0, 0.0, 3.0)
    assert result == pytest.approx(3.0)


def test_variable_names_order() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"b": NormalSpec(mean=0.0, sigma=1.0), "a": NormalSpec(mean=0.0, sigma=1.0)},
        overlay_um={"l2": NormalSpec(mean=0.0, sigma=1.0), "l1": NormalSpec(mean=0.0, sigma=1.0)},
        cd_bias_um={"c": NormalSpec(mean=0.0, sigma=1.0)},
        particle_diameter_delta_um=NormalSpec(mean=0.0, sigma=1.0),
        center_offset_um=NormalSpec(mean=0.0, sigma=1.0),
        embedding_depth_um=UniformSpec(low=0.0, high=1.0),
        planarization_factor={"p": UniformSpec(low=0.0, high=1.0)},
        contact_margin_um=UniformSpec(low=0.0, high=1.0),
    )
    expected = [
        "thickness_delta_um.a",
        "thickness_delta_um.b",
        "overlay_um.l1.x",
        "overlay_um.l1.y",
        "overlay_um.l2.x",
        "overlay_um.l2.y",
        "cd_bias_um.c",
        "particle_diameter_delta_um",
        "center_offset_um.x",
        "center_offset_um.y",
        "embedding_depth_um",
        "planarization_factor.p",
        "contact_margin_um",
    ]
    assert variable_names(config) == expected


def test_variable_names_partial() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"a": NormalSpec(mean=0.0, sigma=1.0)},
    )
    assert variable_names(config) == ["thickness_delta_um.a"]


def test_make_sobol_shape_and_range() -> None:
    result = make_sobol_unit_samples(3, 6, 42)
    assert result.shape == (64, 3)
    assert np.all(result >= 0.0)
    assert np.all(result < 1.0)


def test_make_sobol_reproducible() -> None:
    a = make_sobol_unit_samples(3, 6, 42)
    b = make_sobol_unit_samples(3, 6, 42)
    np.testing.assert_array_equal(a, b)


def test_make_sobol_different_seed() -> None:
    a = make_sobol_unit_samples(3, 6, 42)
    b = make_sobol_unit_samples(3, 6, 43)
    assert not np.allclose(a, b)


@pytest.mark.parametrize("dims, exp", [(0, 6), (-1, 6), (3, -1)])
def test_make_sobol_invalid(dims: int, exp: int) -> None:
    with pytest.raises(InputContractError):
        make_sobol_unit_samples(dims, exp, 42)


def test_build_scenarios_degenerate() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"a": NormalSpec(mean=0.0, sigma=0.0)},
        overlay_um={"l": NormalSpec(mean=0.0, sigma=0.0)},
        cd_bias_um={"c": NormalSpec(mean=0.0, sigma=0.0)},
        particle_diameter_delta_um=NormalSpec(mean=0.0, sigma=0.0),
        center_offset_um=NormalSpec(mean=0.0, sigma=0.0),
        embedding_depth_um=UniformSpec(low=2.0, high=2.0),
        planarization_factor={"p": UniformSpec(low=0.5, high=0.5)},
        contact_margin_um=UniformSpec(low=1.0, high=1.0),
    )
    assert config.is_degenerate()
    scenarios = build_scenarios(config, 6, 7)
    assert len(scenarios) == 1
    s = scenarios[0]
    assert s.scenario_id == 0
    assert s.thickness_delta_um == {"a": 0.0}
    assert s.overlay_xy_um == {"l": (0.0, 0.0)}
    assert s.cd_bias_um == {"c": 0.0}
    assert s.particle_diameter_delta_um == 0.0
    assert s.center_offset_xy_um == (0.0, 0.0)
    assert s.embedding_depth_um == 2.0
    assert s.planarization_factor == {"p": 0.5}
    assert s.contact_margin_um == 1.0


def test_build_scenarios_none_optional() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"a": NormalSpec(mean=0.0, sigma=1.0)},
    )
    scenarios = build_scenarios(config, 6, 7)
    assert len(scenarios) == 64
    for s in scenarios:
        assert s.embedding_depth_um is None
        assert s.contact_margin_um is None
        assert s.center_offset_xy_um == (0.0, 0.0)
        assert s.particle_diameter_delta_um == 0.0


def test_build_scenarios_count_and_ids() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"a": NormalSpec(mean=0.0, sigma=1.0)},
    )
    scenarios = build_scenarios(config, 6, 7)
    assert len(scenarios) == 64
    assert [s.scenario_id for s in scenarios] == list(range(64))


def test_build_scenarios_reproducible() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"a": NormalSpec(mean=0.0, sigma=1.0)},
    )
    a = build_scenarios(config, 6, 7)
    b = build_scenarios(config, 6, 7)
    for sa, sb in zip(a, b, strict=True):
        assert sa == sb


def test_build_scenarios_mean_near_zero() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"a": NormalSpec(mean=0.0, sigma=1.0)},
    )
    scenarios = build_scenarios(config, 6, 7)
    mean = np.mean([s.thickness_delta_um["a"] for s in scenarios])
    assert abs(mean) < 0.2


def test_build_scenarios_keys_match_config() -> None:
    config = ScenarioConfig(
        schema_version="1.0",
        assumption_set="test",
        evidence="assumed",
        thickness_delta_um={"a": NormalSpec(mean=0.0, sigma=1.0), "b": NormalSpec(mean=0.0, sigma=1.0)},
        overlay_um={"l1": NormalSpec(mean=0.0, sigma=1.0)},
        cd_bias_um={"c": NormalSpec(mean=0.0, sigma=1.0)},
        planarization_factor={"p": UniformSpec(low=0.0, high=1.0)},
    )
    scenarios = build_scenarios(config, 6, 7)
    for s in scenarios:
        assert set(s.thickness_delta_um) == {"a", "b"}
        assert set(s.overlay_xy_um) == {"l1"}
        assert set(s.cd_bias_um) == {"c"}
        assert set(s.planarization_factor) == {"p"}
