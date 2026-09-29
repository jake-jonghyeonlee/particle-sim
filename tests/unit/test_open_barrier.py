"""Explicit-open barrier checks against both local and full-image stamp oracles."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest
from scipy.ndimage import distance_transform_edt

from particle_sim.failure import open_failure as of
from particle_sim.failure.open_barrier import (
    Ineligible,
    PreparedBarrier,
    _distance_to_sites,
    prepare_barrier,
)
from particle_sim.io.masks import GridSpec


def _case(electrode, a, b, sy=1.0, sx=1.0, *, regions=False):
    height, width = electrode.shape
    ta = np.zeros((height, width), bool)
    tb = np.zeros((height, width), bool)
    ta[a] = True
    tb[b] = True
    path = SimpleNamespace(path_id="P", electrode="E", terminal_a_yx_px=a,
                           terminal_b_yx_px=b, fallback=None)
    masks = SimpleNamespace(
        electrodes={"E": electrode},
        terminal_regions={"P": (ta, tb)} if regions else {},
        grid=GridSpec(shape_yx=electrode.shape, pixel_size_y_um=sy, pixel_size_x_um=sx),
    )
    manifest = SimpleNamespace(product_id="SYN", contact_margin_um=0.0)
    return manifest, masks, path, ta, tb


def _oracle(electrode, ta, tb, grid, diameter, margin):
    candidates = of.open_candidates(electrode, grid, diameter / 2, margin)
    failure = np.zeros_like(electrode)
    if of.baseline_status(electrode, ta, tb) is not None:
        return failure, 0
    stamp = of.particle_footprint(grid, diameter / 2, margin)
    hy, hx = np.array(stamp.shape) // 2
    for y, x in np.argwhere(candidates):
        top, bottom = max(0, y - hy), min(electrode.shape[0], y + hy + 1)
        left, right = max(0, x - hx), min(electrode.shape[1], x + hx + 1)
        blocked = np.zeros_like(electrode)
        blocked[top:bottom, left:right] = stamp[
            top - (y - hy):stamp.shape[0] - (y + hy + 1 - bottom),
            left - (x - hx):stamp.shape[1] - (x + hx + 1 - right),
        ]
        failure[y, x] = of.causes_open_regions(electrode, blocked, ta, tb)
    return failure, int(candidates.sum())


def _compare(case, diameter, margin, monkeypatch, *, full=True):
    manifest, masks, path, ta, tb = case
    manifest.contact_margin_um = margin
    monkeypatch.setenv("PARTICLE_SIM_OPEN_BACKEND", "local")
    local = of.simulate_open(manifest, masks, path, diameter)
    monkeypatch.setenv("PARTICLE_SIM_OPEN_BACKEND", "auto")
    fast = of.simulate_open(manifest, masks, path, diameter)
    assert fast.path_id == local.path_id
    assert fast.electrode == local.electrode
    assert fast.classification == local.classification
    assert fast.terminal_a_yx == local.terminal_a_yx
    assert fast.terminal_b_yx == local.terminal_b_yx
    assert fast.terminal_kind == local.terminal_kind
    assert fast.candidates == local.candidates
    assert fast.failure.dtype == local.failure.dtype == bool
    assert np.array_equal(fast.failure, local.failure)
    if full:
        expected, count = _oracle(masks.electrodes["E"], ta, tb, masks.grid,
                                  diameter, margin)
        assert fast.candidates == count
        assert np.array_equal(fast.failure, expected)


def test_certificate_and_fallback_reasons():
    e = np.zeros((9, 11), bool)
    e[:, 5] = True
    a = np.zeros_like(e)
    b = np.zeros_like(e)
    a[0, 5] = b[-1, 5] = True
    assert isinstance(prepare_barrier(e, a, b, 1, 1), PreparedBarrier)
    a[1, 5] = True
    assert prepare_barrier(e, a, b, 1, 1) == Ineligible("effective_terminal_not_singleton")
    a[1, 5] = False
    assert prepare_barrier(e, a, b, 0, 1) == Ineligible("unsafe_pixel_arithmetic")
    assert prepare_barrier(e, a, b, 1e-200, 1) == Ineligible("unsafe_pixel_arithmetic")
    a[0, 5] = False
    a[1, 5] = True
    assert prepare_barrier(e, a, b, 1, 1) == Ineligible(
        "terminals_not_opposite_noncorner_sides"
    )


@pytest.mark.parametrize("spacing", [1e110, 1e140])
def test_unsafe_edt_spacing_falls_back_to_local(spacing, monkeypatch):
    rng = np.random.default_rng(8317)
    electrode = rng.random((25, 31)) < .22
    electrode[:, 15] = True
    case = _case(electrode, (0, 15), (24, 15), spacing, spacing)
    assert prepare_barrier(electrode, case[3], case[4], spacing, spacing) == Ineligible(
        "unsafe_pixel_arithmetic"
    )
    for scale in (2.0, 3.1, 5.3):
        _compare(case, scale * spacing, 0, monkeypatch, full=False)


def test_adjacent_side_terminals_select_local(monkeypatch):
    electrode = np.zeros((9, 11), bool)
    electrode[:, 5] = True
    electrode[4, :6] = True
    case = _case(electrode, (0, 5), (4, 0))
    assert prepare_barrier(electrode, case[3], case[4], 1, 1) == Ineligible(
        "terminals_not_opposite_noncorner_sides"
    )
    _compare(case, 2.0, 0, monkeypatch)


@pytest.mark.parametrize("sy,sx", [(1, 1), (.1, .1), (.1, .2), (.2, .1), (.07, .11)])
def test_eligible_shapes_and_spacings(sy, sx, monkeypatch):
    shapes = []
    wire = np.zeros((11, 13), bool)
    wire[:, 6] = True
    shapes.append(wire)
    branch = wire.copy()
    branch[5, 6:12] = True
    shapes.append(branch)
    ring = wire.copy()
    ring[3:8, 4] = True
    ring[3:8, 8] = True
    ring[3, 4:9] = True
    ring[7, 4:9] = True
    shapes.append(ring)
    islands = branch.copy()
    islands[2:4, 10:12] = True
    shapes.append(islands)
    for electrode in shapes:
        case = _case(electrode, (0, 6), (10, 6), sy, sx)
        assert isinstance(prepare_barrier(electrode, case[3], case[4], sy, sx),
                          PreparedBarrier)
        for diameter, margin in [(0.4, 0), (2.0, 0), (3.0, .07), (6.0, .25)]:
            _compare(case, diameter, margin, monkeypatch)
    horizontal = np.zeros((13, 11), bool)
    horizontal[6, :] = True
    horizontal[3:7, 5] = True
    case = _case(horizontal, (6, 0), (6, 10), sy, sx)
    assert isinstance(prepare_barrier(horizontal, case[3], case[4], sy, sx),
                      PreparedBarrier)
    _compare(case, 2.0, .07, monkeypatch)
    large = np.zeros((32, 40), bool)
    large[:, 19] = True
    large[15, 19:39] = True
    _compare(_case(large, (0, 19), (31, 19), sy, sx), 3.0, .25, monkeypatch)


def test_fallback_regions_and_baseline(monkeypatch):
    e = np.zeros((9, 11), bool)
    e[:, 4:7] = True
    cases = [((0, slice(4, 7)), (8, slice(4, 7))),
             ((0, 5), (4, 5)), ((0, 5), (0, 6)), ((0, 5), (4, 4))]
    for a, b in cases:
        case = _case(e, a, b, regions=True)
        assert isinstance(prepare_barrier(e, case[3], case[4], 1, 1), Ineligible)
        _compare(case, 2.0, 0, monkeypatch)
    lost = _case(e, (0, 1), (8, 5))
    _compare(lost, 2.0, 0, monkeypatch)
    separate = e.copy()
    separate[4] = False
    _compare(_case(separate, (0, 5), (8, 5)), 2.0, 0, monkeypatch)
    corner = np.zeros((9, 11), bool)
    corner[:, 0] = True
    corner_case = _case(corner, (0, 0), (8, 0))
    assert isinstance(prepare_barrier(corner, corner_case[3], corner_case[4], 1, 1),
                      Ineligible)
    _compare(corner_case, 2.0, 0, monkeypatch)
    overlap_case = _case(e, (0, slice(4, 7)), (0, slice(5, 7)), regions=True)
    assert isinstance(prepare_barrier(e, overlap_case[3], overlap_case[4], 1, 1),
                      Ineligible)
    _compare(overlap_case, 2.0, 0, monkeypatch)


def test_effective_singleton_with_raw_off_electrode_pixels(monkeypatch):
    e = np.zeros((9, 11), bool)
    e[:, 4:7] = True
    case = _case(e, (0, slice(3, 5)), (8, slice(6, 8)), regions=True)
    prepared = prepare_barrier(e, case[3], case[4], 1, 1)
    assert isinstance(prepared, PreparedBarrier)
    _compare(case, 2, .07, monkeypatch)
    result = of.simulate_open(*case[:3], 2)
    assert result.terminal_a_yx == (0, 3)
    assert result.terminal_b_yx == (8, 6)
    assert result.terminal_kind == "region_mask"


def test_seeded_random_masks_and_reverse_order(monkeypatch):
    rng = np.random.default_rng(28831)
    for i in range(200):
        e = rng.random((7, 9)) < .22
        e[:, 4] = True
        case = _case(e, (0, 4), (6, 4), .1 if i % 2 else 1,
                     .2 if i % 2 else 1)
        for diameter in ((2.0, 1.0) if i % 2 else (1.0, 2.0)):
            _compare(case, diameter, 0, monkeypatch, full=i < 40)


def test_footprint_box_repairs_provisional_true(monkeypatch):
    e = np.ones((90, 100), bool)
    case = _case(e, (0, 50), (89, 50), .1, .1)
    prepared = prepare_barrier(e, case[3], case[4], .1, .1)
    assert isinstance(prepared, PreparedBarrier)
    assert np.hypot(0, 43 * .1) == 4.3
    assert of.particle_footprint(case[1].grid, 4.3, 0).shape[1] // 2 == 42
    assert prepared.critical[43, 50] <= 4.3
    _compare(case, 8.6, 0, monkeypatch)
    assert not of.simulate_open(*case[:3], 8.6).failure[43, 50]


def test_genuine_edt_seed_tie_repairs_provisional_false(monkeypatch):
    # Two right-side sites have squared offsets (9,20) and (15,16).
    # SciPy chooses the latter; np.hypot puts it one ulp above the radius.
    e = np.ones((61, 61), bool)
    e[40, :10] = False
    e[46, 47:] = False
    e[40, 51:] = False
    case = _case(e, (0, 30), (60, 30), .1, .1)
    reach = float(np.hypot(.9, 2.0))
    diameter = 2 * reach
    assert diameter / 2 == reach
    prepared = prepare_barrier(e, case[3], case[4], .1, .1)
    assert isinstance(prepared, PreparedBarrier)
    assert prepared.critical[30, 30] > reach
    assert prepared.ambiguous(reach)[30, 30]
    _compare(case, diameter, 0, monkeypatch)
    result = of.simulate_open(*case[:3], diameter)
    assert result.failure[30, 30]


def test_selected_site_uses_reference_hypot_not_raw_edt():
    sites = np.zeros((101, 101), bool)
    sites[50, 50] = True
    expected = np.hypot(50 * .1, 44 * .1)
    raw = distance_transform_edt(~sites, sampling=(.1, .1))[0, 6]
    assert raw != expected  # This is the one-ulp SciPy/NumPy regression premise.
    assert _distance_to_sites(sites, .1, .1)[0, 6] == expected


def test_roundoff_envelope_keeps_image_extent_term():
    reach = 2.5
    critical = np.full((600, 1), reach + 8 * np.finfo(float).eps * 60)
    prepared = PreparedBarrier(critical, .1, .1)
    assert prepared.ambiguous(reach)[0, 0]
    assert not prepared.ambiguous(reach + .01)[0, 0]


def test_lattice_radii_and_adjacent_floats(monkeypatch):
    e = np.zeros((13, 17), bool)
    e[:, 8] = True
    e[5, 8:15] = True
    case = _case(e, (0, 8), (12, 8), .1, .1)
    for dy, dx in [(3, 4), (7, 24), (15, 20), (9, 20), (15, 16)]:
        radius = float(np.hypot(dy * .1, dx * .1))
        for reach in (np.nextafter(radius, -np.inf), radius,
                      np.nextafter(radius, np.inf)):
            diameter = 2 * reach
            assert diameter / 2 == reach
            _compare(case, diameter, 0, monkeypatch)


def test_candidate_intersection_is_enforced(monkeypatch):
    e = np.zeros((9, 11), bool)
    e[:, 5] = True
    case = _case(e, (0, 5), (8, 5))
    fake_critical = np.zeros(e.shape, dtype=np.float64)
    fake_critical.setflags(write=False)
    monkeypatch.setenv("PARTICLE_SIM_OPEN_BACKEND", "auto")
    monkeypatch.setattr(of, "prepare_barrier", lambda *_: PreparedBarrier(fake_critical, 1, 1))
    monkeypatch.setattr(of, "open_candidates", lambda *_: e.copy())
    result = of.simulate_open(*case[:3], 2)
    assert np.array_equal(result.failure, e)
    assert result.candidates == int(e.sum())


def test_ineligible_certificate_selects_local_backend(monkeypatch):
    e = np.zeros((9, 11), bool)
    e[:, 5] = True
    case = _case(e, (0, 5), (8, 5))
    monkeypatch.setenv("PARTICLE_SIM_OPEN_BACKEND", "auto")
    monkeypatch.setattr(of, "prepare_barrier", lambda *_: Ineligible("probe"))
    original = of._evaluate_local_centres
    calls = []

    def local(*args):
        calls.append(int(args[3].sum()))
        return original(*args)

    monkeypatch.setattr(of, "_evaluate_local_centres", local)
    result = of.simulate_open(*case[:3], 2)
    assert calls == [result.candidates]


def test_zero_reach_invalid_backend_and_force_local(monkeypatch):
    e = np.zeros((9, 11), bool)
    e[:, 5] = True
    case = _case(e, (0, 5), (8, 5))
    _compare(case, 0, 0, monkeypatch)
    _compare(case, 1, -1, monkeypatch)
    monkeypatch.setenv("PARTICLE_SIM_OPEN_BACKEND", "nonsense")
    with pytest.raises(ValueError, match="PARTICLE_SIM_OPEN_BACKEND"):
        of.simulate_open(*case[:3], 2)
    monkeypatch.setenv("PARTICLE_SIM_OPEN_BACKEND", "local")
    case[0].contact_margin_um = 0
    monkeypatch.setattr(of, "prepare_barrier", lambda *_: pytest.fail("prepared local"))
    assert of.simulate_open(*case[:3], 2).candidates > 0


def test_oversized_stamp_uses_local_validation_path(monkeypatch):
    e = np.zeros((9, 11), bool)
    e[:, 5] = True
    case = _case(e, (0, 5), (8, 5))
    monkeypatch.setenv("PARTICLE_SIM_OPEN_BACKEND", "auto")
    monkeypatch.setattr(of, "prepare_barrier", lambda *_: pytest.fail("prepared oversized"))
    assert of.simulate_open(*case[:3], 300).failure.all()
