"""Size monotonicity of open judgement is FALSE for simulate_auto_open.

The premise "open at diameter d1 implies open at every d2 > d1" does not hold
when contact-derived terminals are involved.  The CNO rule (Rule 3 in
auto_open.py) skips a contact terminal that is fully covered by the particle;
a larger particle can cover a contact that a smaller one left alive, removing
it from the open comparison and potentially collapsing a split into a single
group.

This test locks the counterexample so that any future optimisation that assumes
size monotonicity for auto_open will be caught.

For simulate_open (explicit 2-terminal paths) the premise holds: terminal loss
always counts as open, so a larger disc that covers more can only add opens.
"""

from __future__ import annotations

import numpy as np

from particle_sim.failure.auto_open import simulate_auto_open
from particle_sim.failure.open_failure import simulate_open
from particle_sim.io.masks import GridSpec, MaskSet
from particle_sim.schema import OpenPath

# ── helpers ──────────────────────────────────────────────────────────


class _FakeHole:
    def __init__(
        self,
        contact_id: str,
        hole_index: int,
        lower_electrode: str,
        upper_electrode: str,
        ys: np.ndarray,
        xs: np.ndarray,
    ) -> None:
        self.contact_id = contact_id
        self.hole_index = hole_index
        self.lower_electrode = lower_electrode
        self.upper_electrode = upper_electrode
        self.ys = ys
        self.xs = xs


class _FakeTopology:
    def __init__(self, holes: list[_FakeHole]) -> None:
        self.holes = holes
        self.links: dict = {}
        self.lost_holes = 0


CANVAS = 14
SHAPE = (CANVAS, CANVAS)
GRID = GridSpec(shape_yx=SHAPE, pixel_size_y_um=1.0, pixel_size_x_um=1.0)


def _minimal_counterexample() -> MaskSet:
    """14×14 canvas, 2-pixel-wide wire, 4-px contact left, 1-px contact right.

    At centre (7, 4):
      d=2 → particle splits wire, both contacts alive in different components → open.
      d=3 → larger particle fully covers the left contact → CNO → not open.
    """
    electrode = np.zeros(SHAPE, dtype=bool)
    electrode[7:9, 2:12] = True
    h_left = _FakeHole(
        "C", 0, "X", "W",
        ys=np.array([7, 7, 8, 8]),
        xs=np.array([3, 4, 3, 4]),
    )
    h_right = _FakeHole(
        "C", 1, "X", "W",
        ys=np.array([7]),
        xs=np.array([10]),
    )
    return MaskSet(
        grid=GRID,
        layers={},
        electrodes={"W": electrode},
        input_hashes={},
        contact_topology=_FakeTopology([h_left, h_right]),
    )


# ── tests ────────────────────────────────────────────────────────────


class TestAutoOpenMonotonicityFails:
    """simulate_auto_open is NOT monotone in particle diameter."""

    def test_counterexample_exists(self) -> None:
        """There exists a centre that is open at d=2 but not open at d=3."""
        masks = _minimal_counterexample()
        r2 = simulate_auto_open(masks, "W", 2.0, 0.0)
        r3 = simulate_auto_open(masks, "W", 3.0, 0.0)
        violations = r2.failure & ~r3.failure
        assert violations.any(), (
            "Expected at least one centre open at d=2 but not at d=3"
        )

    def test_specific_centre(self) -> None:
        """Centre (7, 4) is the minimal violation point."""
        masks = _minimal_counterexample()
        r2 = simulate_auto_open(masks, "W", 2.0, 0.0)
        r3 = simulate_auto_open(masks, "W", 3.0, 0.0)
        assert r2.failure[7, 4], "centre (7,4) should be open at d=2"
        assert not r3.failure[7, 4], "centre (7,4) should NOT be open at d=3"


class TestSimulateOpenMonotonicityHolds:
    """simulate_open (explicit 2-terminal) IS monotone in particle diameter."""

    def test_no_violation_on_wire(self) -> None:
        """No centre is open at d1 but not at d2 > d1."""
        electrode = np.zeros(SHAPE, dtype=bool)
        electrode[7:9, 2:12] = True
        masks = MaskSet(
            grid=GRID,
            layers={},
            electrodes={"W": electrode},
            input_hashes={},
            contact_topology=None,
        )
        path = OpenPath(
            path_id="P",
            electrode="W",
            terminal_a_yx_px=(7, 3),
            terminal_b_yx_px=(7, 10),
        )

        class _Manifest:
            product_id = "test"
            contact_margin_um = 0.0

        manifest = _Manifest()
        diameters = [2.0, 3.0, 4.0, 6.0, 8.0]
        failures: dict[float, np.ndarray] = {}
        for d in diameters:
            result = simulate_open(manifest, masks, path, d, contact_margin_um=0.0)
            failures[d] = result.failure.copy()

        for i, d1 in enumerate(diameters):
            for d2 in diameters[i + 1 :]:
                violation = failures[d1] & ~failures[d2]
                assert not violation.any(), (
                    f"simulate_open monotonicity violated: "
                    f"d1={d1} d2={d2} violations={int(violation.sum())}"
                )
