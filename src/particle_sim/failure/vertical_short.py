"""Vertical (and diagonal) short between nets on different layers."""

from __future__ import annotations

from particle_sim.failure.lateral_short import PairResult, build_pair_result
from particle_sim.geometry.contact import ContactResult
from particle_sim.schema import ElectrodePair, ProductManifest

# Particle conductivity, insulator continuity and step coverage are TODO_PHYSICS,
# so a geometric 3-D contact is never reported as a confirmed short.
VERTICAL_CLASSIFICATION = "potential_vertical_short"


def vertical_short(
    manifest: ProductManifest, pair: ElectrodePair, contacts: dict[str, ContactResult]
) -> PairResult:
    """Particle touches two different nets that sit at different heights."""
    return build_pair_result(
        manifest,
        pair,
        contacts,
        expected_mode="vertical_short",
        same_layer=False,
        classification=VERTICAL_CLASSIFICATION,
    )
