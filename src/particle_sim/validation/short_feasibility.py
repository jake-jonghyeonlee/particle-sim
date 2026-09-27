"""Conservative structural-zero warnings for nominal short geometry."""

from particle_sim.io.stack import LayerZ
from particle_sim.schema import ProductManifest


def structural_zero_short_warnings(
    manifest: ProductManifest, stack_z: dict[str, LayerZ]
) -> list[str]:
    """Warn only for proven impossible declared shorts; do not alter failure units.

    At/below a flat landing plane, zero-margin unembedded spheres can only
    touch at one XY point. Validated same-layer masks are disjoint; different
    earlier layers are buried. Later layers can wrap particles, so conservatively
    omit those pairs unless the exclude policy prevents their contact entirely.
    """
    if manifest.contact_margin_um != 0 or manifest.embedding_depth_um != 0:
        return []
    electrodes = {e.electrode_id: e for e in manifest.electrodes}
    pairs = [p for p in manifest.electrode_pairs
             if p.failure_mode in ("lateral_short", "vertical_short")]
    warnings = []
    for step in manifest.particle_source_steps:
        if step.surface_model != "flat":
            continue
        landing_order = stack_z[step.lands_after_layer].process_order
        count = sum(
            all(stack_z[electrodes[e].layer_name].process_order <= landing_order
                for e in (pair.electrode_a, pair.electrode_b))
            or manifest.post_step_policy.mode == "exclude"
            for pair in pairs
        )
        if count:
            warnings.append(
                f"STRUCTURAL_ZERO_SHORT: step {step.step_id}: flat surface with zero "
                f"contact margin and embedding; {count} of {len(pairs)} short pairs "
                "cannot be contacted simultaneously in nominal geometry "
                "(at/below landing plane: point contact or buried electrode; "
                "later layers counted only under exclude policy)"
            )
    return warnings
