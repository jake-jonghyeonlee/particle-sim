# arch: mask and stack I/O

`load_product_masks` reads 2-D binary PNG/TIFF files without resizing or thresholding. All layer,
electrode, and terminal masks must share one shape and match the configured ROI within half a pixel.

Mask relations:

- same-layer electrode masks must not overlap;
- every electrode mask must be a nonempty subset of its layer mask;
- the electrode union may be smaller than the layer mask because TR or unevaluated geometry can remain;
- residual pixels per layer are returned as `MaskSet.layer_unassigned_pixels`;
- terminal masks must be nonempty, disjoint subsets of their open-path electrode.

`MaskSet` stores `layers`, `electrodes`, terminal region pairs, input hashes and `GridSpec`. `GridSpec`
keeps `(y,x)` shape, separate x/y pixel sizes, pixel area and ROI area.

`build_stack_z` sorts by process order and accumulates nominal `z_bottom_um`/`z_top_um`. Explicit
`z_bottom_um` values are cross-checks. Landing, prior and later layer helpers use this ordered stack.
