# phenoshot 0.1.0

## Breaking changes

* Object outlines are now drawn along the pixels that are actually measured, so
  background gaps enclosed between leaves appear as holes instead of being
  visually swallowed by the outline. The reported area is unchanged — it already
  excluded those gaps — but the annotated images look different. Set
  `outline_mode = "outer"` to restore the previous drawing.

## New features

* `weed_filter = "width"` removes thin non-target structures such as grass weed
  leaves using a distance-transform width filter. `leaf_min_width_cm` sets the
  threshold in centimeters so it stays valid across image resolutions, and
  `debug_weed_mask = TRUE` writes `<stem>_weed.png` showing exactly which pixels
  were discarded.
* `scale_marker = "aruco"` derives the pixel-to-centimeter scale from an ArUco
  marker of known size detected in the original image, removing the assumption
  that camera distance is identical for every capture. Guarded by `marker_id`,
  `marker_min_px` and `marker_max_dev`, which reject false positives and fall
  back to `image_real_cm` with a warning.
* `camera_height_cm` and `canopy_height_cm` correct the magnification caused by
  leaves sitting above the plane on which `image_real_cm` was measured.
  `canopy_height_cm` accepts `c(min, max)` to also report an uncertainty range.
* `object_order = "position"` numbers objects in reading order rather than by
  area, keeping Object ID tied to a physical location across dates.
* `params_json = TRUE` (default) writes `phenoshot_params_<timestamp>.json`
  recording every parameter, package versions, and the scale used per image.
* `k_close_iter` exposes the morphological closing iterations, previously fixed
  at 2. Set `0L` for the strictest area estimate.
* `outline_min_hole_cm2` controls the smallest gap that gets outlined.

## Output changes

* New columns: `Object Gap Area (cm2)`, `Object Area Min (cm2)`,
  `Object Area Max (cm2)`, `Canopy Height (cm)`, `Scale Source`.

## Notes

* `alpha_threshold` defaults to `10L` for backward compatibility, but `128L` is
  recommended: the lower value counts most of the semi-transparent halo that
  background removal leaves around each leaf edge, which inflates area for
  objects with long perimeters. Whichever value you choose, keep it constant
  across an experiment.
