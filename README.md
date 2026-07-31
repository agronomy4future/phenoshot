# phenoshot

**AI-based plant image analysis for morphological trait measurement**

`phenoshot` is an R package that measures plant morphological traits from
photographs. It combines AI background removal (Photoroom API) with
OpenCV-based object detection to report **area, perimeter, length, and width**
in real-world units (cm), and writes an annotated image so every measured
object can be visually verified.

It can also discard thin non-target structures such as grass weed leaves, and
correct the pixel-to-centimetre scale either from a fiducial marker or from the
camera and canopy heights.

Supported input formats: **JPG, PNG, WEBP, HEIC**.

---

## Installation

```r
# Install the package
if(!require(remotes)) install.packages("remotes")
if (!requireNamespace("phenoshot", quietly= TRUE)) {
  remotes::install_github("agronomy4future/phenoshot", force= TRUE)
}
library(remotes)
library(phenoshot)
```

Python dependencies (`numpy`, `pandas`, `requests`, `pillow`, `pillow-heif`,
`opencv-python-headless`) are installed automatically on first run through
`reticulate`. No manual Python setup is required.

> `scale_marker = "aruco"` additionally needs `cv2.aruco`, which ships only with
> the contrib build. Install it yourself if you use that option:
> `reticulate::py_install("opencv-contrib-python-headless")`. Every other
> feature works with the default install.

---

## API key

Background removal uses the [Photoroom Remove Background API](https://www.photoroom.com/api).
Get a key from the Photoroom dashboard and pass it as `photoroom_api_key`.

A key is **not** required when pre-processed `_nobg.png` files already exist —
see *Mode 1* below.

Keep the key out of your scripts. Put it in `~/.Renviron`:

```
PHOTOROOM_API_KEY=your_key_here
```

```r
phenoshot(..., photoroom_api_key = Sys.getenv("PHOTOROOM_API_KEY"))
```

> **Sandbox keys**: prefixing your key with `sandbox_` enables free testing, but
> results are watermarked. The watermark contaminates the alpha mask and
> corrupts area measurements, so sandbox keys must not be used for real
> analysis. `phenoshot()` raises a warning if it detects one.

---

## Usage

### Mode 1 — Pre-processed images (no API key)

If your input folder already contains `_nobg.png` files (background removed
elsewhere), `phenoshot` uses them directly and makes no API calls.

```r
phenoshot(
  input_folder  = "./Input",
  output_folder = "./Output",
  image_real_cm = c(20, 20)
)
```

### Mode 2 — Original photographs (Photoroom API)

```r
phenoshot(
  input_folder      = "./Input",
  output_folder     = "./Output",
  image_real_cm     = c(75, 75),
  photoroom_api_key = Sys.getenv("PHOTOROOM_API_KEY"),
  photoroom_size    = "full",
  distinct_colors   = TRUE,
  fill_opacity      = 0.25
)
```

### Mode 3 — Field canopy with a ground quadrat

```r
phenoshot(
  input_folder        = "./Input",
  output_folder       = "./Output",
  image_real_cm       = c(75, 75),   # quadrat inner size, at the ground
  photoroom_api_key   = Sys.getenv("PHOTOROOM_API_KEY"),

  camera_height_cm    = 100,         # lens to ground, measured
  canopy_height_cm    = c(10, 30),   # leaf layer above the ground

  alpha_threshold     = 128L,
  object_min_area_cm2 = 5,
  weed_filter         = "width",
  leaf_min_width_cm   = 1.5,
  object_order        = "position"
)
```

The function returns a data frame and also writes it to CSV.

---

## Scale calibration

Everything reported in centimetres depends on one number: how many pixels
correspond to one centimetre. There are three ways to establish it.

### 1. Fixed frame (`image_real_cm`)

`image_real_cm = c(width_cm, height_cm)` defines the real-world size of the
**full image frame**. Photograph a known area — for example a 75 × 75 cm quadrat
filling the frame — and keep camera height constant across the experiment.

### 2. Camera and canopy height

A quadrat laid on the **ground** calibrates the ground plane, but leaves sit
above it, closer to the camera, and therefore appear magnified. With the camera
1 m up, a canopy 30 cm tall is measured roughly twice its true area.

```r
camera_height_cm = 100,        # lens to the calibrated plane
canopy_height_cm = 25          # leaf layer above that plane
```

The correction is a single scalar, so it cannot fully resolve a canopy spread
across many heights. Pass a range to have the residual uncertainty reported
alongside the estimate:

```r
canopy_height_cm = c(10, 30)   # corrects at the midpoint,
                               # also fills Object Area Min/Max (cm2)
```

Raising the camera reduces both the size of the correction and its sensitivity
to a mis-measured canopy height. Update `camera_height_cm` whenever you change
the rig.

### 3. Fiducial marker (`scale_marker = "aruco"`)

Place a printed ArUco marker of known size in the frame and the scale is derived
per image, so varying camera distance stops mattering.

```r
scale_marker   = "aruco",
marker_size_cm = 5,            # measure the printed square with a ruler
marker_id      = 0L            # strongly recommended
```

The marker is searched in the **original** image, because background removal
erases it from the transparent PNG. Detection is guarded three ways —
`marker_id`, `marker_min_px`, and `marker_max_dev` — and any rejection falls
back to `image_real_cm` with a `scale WARNING` on the console. Without those
guards, leaf damage or speckle can be misread as a tiny marker and inflate area
by an order of magnitude.

Put the marker at **canopy height**, not on the ground. When a marker supplies
the scale, the `canopy_height_cm` correction is skipped, since the marker
already measures the correct plane.

---

## Removing thin weed leaves

Background removal separates plants from soil, but it cannot tell a crop leaf
from a grass weed leaf growing through the canopy. `weed_filter = "width"`
discards structures narrower than a given physical width, including thin ones
touching a crop leaf.

```r
weed_filter       = "width",
leaf_min_width_cm = 1.5,
debug_weed_mask   = TRUE     # writes <name>_weed.png with the discarded pixels
```

Broad crop leaflets survive; grass blades and petioles do not. Tune the width by
looking at `<name>_weed.png`: if crop tissue appears there, lower the value. The
threshold is given in centimetres, so it stays valid across image resolutions.

---

## Output

Files written to `output_folder`:

| File | Description |
|---|---|
| `<name>_nobg.png` | Transparent PNG returned by the API |
| `<name>_nobg_processed.jpg` | Annotated image for visual verification |
| `<name>_weed.png` | Pixels discarded by `weed_filter` (when `debug_weed_mask = TRUE`) |
| `image_processed.csv` | All measurements, one row per object |
| `phenoshot_params_<timestamp>.json` | Every parameter, package versions, and the scale used per image |

CSV columns:

`File Name`, `Object ID`, `Image Path`, `Object Area (cm2)`,
`Object Area (px)`, `Object Gap Area (cm2)`, `Object Area Min (cm2)`,
`Object Area Max (cm2)`, `Object Perimeter (cm)`, `Object Length (cm)`,
`Object Width (cm)`, `Object % of Image`, `Pixel Area (cm2/px)`,
`Canopy Height (cm)`, `Scale Source`, `Num objects in Image`

| Column | Meaning |
|---|---|
| `Object Gap Area (cm2)` | Background enclosed by the object's outline but **excluded** from `Object Area (cm2)` |
| `Object Area Min/Max (cm2)` | Bracket the estimate when `canopy_height_cm` is a range; equal to `Object Area (cm2)` otherwise |
| `Canopy Height (cm)` | Height actually used for the correction |
| `Scale Source` | `aruco` or `image_real_cm` |

Objects are numbered by area, descending. Set `object_order = "position"` to
number them in reading order instead, which keeps `Object ID` tied to a physical
location across dates — what time-series analysis needs.

### Verifying the measurement

The annotated image is not decorative. Each detected object is drawn in its own
colour, with a semi-transparent fill covering **exactly the pixels counted
toward its area**. By default the outline traces those same pixels, so gaps of
background between leaves appear as holes rather than being swallowed by a
single enclosing contour — what you see is what was measured.

Set `outline_mode = "outer"` for the older external-contour drawing,
`fill_opacity = 0` to keep outlines only, or `distinct_colors = FALSE` to draw
everything in a single colour.

---

## Arguments

### Core

| Argument | Default | Description |
|---|---|---|
| `input_folder` | — | Folder with input images or `_nobg.png` files |
| `output_folder` | — | Folder for results; created if absent |
| `image_real_cm` | `c(20, 20)` | Real image dimensions `c(width_cm, height_cm)` |
| `photoroom_api_key` | `""` | Photoroom API key |
| `photoroom_size` | `NULL` | Output resolution: `"preview"` (0.25 MP), `"medium"` (1.5 MP), `"hd"` (4 MP), `"full"` (36 MP) |

### Detection and filtering

| Argument | Default | Description |
|---|---|---|
| `alpha_threshold` | `10L` | Alpha cutoff (0–255) separating object from background |
| `k_open` | `c(3L, 3L)` | Kernel for morphological opening (noise removal) |
| `k_close` | `c(5L, 5L)` | Kernel for morphological closing (hole filling) |
| `k_close_iter` | `2L` | Closing iterations; `0L` disables closing for the strictest area |
| `min_component_area_px` | `500L` | Minimum component size in pixels, before calibration |
| `object_min_area_cm2` | `0.3` | Minimum object area in cm² |
| `rel_min_frac_of_largest` | `0.05` | Fallback threshold relative to the largest object |
| `max_keep` | `30L` | Maximum objects retained per image |
| `label_uses_aabb` | `TRUE` | `TRUE` = axis-aligned bounding box; `FALSE` = rotated minimum-area rectangle |
| `object_order` | `"area"` | `"area"` = largest first; `"position"` = reading order |

### Weed removal

| Argument | Default | Description |
|---|---|---|
| `weed_filter` | `"none"` | `"width"` removes structures narrower than `leaf_min_width_cm` |
| `leaf_min_width_cm` | `1.5` | Minimum local width (cm) a structure must reach to be kept |
| `debug_weed_mask` | `FALSE` | Write `<name>_weed.png` showing the discarded pixels |

### Scale

| Argument | Default | Description |
|---|---|---|
| `camera_height_cm` | `NULL` | Lens height above the calibrated plane; `NULL` disables the canopy correction |
| `canopy_height_cm` | `0` | Leaf layer height; length 1, or `c(min, max)` to also report Min/Max |
| `scale_marker` | `"none"` | `"aruco"` derives the scale per image from a detected marker |
| `marker_size_cm` | `5` | Printed marker side length; measure it, don't assume |
| `aruco_dict` | `"DICT_4X4_50"` | Predefined ArUco dictionary name |
| `marker_id` | `NULL` | Accept only this marker ID |
| `marker_min_px` | `40L` | Reject detections smaller than this as false positives |
| `marker_max_dev` | `2` | Reject a marker-derived field of view differing from `image_real_cm` by more than this factor |

### Annotation

| Argument | Default | Description |
|---|---|---|
| `distinct_colors` | `TRUE` | Colour each object differently |
| `contour_color` | `c(0L, 180L, 255L)` | BGR colour used when `distinct_colors = FALSE` |
| `outline_thickness` | `0L` | Outline width in px; `0L` auto-scales to resolution |
| `outline_mode` | `"exact"` | `"exact"` traces the counted pixels; `"outer"` draws the external contour only |
| `outline_min_hole_cm2` | `1` | Smallest internal gap that gets outlined |
| `fill_opacity` | `0.25` | Opacity of the fill over counted pixels; `0` disables |
| `show_crosshair` | `FALSE` | Draw centre crosshairs |

### Reproducibility

| Argument | Default | Description |
|---|---|---|
| `params_json` | `TRUE` | Write `phenoshot_params_<timestamp>.json` with every setting and per-image scale |

---

## Notes

**`alpha_threshold = 128L` is recommended.** The default `10L` is kept for
backward compatibility, but it counts almost all of the semi-transparent halo
that background removal leaves along every leaf edge. On a canopy with a long
perimeter this can inflate area by several percent, and the inflation grows with
edge complexity — a systematic bias between genotypes, not random noise.
Whatever value you pick, keep it constant across an experiment.

**`Object % of Image` needs no calibration.** It is a pixel ratio, so it is
unaffected by `camera_height_cm`, `canopy_height_cm`, or camera distance. For
canopy work it is the most robust metric available, and multiplying it by the
quadrat area gives cm² over a fixed ground plot:

```r
df$cover_pct     <- df$`Object % of Image`
df$plot_area_cm2 <- df$cover_pct / 100 * (75 * 75)
```

Report this as **canopy cover** or **projected canopy area** — not leaf area or
LAI. Overlapping leaves are counted once.

**Existing `_nobg.png` files are reused.** If `output_folder` already contains a
`_nobg.png` for an image, that file is used and no API call is made. To
re-process with different API settings, delete the old files or write to a new
output folder.

**Length and width depend on `label_uses_aabb`.** With `TRUE` (default), width
and length come from the axis-aligned bounding box, so a tilted object reports
larger values than its true dimensions. Set `FALSE` to use the rotated
minimum-area rectangle, which is orientation-independent.

**Cropping is never requested.** The Photoroom `crop` option is deliberately
not used, because cropping changes the image framing and would invalidate the
`image_real_cm` calibration.

**HEIC is uploaded unchanged.** The Remove Background API accepts HEIC
directly, so no re-encoding is performed and edge detail is preserved. Files
above 45 MB are converted to JPEG before upload.

**Touching objects merge.** Plants whose canopies overlap are detected as one
contour and reported as a single object. Check the annotated image and adjust
spacing at photography time if this matters.

**The field of view narrows with height.** When the quadrat is on the ground and
the canopy is above it, the image covers less real area at canopy level than at
the ground — 52.5 cm across a 75 cm quadrat for a 30 cm canopy photographed from
1 m. Corrected `Object Area (cm2)` is the true area of the leaves visible in the
image, not the leaf area over the whole quadrat. Use `Object % of Image` when
you need a value tied to a fixed ground area.

---

## Contact

- GitHub: [agronomy4future/phenoshot](https://github.com/agronomy4future/phenoshot)
- Email: kimjk@agronomy4future.com

© J.K. Kim
