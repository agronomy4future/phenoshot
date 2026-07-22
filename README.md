# phenoshot

**AI-based plant image analysis for morphological trait measurement**

`phenoshot` is an R package that measures plant morphological traits from
photographs. It combines AI background removal (Photoroom API) with
OpenCV-based object detection to report **area, perimeter, length, and width**
in real-world units (cm), and writes an annotated image so every measured
object can be visually verified.

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

---

## API key

Background removal uses the [Photoroom Remove Background API](https://www.photoroom.com/api).
Get a key from the Photoroom dashboard and pass it as `photoroom_api_key`.

A key is **not** required when pre-processed `_nobg.png` files already exist —
see *Mode 1* below.

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
  photoroom_api_key = "your_api_key_here",
  photoroom_size    = "full",
  distinct_colors   = TRUE,
  fill_opacity      = 0.25
)
```

The function returns a data frame and also writes it to CSV.

---

## Calibration

`image_real_cm = c(width_cm, height_cm)` defines the real-world size of the
**full image frame**, and is what converts pixels to centimetres. Measurement
accuracy depends entirely on this value being correct, so photograph a fixed,
known frame area (for example a 75 × 75 cm quadrat filling the frame) and keep
camera height constant across the experiment.

---

## Output

Files written to `output_folder`:

| File | Description |
|---|---|
| `<name>_nobg.png` | Transparent PNG returned by the API |
| `<name>_nobg_processed.jpg` | Annotated image for visual verification |
| `image_processed.csv` | All measurements, one row per object |

CSV columns:

`File Name`, `Object ID`, `Image Path`, `Object Area (cm2)`,
`Object Area (px)`, `Object Perimeter (cm)`, `Object Length (cm)`,
`Object Width (cm)`, `Object % of Image`, `Pixel Area (cm2/px)`,
`Num objects in Image`

Objects are numbered by area, descending — `Obj 1` is the largest.

### Verifying the measurement

The annotated image is not decorative. Each detected object is drawn in its own
colour with a thick outline tracing its **actual contour**, and a
semi-transparent fill covering **exactly the pixels counted toward its area**.
Labels are coloured to match their outline. This makes it immediately visible
whether the segmentation captured the object correctly, whether neighbouring
plants were merged, and how cleanly the background removal cut the edges.

Set `fill_opacity = 0` to keep outlines only, or `distinct_colors = FALSE` to
draw everything in a single colour.

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
| `min_component_area_px` | `500L` | Minimum component size in pixels, before calibration |
| `object_min_area_cm2` | `0.3` | Minimum object area in cm² |
| `rel_min_frac_of_largest` | `0.05` | Fallback threshold relative to the largest object |
| `max_keep` | `30L` | Maximum objects retained per image |
| `label_uses_aabb` | `TRUE` | `TRUE` = axis-aligned bounding box; `FALSE` = rotated minimum-area rectangle |

### Annotation

| Argument | Default | Description |
|---|---|---|
| `distinct_colors` | `TRUE` | Colour each object differently |
| `contour_color` | `c(0L, 180L, 255L)` | BGR colour used when `distinct_colors = FALSE` |
| `outline_thickness` | `0L` | Outline width in px; `0L` auto-scales to resolution |
| `fill_opacity` | `0.25` | Opacity of the fill over counted pixels; `0` disables |
| `show_crosshair` | `FALSE` | Draw centre crosshairs |

---

## Notes

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

---

## Contact

- GitHub: [agronomy4future/phenoshot](https://github.com/agronomy4future/phenoshot)
- Email: kimjk@agronomy4future.com

© J.K. Kim
