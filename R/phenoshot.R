# --- dependencies -------------------------------------------------------
if (!requireNamespace("reticulate", quietly = TRUE)) install.packages("reticulate")
library(reticulate)

ensure_pydeps= function(verbose = TRUE) {
  invisible(reticulate::py_available(initialize = TRUE))
  try({
    reticulate::py_run_string("
import sys, subprocess
subprocess.check_call([sys.executable, '-m', 'pip', 'install', '--upgrade', 'pip', 'setuptools', 'wheel'])
")
  }, silent = TRUE)
  for (pkg in c("numpy", "pandas", "requests", "pillow")) {
    if (!reticulate::py_module_available(pkg)) {
      if (verbose) message("Installing '", pkg, "' ...")
      reticulate::py_install(pkg, pip = TRUE)
    }
  }
  if (!reticulate::py_module_available("pillow_heif")) {
    if (verbose) message("Installing 'pillow-heif' for HEIC support ...")
    reticulate::py_install("pillow-heif", pip = TRUE)
  }
  if (!reticulate::py_module_available("cv2")) {
    try(reticulate::py_install("opencv-python-headless", pip = TRUE), silent = TRUE)
  }
  if (!reticulate::py_module_available("cv2")) stop("cv2 not available.")
  if (verbose) {
    ver= reticulate::py_eval("(__import__('cv2')).__version__")
    message("Dependencies ready. cv2 version: ", ver)
  }
  invisible(TRUE)
}
#' AI-Based Plant Image Analysis for Morphological Trait Measurement
#'
#' Automated pipeline for plant phenotyping through image analysis.
#' Uses AI-based background removal combined with OpenCV-based object
#' detection to measure morphological traits including area, perimeter,
#' length, and width. Optional steps discard thin non-target structures such
#' as grass weed leaves, and correct the pixel-to-centimeter scale either from
#' a detected ArUco marker or from the camera and canopy heights. Supports
#' JPG, PNG, WEBP, and HEIC image formats.
#'
#' @param input_folder Character. Path to folder containing input images
#'   (JPG, PNG, WEBP, HEIC) or pre-processed `_nobg.png` files.
#' @param output_folder Character. Path to folder for saving results.
#'   Created automatically if it does not exist. Output files include
#'   `_nobg.png`, `_nobg_processed.jpg`, and `image_processed.csv`.
#' @param image_real_cm Numeric vector of length 2. Real-world image dimensions
#'   in centimeters as `c(width_cm, height_cm)`. Used to convert pixel
#'   measurements to cm. Default: `c(20, 20)`.
#' @param photoroom_api_key Character. Photoroom API key. Required only when
#'   `_nobg.png` files are not pre-existing. Prefixing the key with
#'   `sandbox_` enables Photoroom's free sandbox mode, but results are
#'   watermarked and must not be used for measurement. Default: `""`.
#' @param photoroom_size Character or NULL. Optional Photoroom `size`
#'   parameter controlling output resolution: `"preview"` (0.25 MP),
#'   `"medium"` (1.5 MP), `"hd"` (4 MP), or `"full"` (36 MP). `"full"` is
#'   recommended for measurement work. `NULL` uses the service default.
#'   Default: `NULL`.
#' @param min_component_area_px Integer. Minimum connected-component area in
#'   pixels to retain before calibration. Default: `500L`.
#' @param k_open Integer vector of length 2. Kernel size `c(w, h)` for
#'   morphological opening (noise removal). Default: `c(3L, 3L)`.
#' @param k_close Integer vector of length 2. Kernel size `c(w, h)` for
#'   morphological closing (hole filling). Default: `c(5L, 5L)`.
#' @param object_min_area_cm2 Numeric. Minimum object area in cm2 after
#'   calibration. Default: `0.3`.
#' @param rel_min_frac_of_largest Numeric between 0 and 1. Minimum size of an
#'   object relative to the largest detected object. Default: `0.05`.
#' @param max_keep Integer. Maximum number of objects to retain per image,
#'   sorted by area descending. Default: `30L`.
#' @param label_uses_aabb Logical. If `TRUE`, length and width are derived from
#'   the axis-aligned bounding box. If `FALSE`, from the rotated minimum-area
#'   rectangle. Default: `TRUE`.
#' @param contour_color Integer vector of length 3. Contour color in BGR format
#'   `c(B, G, R)`, used only when `distinct_colors = FALSE`.
#'   Default: `c(0L, 180L, 255L)` (orange).
#' @param alpha_threshold Integer. Alpha channel threshold between 0 and 255
#'   for object detection. Default: `10L`.
#' @param fill_opacity Numeric between 0 and 1. Opacity of semi-transparent
#'   fill applied to the exact pixels counted for each object. 0 means no
#'   fill. Default: `0.25`.
#' @param distinct_colors Logical. If `TRUE`, each measured object is drawn in a
#'   different color from an internal palette so adjacent objects can be told
#'   apart. If `FALSE`, all objects use `contour_color`. Default: `TRUE`.
#' @param outline_thickness Integer. Thickness (px) of the per-object contour
#'   outline. `0L` = auto-scale to image resolution. Default: `0L`.
#' @param show_crosshair Logical. If `TRUE`, also draws the horizontal/vertical
#'   center crosshair for each object. Default: `FALSE`.
#' @param weed_filter Character. Strategy for discarding thin, non-target
#'   structures such as grass-type weed leaves that survive background removal.
#'   `"none"` disables the step and reproduces the behaviour of earlier
#'   versions. `"width"` applies a distance-transform width filter that removes
#'   any structure narrower than `leaf_min_width_cm`, including thin structures
#'   physically touching a target object. Default: `"none"`.
#' @param leaf_min_width_cm Numeric. Minimum local width, in centimeters, that a
#'   structure must reach to be retained when `weed_filter = "width"`. Narrower
#'   structures such as grass weed leaves and petioles are removed. Specified in
#'   cm rather than pixels so the setting stays valid across image resolutions.
#'   Default: `1.5`.
#' @param debug_weed_mask Logical. If `TRUE`, writes `<stem>_weed.png`
#'   containing the pixels discarded by `weed_filter`, so `leaf_min_width_cm`
#'   can be tuned by visual inspection. Default: `FALSE`.
#' @param outline_mode Character. How the per-object outline is drawn.
#'   `"exact"` traces the boundary of the pixels actually counted, so internal
#'   background gaps are outlined as holes and the drawing matches the reported
#'   area. `"outer"` draws only the external contour, which visually encloses
#'   background gaps that are not counted; this reproduces the behaviour of
#'   earlier versions. Default: `"exact"`.
#' @param outline_min_hole_cm2 Numeric. Minimum area (cm2) an internal gap must
#'   have to be outlined when `outline_mode = "exact"`. Prevents insect-damage
#'   holes and single-pixel gaps from cluttering the overlay. Default: `1`.
#' @param k_close_iter Integer. Number of iterations for the morphological
#'   closing that fills small gaps. Higher values merge more background between
#'   touching objects into the measured area; `0L` disables closing entirely for
#'   the strictest area estimate. Default: `2L`.
#' @param scale_marker Character. Source of the pixel-to-centimeter scale.
#'   `"none"` derives it from `image_real_cm`, which assumes the camera
#'   distance is identical for every image. `"aruco"` detects an ArUco marker of
#'   known size in the ORIGINAL image and derives the scale per image, so
#'   variation in camera height no longer biases area. Requires
#'   `opencv-contrib-python`; falls back to `image_real_cm` with a message when
#'   no marker is found, or when the input folder holds `_nobg.png` files (the
#'   marker has already been erased by background removal). Default: `"none"`.
#' @param marker_size_cm Numeric. Physical side length of the ArUco marker in
#'   centimeters. Default: `5`.
#' @param aruco_dict Character. Name of the predefined ArUco dictionary, e.g.
#'   `"DICT_4X4_50"`, `"DICT_5X5_100"`, `"DICT_6X6_250"`. Default:
#'   `"DICT_4X4_50"`.
#' @param marker_id Integer or `NULL`. If given, only a marker carrying this ID
#'   is accepted, which is the strongest protection against leaf damage or
#'   speckle being misread as a marker. `NULL` accepts any ID. Default: `NULL`.
#' @param marker_min_px Integer. Minimum marker side length in pixels. Smaller
#'   detections are rejected as false positives. Default: `40L`.
#' @param marker_max_dev Numeric. A marker-derived field of view is rejected
#'   when it differs from `image_real_cm` by more than this factor in either
#'   direction, and `image_real_cm` is used instead. Set larger only if the
#'   fallback value is a rough guess. Default: `2`.
#' @param camera_height_cm Numeric or `NULL`. Height of the camera above the
#'   plane on which `image_real_cm` was measured, usually the ground. Supplying
#'   it together with `canopy_height_cm` corrects the magnification caused by
#'   leaves sitting above that plane. `NULL` disables the correction. The
#'   correction is skipped when a marker supplied the scale, since a marker
#'   placed at canopy height already measures the correct plane.
#'   Default: `NULL`.
#' @param canopy_height_cm Numeric of length 1 or 2. Height of the leaf layer
#'   above the ground at the time of capture. A single value applies one
#'   correction. A length-2 vector `c(min, max)` corrects using the midpoint and
#'   additionally reports Object Area Min (cm2) and Object Area Max (cm2), which
#'   bracket the estimate for a canopy spanning that range. Default: `0`.
#' @param object_order Character. Rule for assigning Object ID. `"area"` numbers
#'   objects from largest to smallest, so IDs change whenever an object grows.
#'   `"position"` numbers them in reading order (top-left to bottom-right),
#'   keeping ID to physical location stable across dates, which is what
#'   time-series analysis needs. Default: `"area"`.
#' @param params_json Logical. If `TRUE`, writes
#'   `phenoshot_params_<timestamp>.json` into `output_folder` recording every
#'   parameter, package versions, and the scale actually used for each image.
#'   Default: `TRUE`.
#'
#' @return A data frame with one row per detected object containing columns
#'   File Name, Object ID, Image Path, Object Area (cm2), Object Area (px),
#'   Object Gap Area (cm2), Object Area Min (cm2), Object Area Max (cm2),
#'   Object Perimeter (cm), Object Length (cm), Object Width (cm),
#'   Object percent of Image, Pixel Area (cm2/px), Canopy Height (cm),
#'   Scale Source, and Num objects in Image.
#'
#'   Object Gap Area (cm2) is the background enclosed by the object's external
#'   contour but excluded from Object Area (cm2). Object Area Min (cm2) and
#'   Object Area Max (cm2) bracket the estimate when `canopy_height_cm` is given
#'   as a range, and equal Object Area (cm2) otherwise. Canopy Height (cm) is
#'   the height actually used for the correction, and Scale Source records
#'   whether the scale came from a detected marker or from image_real_cm.
#'
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # Install the package
#' if(!require(remotes)) install.packages("remotes")
#' if (!requireNamespace("phenoshot", quietly= TRUE)) {
#'   remotes::install_github("agronomy4future/phenoshot", force= TRUE)
#' }
#' library(remotes)
#' library(phenoshot)
#'
#' # Case 1: Using _nobg.png files in Input folder (no API key required)
#' phenoshot(
#'   input_folder  = "./Input",
#'   output_folder = "./Output",
#'   image_real_cm = c(20, 20)
#' )
#'
#' # Case 2: Using original images with the Photoroom API
#' phenoshot(
#'   input_folder      = "./Input",
#'   output_folder     = "./Output",
#'   image_real_cm     = c(20, 20),
#'   photoroom_api_key = "your_api_key_here",
#'   photoroom_size    = "full",
#'   min_component_area_px = 500L,
#'   object_min_area_cm2   = 0.3,
#'   max_keep              = 50L,
#'   distinct_colors       = TRUE,
#'   outline_thickness     = 0L,
#'   fill_opacity          = 0.25,
#'   show_crosshair        = FALSE
#' )
#'
#' # Case 3: Removing thin grass-type weed leaves from a crop canopy
#' phenoshot(
#'   input_folder      = "./Input",
#'   output_folder     = "./Output",
#'   image_real_cm     = c(60, 60),
#'   weed_filter       = "width",
#'   leaf_min_width_cm = 1.5,
#'   debug_weed_mask   = TRUE   # inspect <stem>_weed.png, then tune the width
#' )
#'
#' # Case 4: Frame on the ground, camera above it, canopy in between
#' phenoshot(
#'   input_folder      = "./Input",
#'   output_folder     = "./Output",
#'   image_real_cm     = c(75, 75),   # frame inner size, measured at the ground
#'   camera_height_cm  = 100,         # lens to ground, measured
#'   canopy_height_cm  = c(10, 30),   # leaf layer above the ground
#'   alpha_threshold   = 128L,
#'   object_order      = "position"
#' )
#'
#' # Github: https://github.com/agronomy4future/phenoshot
#' # All Rights Reserved (c) J.K Kim (kimjk@agronomy4future.com)
#' }
phenoshot= function(
    input_folder,
    output_folder,
    image_real_cm           = c(20, 20),
    photoroom_api_key       = "",
    photoroom_size          = NULL,
    min_component_area_px   = 500L,
    k_open                  = c(3L, 3L),
    k_close                 = c(5L, 5L),
    object_min_area_cm2     = 0.3,
    rel_min_frac_of_largest = 0.05,
    max_keep                = 30L,
    label_uses_aabb         = TRUE,
    contour_color           = c(0L, 180L, 255L),
    alpha_threshold         = 10L,
    fill_opacity            = 0.25,
    distinct_colors         = TRUE,
    outline_thickness       = 0L,
    show_crosshair          = FALSE,
    weed_filter             = c("none", "width"),
    leaf_min_width_cm       = 1.5,
    debug_weed_mask         = FALSE,
    outline_mode            = c("exact", "outer"),
    outline_min_hole_cm2    = 1,
    k_close_iter            = 2L,
    scale_marker            = c("none", "aruco"),
    marker_size_cm          = 5,
    aruco_dict              = "DICT_4X4_50",
    marker_id               = NULL,
    marker_min_px           = 40L,
    marker_max_dev          = 2,
    camera_height_cm        = NULL,
    canopy_height_cm        = 0,
    object_order            = c("area", "position"),
    params_json             = TRUE
) {
  if (length(image_real_cm) == 1) image_real_cm= rep(image_real_cm, 2)
  stopifnot(length(image_real_cm) == 2, all(is.finite(image_real_cm)), all(image_real_cm > 0))

  weed_filter= match.arg(weed_filter)
  outline_mode= match.arg(outline_mode)
  stopifnot(length(leaf_min_width_cm) == 1, is.finite(leaf_min_width_cm), leaf_min_width_cm >= 0)
  stopifnot(length(outline_min_hole_cm2) == 1, is.finite(outline_min_hole_cm2),
            outline_min_hole_cm2 >= 0)
  stopifnot(length(k_close_iter) == 1, is.finite(k_close_iter), k_close_iter >= 0)
  scale_marker= match.arg(scale_marker)
  object_order= match.arg(object_order)
  stopifnot(length(marker_size_cm) == 1, is.finite(marker_size_cm), marker_size_cm > 0)
  stopifnot(length(aruco_dict) == 1, is.character(aruco_dict))
  stopifnot(is.null(marker_id) || (length(marker_id) == 1 && is.finite(marker_id)))
  stopifnot(length(marker_min_px) == 1, is.finite(marker_min_px), marker_min_px > 0)
  stopifnot(length(marker_max_dev) == 1, is.finite(marker_max_dev), marker_max_dev > 1)
  if (!is.null(camera_height_cm)) {
    stopifnot(length(camera_height_cm) == 1, is.finite(camera_height_cm), camera_height_cm > 0)
    stopifnot(length(canopy_height_cm) %in% 1:2, all(is.finite(canopy_height_cm)),
              all(canopy_height_cm >= 0))
    if (max(canopy_height_cm) >= camera_height_cm)
      stop("canopy_height_cm must be below camera_height_cm.", call. = FALSE)
  }
  if (weed_filter == "width" && leaf_min_width_cm >= min(image_real_cm) / 2) {
    warning("leaf_min_width_cm (", leaf_min_width_cm, ") is large relative to ",
            "image_real_cm; most or all objects may be removed.", call. = FALSE)
  }

  if (grepl("^sandbox_", photoroom_api_key)) {
    warning("Photoroom sandbox key detected: results are watermarked and are ",
            "not suitable for measurement. Use a live key for real analysis.",
            call. = FALSE)
  }

  if (!dir.exists(output_folder)) dir.create(output_folder, recursive = TRUE, showWarnings = FALSE)

  ensure_pydeps(verbose = FALSE)

  py_code = "
import cv2
import numpy as np
import os, glob, pandas as pd
import json
from datetime import datetime
import requests
from PIL import Image
import io

try:
    from pillow_heif import register_heif_opener
    register_heif_opener()
except Exception:
    pass

# Photoroom Remove Background API
PHOTOROOM_URL     = 'https://sdk.photoroom.com/v1/segment'
PHOTOROOM_HEADER  = 'x-api-key'
PHOTOROOM_MAX_MB  = 45     # fall back to JPEG above this upload size
QUOTA_HDRS        = ['x-credits-remaining', 'x-remaining-credits', 'X-Credits-Remaining']

# Visually distinct BGR colors used when distinct_colors=True
PALETTE = [
    (0,   0,   255),   # red
    (0,   140, 255),   # orange
    (0,   215, 255),   # amber
    (0,   220, 0),     # green
    (220, 220, 0),     # cyan
    (255, 120, 0),     # azure
    (255, 0,   0),     # blue
    (255, 0,   200),   # magenta
    (180, 0,   255),   # pink
    (0,   180, 140),   # olive
    (120, 90,  255),   # coral
    (0,   255, 120),   # spring green
]

def _discover_unique_images(folder):
    patterns = [
        '*.jpg','*.jpeg','*.png','*.webp',
        '*.JPG','*.JPEG','*.PNG','*.WEBP',
        '*.heic','*.HEIC','*.heif','*.HEIF'
    ]
    all_paths = []
    for p in patterns:
        all_paths.extend(glob.glob(os.path.join(folder, p)))
    seen = {}
    for p in all_paths:
        key = os.path.normcase(os.path.abspath(p))
        if key not in seen:
            seen[key] = p
    unique_paths = sorted(seen.values())

    # if _nobg.png exists -> use only those (no API / no original needed)
    # otherwise -> use original images
    nobg_files = [p for p in unique_paths
                  if '_nobg' in os.path.splitext(os.path.basename(p))[0]]
    orig_files = [p for p in unique_paths
                  if '_nobg' not in os.path.splitext(os.path.basename(p))[0]]

    if nobg_files:
        print(f'Found {len(nobg_files)} _nobg image(s) -> API skipped.')
        return nobg_files, True   # True = nobg mode
    print(f'Found {len(orig_files)} original image(s).')
    return orig_files, False      # False = original mode

def _open_image(path):
    return Image.open(path).convert('RGB')

def _build_files(image_path):
    '''Prepare the multipart payload for Photoroom.

    Photoroom accepts PNG, JPEG, WEBP and HEIC, so the original file is
    uploaded untouched. This avoids a lossy downscale + re-encode that
    would blur object edges and bias the area measurement. Only files
    above the upload size limit are converted to JPEG.
    '''
    ext     = os.path.splitext(image_path)[1].lower()
    size_mb = os.path.getsize(image_path) / (1024 * 1024)

    if size_mb <= PHOTOROOM_MAX_MB:
        ctype = {'.heic': 'image/heic', '.heif': 'image/heic',
                 '.png' : 'image/png',  '.webp': 'image/webp',
                 '.jpg' : 'image/jpeg', '.jpeg': 'image/jpeg'}.get(ext)
        if ext in ('.heic', '.heif'):
            print(f'  Uploading HEIC natively ({size_mb:.1f} MB, no re-encode).')
        with open(image_path, 'rb') as f:
            data = f.read()
        name = os.path.basename(image_path)
        if ctype:
            return {'image_file': (name, data, ctype)}
        return {'image_file': (name, data)}

    # oversized -> downscale and re-encode as JPEG
    print(f'  Image too large ({size_mb:.1f} MB) -> converting to JPEG.')
    img = _open_image(image_path)
    ow, oh = img.size
    max_px = 5000
    if max(ow, oh) > max_px:
        scale = max_px / max(ow, oh)
        img = img.resize((int(ow*scale), int(oh*scale)), Image.LANCZOS)
    buf = io.BytesIO()
    img.save(buf, format='JPEG', quality=92)
    buf.seek(0)
    return {'image_file': ('image.jpg', buf, 'image/jpeg')}

def _remove_background(image_path, api_key, nobg_save_path, photoroom_size):
    '''Call the Photoroom Remove Background API and save the RGBA PNG.'''
    files = _build_files(image_path)

    # png output keeps the alpha channel, which the pipeline relies on.
    # NOTE: 'crop' is deliberately NOT used -- cropping would change the
    # framing and invalidate the image_real_cm calibration.
    data = {'format': 'png'}
    if photoroom_size:
        data['size'] = str(photoroom_size)

    response = requests.post(
        PHOTOROOM_URL,
        files=files,
        data=data,
        headers={PHOTOROOM_HEADER: api_key},
        timeout=120
    )
    if response.status_code != 200:
        raise RuntimeError(
            f'Photoroom API error {response.status_code}: {response.text[:300]}'
        )

    quota = None
    for hdr in QUOTA_HDRS:
        if hdr in response.headers:
            quota = f'{hdr}={response.headers[hdr]}'
            break
    print('  Photoroom OK.' + (f' ({quota})' if quota else ''))

    with open(nobg_save_path, 'wb') as f:
        f.write(response.content)
    print(f'  nobg saved: {nobg_save_path}')
    return response.content

def _perimeter_cm_from_cnt(cnt, sx, sy):
    c = cnt.reshape(-1, 2).astype(np.float64)
    c_scaled = np.empty_like(c)
    c_scaled[:,0] = c[:,0] * sx
    c_scaled[:,1] = c[:,1] * sy
    d = np.diff(np.vstack([c_scaled, c_scaled[0]]), axis=0)
    return float(np.sqrt((d[:,0]**2) + (d[:,1]**2)).sum())

def _detect_scale_px_per_cm(bgr_or_rgb, dict_name, marker_size_cm,
                            marker_id=None, marker_min_px=40):
    '''Return (px_per_cm, note). px_per_cm is None when no marker is usable.'''
    if not hasattr(cv2, 'aruco'):
        return None, 'cv2.aruco unavailable (install opencv-contrib-python)'
    key = getattr(cv2.aruco, str(dict_name), None)
    if key is None:
        return None, f'unknown dictionary {dict_name}'
    try:
        dic = cv2.aruco.getPredefinedDictionary(key)
    except Exception as e:
        return None, f'dictionary error: {e}'
    gray = cv2.cvtColor(bgr_or_rgb, cv2.COLOR_RGB2GRAY)
    corners, ids = None, None
    try:                                        # OpenCV >= 4.7
        det = cv2.aruco.ArucoDetector(dic, cv2.aruco.DetectorParameters())
        corners, ids, _ = det.detectMarkers(gray)
    except AttributeError:                      # OpenCV < 4.7
        try:
            corners, ids, _ = cv2.aruco.detectMarkers(
                gray, dic, parameters=cv2.aruco.DetectorParameters_create())
        except Exception as e:
            return None, f'detector error: {e}'
    if not corners:
        return None, 'no marker detected'
    if float(marker_size_cm) <= 0:
        return None, 'marker_size_cm must be > 0'

    id_list = [] if ids is None else [int(v) for v in np.asarray(ids).ravel()]
    n_found = len(corners)

    # keep only the requested ID, when one is specified
    if marker_id is not None and id_list:
        want = int(marker_id)
        pairs = [(c, i) for c, i in zip(corners, id_list) if i == want]
        if not pairs:
            return None, f'marker id {want} not among detected ids {id_list}'
        corners = [p[0] for p in pairs]

    # median side length per marker, then reject implausibly small patterns
    per_marker = []
    for c in corners:
        p = np.asarray(c).reshape(4, 2).astype(np.float64)
        s = [float(np.linalg.norm(p[i] - p[(i + 1) % 4])) for i in range(4)]
        per_marker.append((float(np.median(s)), float(np.max(s)), float(np.min(s))))
    med_side = float(np.median([m[0] for m in per_marker]))
    if med_side < float(marker_min_px):
        return None, (f'rejected: marker side {med_side:.0f} px < marker_min_px '
                      f'({marker_min_px}) -- almost certainly a false positive '
                      f'(found {n_found}, ids {id_list})')

    ppc    = med_side / float(marker_size_cm)
    hi, lo = max(m[1] for m in per_marker), min(m[2] for m in per_marker)
    spread = (hi - lo) / max(1e-9, med_side)
    note = (f'{len(corners)} marker(s) ids {id_list}, side {med_side:.0f} px, '
            f'{ppc:.2f} px/cm, spread {100*spread:.1f}%')
    if spread > 0.15:
        note += ' -- WARNING: large spread, camera may be tilted'
    return ppc, note


def process_images(input_folder, output_folder,
                   image_real_cm_W, image_real_cm_H,
                   photoroom_api_key, photoroom_size,
                   min_component_area_px, k_open, k_close,
                   object_min_area_cm2, rel_min_frac_of_largest, max_keep,
                   label_uses_aabb, contour_color, alpha_threshold, fill_opacity,
                   distinct_colors, outline_thickness, show_crosshair,
                   weed_filter='none', leaf_min_width_cm=0.0,
                   debug_weed_mask=False, outline_mode='exact',
                   outline_min_hole_cm2=1.0, k_close_iter=2,
                   scale_marker='none', marker_size_cm=5.0,
                   aruco_dict='DICT_4X4_50', marker_id=None,
                   marker_min_px=40, marker_max_dev=2.0,
                   object_order='area', params_json=True,
                   camera_height_cm=None, canopy_height_cm=0.0):

    if not os.path.isdir(input_folder):
        print(f'Input folder does not exist: {input_folder}')
        return pd.DataFrame()

    os.makedirs(output_folder, exist_ok=True)
    image_paths, is_nobg_mode = _discover_unique_images(input_folder)
    if not image_paths:
        print('No images found.')
        return pd.DataFrame()

    col  = tuple(int(c) for c in contour_color)
    rows = []
    run_images = []

    for path in image_paths:
        filename = os.path.basename(path)
        base     = os.path.splitext(filename)[0]

        if is_nobg_mode:
            # direct _nobg.png input mode
            # C1_nobg.png -> stem = C1
            stem      = base.replace('_nobg', '')
            nobg_path = path   # the file itself is the nobg
        else:
            # original image mode
            stem      = base
            nobg_path = os.path.join(output_folder, stem + '_nobg.png')

        proc_path = os.path.join(output_folder, stem + '_nobg_processed.jpg')
        print(f'Processing: {filename}')

        # -- 1. obtain the transparent PNG ---------------------------
        if is_nobg_mode:
            with open(nobg_path, 'rb') as f:
                png_bytes = f.read()
            print(f'  Using _nobg directly: {nobg_path}')
        elif os.path.exists(nobg_path):
            print(f'  Using existing nobg: {nobg_path}')
            with open(nobg_path, 'rb') as f:
                png_bytes = f.read()
        elif photoroom_api_key:
            try:
                png_bytes = _remove_background(path, photoroom_api_key,
                                               nobg_path, photoroom_size)
            except Exception as e:
                print(f'  Photoroom failed: {e} -- skipping')
                continue
        else:
            print('  No _nobg.png and no API key -- skipping')
            continue

        # -- 2. alpha channel -> mask --------------------------------
        rgba  = np.array(Image.open(io.BytesIO(png_bytes)).convert('RGBA'))
        alpha = rgba[:, :, 3]
        rgb   = rgba[:, :, :3]

        # -- 2a. scale: fiducial marker overrides image_real_cm ------
        # The marker is searched in the ORIGINAL image, because background
        # removal deletes it from the transparent PNG.
        img_cm_W     = float(image_real_cm_W)
        img_cm_H     = float(image_real_cm_H)
        scale_source = 'image_real_cm'
        scale_note   = 'marker detection disabled'
        if str(scale_marker) == 'aruco':
            if is_nobg_mode:
                scale_note = 'nobg input mode: original image unavailable'
            else:
                try:
                    orig = np.array(Image.open(path).convert('RGB'))
                except Exception as e:
                    orig, scale_note = None, f'could not read original: {e}'
                if orig is not None:
                    ppc, scale_note = _detect_scale_px_per_cm(
                        orig, aruco_dict, marker_size_cm, marker_id, marker_min_px)
                    if ppc:
                        mh, mw = orig.shape[:2]
                        rh, rw = rgba.shape[:2]
                        # Photoroom may resize; rescale px/cm to the returned PNG
                        ppc_r = ppc * (0.5 * (rw / float(mw) + rh / float(mh)))
                        cand_W, cand_H = rw / ppc_r, rh / ppc_r
                        # sanity check against the declared field of view
                        dev = max(cand_W / float(image_real_cm_W),
                                  float(image_real_cm_W) / cand_W)
                        if dev > float(marker_max_dev):
                            scale_note = (f'rejected: implied field of view '
                                          f'{cand_W:.0f} x {cand_H:.0f} cm differs from '
                                          f'image_real_cm by {dev:.1f}x '
                                          f'(limit {marker_max_dev}x) -- {scale_note}')
                        else:
                            img_cm_W, img_cm_H = cand_W, cand_H
                            scale_source = 'aruco'
            if scale_source == 'aruco':
                print(f'  scale: {scale_note} -> field of view {img_cm_W:.1f} x {img_cm_H:.1f} cm')
            else:
                print(f'  scale WARNING: {scale_note}')
                print(f'  scale: falling back to image_real_cm '
                      f'({image_real_cm_W:.1f} x {image_real_cm_H:.1f} cm)')

        # -- 2b. canopy height correction ----------------------------
        # The frame sits on the ground but the leaves are above it, so the
        # leaf plane is closer to the camera and appears magnified. Skipped
        # when a marker supplied the scale -- place the marker at canopy
        # height and it already measures the correct plane.
        canopy_mid   = 0.0
        canopy_lo_r  = 1.0          # area multiplier at the TALLEST canopy
        canopy_hi_r  = 1.0          # area multiplier at the SHORTEST canopy
        canopy_note  = 'not applied'
        if camera_height_cm is not None and scale_source == 'image_real_cm':
            H = float(camera_height_cm)
            try:
                ch = [float(x) for x in canopy_height_cm]
            except TypeError:
                ch = [float(canopy_height_cm)]
            ch = [c for c in ch if c == c]
            h_lo, h_hi = min(ch), max(ch)
            if H <= 0 or h_hi >= H:
                canopy_note = f'skipped: canopy {h_hi} cm not below camera {H} cm'
            elif h_hi <= 0:
                canopy_note = 'canopy height 0 -- no correction needed'
            else:
                canopy_mid = 0.5 * (h_lo + h_hi)
                f_mid      = (H - canopy_mid) / H
                img_cm_W  *= f_mid
                img_cm_H  *= f_mid
                canopy_lo_r = ((H - h_hi) / (H - canopy_mid)) ** 2
                canopy_hi_r = ((H - h_lo) / (H - canopy_mid)) ** 2
                canopy_note = (f'camera {H:.0f} cm, canopy {h_lo:g}-{h_hi:g} cm '
                               f'(mid {canopy_mid:g}) -> area x{f_mid**2:.3f}')
                print(f'  canopy: {canopy_note}, '
                      f'field of view {img_cm_W:.2f} x {img_cm_H:.2f} cm')

        mask = (alpha > int(alpha_threshold)).astype(np.uint8) * 255
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN,
                                np.ones(tuple(k_open),  np.uint8), iterations=1)

        # -- 2b. thin-structure (grass weed) removal -----------------
        # Runs before MORPH_CLOSE, which would otherwise bridge weed leaves
        # onto neighbouring target objects and make them inseparable.
        if str(weed_filter) == 'width' and float(leaf_min_width_cm) > 0:
            hh, ww    = mask.shape[:2]
            px_per_cm = 0.5 * (ww / img_cm_W + hh / img_cm_H)
            r         = max(1, int(round(0.5 * float(leaf_min_width_cm) * px_per_cm)))
            dist = cv2.distanceTransform(mask, cv2.DIST_L2, 5)
            core = (dist >= r).astype(np.uint8) * 255
            if cv2.countNonZero(core) == 0:
                print(f'  weed_filter: nothing reaches {leaf_min_width_cm} cm wide -- step skipped')
            else:
                ker    = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (2*r+1, 2*r+1))
                kept   = cv2.bitwise_and(cv2.dilate(core, ker), mask)
                before = cv2.countNonZero(mask)
                after  = cv2.countNonZero(kept)
                if debug_weed_mask:
                    weed_path = os.path.join(output_folder, stem + '_weed.png')
                    cv2.imwrite(weed_path, cv2.subtract(mask, kept))
                    print(f'  weed_filter: discarded pixels saved to {weed_path}')
                pct = (100.0 * (before - after) / before) if before else 0.0
                print(f'  weed_filter: removed {before-after} px ({pct:.1f}%), radius {r} px')
                mask = kept

        if int(k_close_iter) > 0:
            mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE,
                                    np.ones(tuple(k_close), np.uint8),
                                    iterations=int(k_close_iter))

        h, w = mask.shape[:2]
        area_per_pixel_cm2 = (img_cm_W / w) * (img_cm_H / h)
        sx = img_cm_W / w
        sy = img_cm_H / h

        font_scale = max(0.6, round(max(w, h) / 1600, 2))
        line_thick = max(1, int(max(w, h) / 1500))

        # per-object outline thickness (auto-scale when outline_thickness <= 0)
        ot = int(outline_thickness)
        ol_thick = ot if ot > 0 else max(2, int(max(w, h) / 500))

        # -- 3. contour detection ------------------------------------
        contours, _ = cv2.findContours(mask, cv2.RETR_EXTERNAL, cv2.CHAIN_APPROX_SIMPLE)
        contours = [c for c in contours if cv2.contourArea(c) >= float(min_component_area_px)]

        candidates = [
            (c, float(cv2.contourArea(c)), float(cv2.contourArea(c)) * area_per_pixel_cm2)
            for c in contours
        ]
        kept = [t for t in candidates if t[2] >= object_min_area_cm2]
        if not kept and candidates:
            largest_cm2 = max(t[2] for t in candidates)
            kept = [t for t in candidates if t[2] >= rel_min_frac_of_largest * largest_cm2]
        kept.sort(key=lambda t: t[2], reverse=True)
        kept = kept[:int(max_keep)]
        contours = [t[0] for t in kept]

        # -- 3b. object numbering ------------------------------------
        # 'area' reshuffles IDs whenever an object grows; 'position' keeps
        # ID -> physical location stable across dates for time-series work.
        if str(object_order) == 'position':
            band = max(1, int(round(mask.shape[0] * 0.10)))
            def _pos_key(c):
                bx, by, bw_, bh_ = cv2.boundingRect(c)
                return ((by + bh_ // 2) // band, bx + bw_ // 2)
            contours.sort(key=_pos_key)

        print(f'  Detected: {len(contours)} objects')

        # -- 4. composite over white background ----------------------
        comp = np.ones_like(rgb, dtype=np.uint8) * 255
        comp[alpha > 0] = rgb[alpha > 0]
        annotated = cv2.cvtColor(comp, cv2.COLOR_RGB2BGR)

        if contours:
            for idx, cnt in enumerate(contours, start=1):
                # exact pixels counted toward this object's area
                filled = np.zeros(mask.shape, dtype=np.uint8)
                cv2.drawContours(filled, [cnt], -1, 255, thickness=-1)
                actual_mask = cv2.bitwise_and(mask, filled)
                area_px  = float(cv2.countNonZero(actual_mask))
                area_cm2 = area_px * area_per_pixel_cm2
                gap_cm2  = (float(cv2.countNonZero(filled)) - area_px) * area_per_pixel_cm2
                perim_cm = _perimeter_cm_from_cnt(cnt, sx, sy)

                x, y, bw, bh = cv2.boundingRect(cnt)
                cx = x + bw // 2
                cy = y + bh // 2

                if label_uses_aabb:
                    length_cm = bh * sy
                    width_cm  = bw * sx
                else:
                    pts    = cnt.reshape(-1, 2).astype(np.float32)
                    pts_cm = np.column_stack([pts[:,0]*sx, pts[:,1]*sy]).astype(np.float32)
                    rc     = cv2.minAreaRect(pts_cm)
                    length_cm = float(max(rc[1]))
                    width_cm  = float(min(rc[1]))

                # pick this object's color
                obj_col = PALETTE[(idx - 1) % len(PALETTE)] if bool(distinct_colors) else col

                # (a) semi-transparent fill over the EXACT counted pixels
                if float(fill_opacity) > 0:
                    overlay = annotated.copy()
                    overlay[actual_mask > 0] = obj_col
                    cv2.addWeighted(overlay, float(fill_opacity),
                                    annotated, 1.0 - float(fill_opacity),
                                    0, annotated)

                # (b) outline tracing this object
                if str(outline_mode) == 'exact':
                    # boundary of the pixels actually counted, holes included,
                    # so the drawing matches Object Area (cm2)
                    sub, _ = cv2.findContours(actual_mask, cv2.RETR_CCOMP,
                                              cv2.CHAIN_APPROX_SIMPLE)
                    min_hole_px = float(outline_min_hole_cm2) / area_per_pixel_cm2
                    sub = [sc for sc in sub if cv2.contourArea(sc) >= min_hole_px]
                    cv2.drawContours(annotated, sub, -1, obj_col, ol_thick, cv2.LINE_AA)
                else:
                    cv2.drawContours(annotated, [cnt], -1, obj_col, ol_thick, cv2.LINE_AA)

                # (c) optional center crosshair
                if bool(show_crosshair):
                    cv2.line(annotated, (x, cy), (x+bw, cy), obj_col, line_thick)
                    cv2.line(annotated, (cx, y), (cx, y+bh), obj_col, line_thick)

                # (d) label colored to match the outline (black halo for legibility)
                label = f'Obj {idx} ({area_cm2:.2f} cm2)'
                ty = max(int(font_scale*20)+4, y - 6)
                cv2.putText(annotated, label, (x, ty),
                            cv2.FONT_HERSHEY_SIMPLEX, font_scale,
                            (0,0,0), line_thick+3, cv2.LINE_AA)
                cv2.putText(annotated, label, (x, ty),
                            cv2.FONT_HERSHEY_SIMPLEX, font_scale,
                            obj_col, line_thick+1, cv2.LINE_AA)

                rows.append({
                    'File Name'            : filename,
                    'Object ID'            : idx,
                    'Image Path'           : path,
                    'Object Area (cm2)'    : round(area_cm2, 2),
                    'Object Area (px)'     : int(round(area_px)),
                    'Object Gap Area (cm2)': round(gap_cm2, 2),
                    'Object Area Min (cm2)': round(area_cm2 * canopy_lo_r, 2),
                    'Object Area Max (cm2)': round(area_cm2 * canopy_hi_r, 2),
                    'Canopy Height (cm)'   : round(canopy_mid, 2),
                    'Object Perimeter (cm)': round(perim_cm, 2),
                    'Object Length (cm)'   : round(float(length_cm), 2),
                    'Object Width (cm)'    : round(float(width_cm), 2),
                    'Object % of Image'    : round(100.0 * area_px / (w * h), 2),
                    'Pixel Area (cm2/px)'  : round(area_per_pixel_cm2, 8),
                    'Scale Source'         : scale_source,
                    'Num objects in Image' : len(contours)
                })
        else:
            cv2.putText(annotated, 'No objects detected', (20, 60),
                        cv2.FONT_HERSHEY_SIMPLEX, font_scale, (0,0,255), 2, cv2.LINE_AA)

        cv2.imwrite(proc_path, annotated, [cv2.IMWRITE_JPEG_QUALITY, 95])
        print(f'  Saved: {proc_path}')

        run_images.append({
            'file'          : filename,
            'stem'          : stem,
            'width_px'      : int(w),
            'height_px'     : int(h),
            'fov_cm'        : [round(img_cm_W, 3), round(img_cm_H, 3)],
            'px_per_cm'     : round(0.5 * (w / img_cm_W + h / img_cm_H), 4),
            'scale_source'  : scale_source,
            'scale_note'    : scale_note,
            'canopy_note'   : canopy_note,
            'n_objects'     : len(contours)
        })

    df = pd.DataFrame(rows)
    csv_path = os.path.join(output_folder, 'image_processed.csv')
    if len(df) > 0:
        try:
            df.to_csv(csv_path, index=False, encoding='utf-8-sig')
            print(f'CSV saved: {csv_path}')
        except Exception as e:
            print(f'Failed to save CSV: {e}')

    if bool(params_json):
        stamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        meta = {
            'run_timestamp' : datetime.now().isoformat(timespec='seconds'),
            'versions'      : {'opencv': cv2.__version__, 'numpy': np.__version__,
                               'pandas': pd.__version__},
            'input_folder'  : input_folder,
            'output_folder' : output_folder,
            'csv_path'      : csv_path if len(df) > 0 else None,
            'parameters'    : {
                'image_real_cm'          : [float(image_real_cm_W), float(image_real_cm_H)],
                'min_component_area_px'  : int(min_component_area_px),
                'k_open'                 : [int(v) for v in k_open],
                'k_close'                : [int(v) for v in k_close],
                'k_close_iter'           : int(k_close_iter),
                'object_min_area_cm2'    : float(object_min_area_cm2),
                'rel_min_frac_of_largest': float(rel_min_frac_of_largest),
                'max_keep'               : int(max_keep),
                'label_uses_aabb'        : bool(label_uses_aabb),
                'alpha_threshold'        : int(alpha_threshold),
                'fill_opacity'           : float(fill_opacity),
                'outline_mode'           : str(outline_mode),
                'outline_min_hole_cm2'   : float(outline_min_hole_cm2),
                'weed_filter'            : str(weed_filter),
                'leaf_min_width_cm'      : float(leaf_min_width_cm),
                'scale_marker'           : str(scale_marker),
                'marker_size_cm'         : float(marker_size_cm),
                'aruco_dict'             : str(aruco_dict),
                'marker_id'              : (None if marker_id is None else int(marker_id)),
                'marker_min_px'          : int(marker_min_px),
                'marker_max_dev'         : float(marker_max_dev),
                'object_order'           : str(object_order),
                'camera_height_cm'       : (None if camera_height_cm is None
                                            else float(camera_height_cm))
            },
            'images'        : run_images
        }
        json_path = os.path.join(output_folder, f'phenoshot_params_{stamp}.json')
        try:
            with open(json_path, 'w', encoding='utf-8') as f:
                json.dump(meta, f, indent=2, ensure_ascii=False)
            print(f'Parameters saved: {json_path}')
        except Exception as e:
            print(f'Failed to save parameter log: {e}')

    return df
"
reticulate::py_run_string(py_code)

df_py= reticulate::py$process_images(
  input_folder            = normalizePath(input_folder,  winslash="\\", mustWork=FALSE),
  output_folder           = normalizePath(output_folder, winslash="\\", mustWork=FALSE),
  image_real_cm_W         = as.numeric(image_real_cm[1]),
  image_real_cm_H         = as.numeric(image_real_cm[2]),
  photoroom_api_key       = as.character(photoroom_api_key),
  photoroom_size          = if (is.null(photoroom_size)) NULL else as.character(photoroom_size),
  min_component_area_px   = as.integer(min_component_area_px),
  k_open                  = as.integer(k_open),
  k_close                 = as.integer(k_close),
  object_min_area_cm2     = as.numeric(object_min_area_cm2),
  rel_min_frac_of_largest = as.numeric(rel_min_frac_of_largest),
  max_keep                = as.integer(max_keep),
  label_uses_aabb         = isTRUE(label_uses_aabb),
  contour_color           = as.integer(contour_color),
  alpha_threshold         = as.integer(alpha_threshold),
  fill_opacity            = as.numeric(fill_opacity),
  distinct_colors         = isTRUE(distinct_colors),
  outline_thickness       = as.integer(outline_thickness),
  show_crosshair          = isTRUE(show_crosshair),
  weed_filter             = as.character(weed_filter),
  leaf_min_width_cm       = as.numeric(leaf_min_width_cm),
  debug_weed_mask         = isTRUE(debug_weed_mask),
  outline_mode            = as.character(outline_mode),
  outline_min_hole_cm2    = as.numeric(outline_min_hole_cm2),
  k_close_iter            = as.integer(k_close_iter),
  scale_marker            = as.character(scale_marker),
  marker_size_cm          = as.numeric(marker_size_cm),
  aruco_dict              = as.character(aruco_dict),
  marker_id               = if (is.null(marker_id)) NULL else as.integer(marker_id),
  marker_min_px           = as.integer(marker_min_px),
  marker_max_dev          = as.numeric(marker_max_dev),
  camera_height_cm        = if (is.null(camera_height_cm)) NULL else as.numeric(camera_height_cm),
  canopy_height_cm        = as.numeric(canopy_height_cm),
  object_order            = as.character(object_order),
  params_json             = isTRUE(params_json)
)

out= reticulate::py_to_r(df_py)
if (is.null(out)) out= data.frame()
out
}
