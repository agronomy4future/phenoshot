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
#' length, and width. Supports JPG, PNG, WEBP, and HEIC image formats.
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
#'
#' @return A data frame with one row per detected object containing columns
#'   File Name, Object ID, Image Path, Object Area (cm2), Object Area (px),
#'   Object Perimeter (cm), Object Length (cm), Object Width (cm),
#'   Object percent of Image, Pixel Area (cm2/px), and Num objects in Image.
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
    show_crosshair          = FALSE
) {
  if (length(image_real_cm) == 1) image_real_cm= rep(image_real_cm, 2)
  stopifnot(length(image_real_cm) == 2, all(is.finite(image_real_cm)), all(image_real_cm > 0))

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

def process_images(input_folder, output_folder,
                   image_real_cm_W, image_real_cm_H,
                   photoroom_api_key, photoroom_size,
                   min_component_area_px, k_open, k_close,
                   object_min_area_cm2, rel_min_frac_of_largest, max_keep,
                   label_uses_aabb, contour_color, alpha_threshold, fill_opacity,
                   distinct_colors, outline_thickness, show_crosshair):

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

        mask = (alpha > int(alpha_threshold)).astype(np.uint8) * 255
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN,
                                np.ones(tuple(k_open),  np.uint8), iterations=1)
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE,
                                np.ones(tuple(k_close), np.uint8), iterations=2)

        h, w = mask.shape[:2]
        area_per_pixel_cm2 = (image_real_cm_W / w) * (image_real_cm_H / h)
        sx = image_real_cm_W / w
        sy = image_real_cm_H / h

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

                # (b) thick outline tracing this object's real contour
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
                    'Object Perimeter (cm)': round(perim_cm, 2),
                    'Object Length (cm)'   : round(float(length_cm), 2),
                    'Object Width (cm)'    : round(float(width_cm), 2),
                    'Object % of Image'    : round(100.0 * area_px / (w * h), 2),
                    'Pixel Area (cm2/px)'  : round(area_per_pixel_cm2, 8),
                    'Num objects in Image' : len(contours)
                })
        else:
            cv2.putText(annotated, 'No objects detected', (20, 60),
                        cv2.FONT_HERSHEY_SIMPLEX, font_scale, (0,0,255), 2, cv2.LINE_AA)

        cv2.imwrite(proc_path, annotated, [cv2.IMWRITE_JPEG_QUALITY, 95])
        print(f'  Saved: {proc_path}')

    df = pd.DataFrame(rows)
    csv_path = os.path.join(output_folder, 'image_processed.csv')
    if len(df) > 0:
        try:
            df.to_csv(csv_path, index=False, encoding='utf-8-sig')
            print(f'CSV saved: {csv_path}')
        except Exception as e:
            print(f'Failed to save CSV: {e}')

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
  show_crosshair          = isTRUE(show_crosshair)
)

out= reticulate::py_to_r(df_py)
if (is.null(out)) out= data.frame()
out
}
