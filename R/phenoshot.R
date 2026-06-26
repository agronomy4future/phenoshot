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
#' length, and width. Supports JPG, PNG, and HEIC image formats.
#' To use the background removal feature, an API key from
#' \url{https://www.remove.bg/api} is required. If pre-processed
#' \code{_nobg.png} files are already available in the input folder,
#' the API key can be omitted.
#'
#' @param input_folder Character. Path to folder containing input images
#'   (JPG, PNG, HEIC) or pre-processed `_nobg.png` files.
#'   If `_nobg.png` files are present, API key is not required.
#' @param output_folder Character. Path to folder for saving results.
#'   Created automatically if it does not exist. Output files include
#'   `_nobg.png`, `_nobg_processed.jpg`, and `image_processed.csv`.
#' @param image_real_cm Numeric vector of length 2. Real-world image dimensions
#'   in centimeters as `c(width_cm, height_cm)`. Used to convert pixel
#'   measurements to cm. Default: `c(20, 20)`.
#' @param removebg_api_key Character. API key for remove.bg background removal
#'   service. Required only when `_nobg.png` files are not pre-existing.
#'   Default: `""`.
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
#'   `c(B, G, R)`. Default: `c(0L, 180L, 255L)` (orange).
#' @param alpha_threshold Integer. Alpha channel threshold between 0 and 255
#'   for object detection. Default: `10L`.
#' @param fill_opacity Numeric between 0 and 1. Opacity of semi-transparent
#'   fill applied to detected objects. 0 means no fill. Default: `0.15`.
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
#' # Case 2: Using original images with remove.bg API
#' phenoshot(
#'   input_folder     = "./Input",
#'   output_folder    = "./Output",
#'   image_real_cm    = c(20, 20),
#'   removebg_api_key = "your_api_key_here",
#'   min_component_area_px   = 500L,
#'   object_min_area_cm2     = 0.3,
#'   max_keep                = 50L,
#'   contour_color           = c(0L, 180L, 255L),
#'   alpha_threshold         = 10L,
#'   fill_opacity            = 0.15
#' )
#'
#' ■ Github: https://github.com/agronomy4future/phenoshot
#' - All Rights Reserved © J.K Kim (kimjk@agronomy4future.com)
#' }
phenoshot= function(
    input_folder,
    output_folder,
    image_real_cm           = c(20, 20),
    removebg_api_key        = "",
    min_component_area_px   = 500L,
    k_open                  = c(3L, 3L),
    k_close                 = c(5L, 5L),
    object_min_area_cm2     = 0.3,
    rel_min_frac_of_largest = 0.05,
    max_keep                = 30L,
    label_uses_aabb         = TRUE,
    contour_color           = c(0L, 180L, 255L),
    alpha_threshold         = 10L,
    fill_opacity            = 0.15
) {
  if (length(image_real_cm) == 1) image_real_cm <- rep(image_real_cm, 2)
  stopifnot(length(image_real_cm) == 2, all(is.finite(image_real_cm)), all(image_real_cm > 0))
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

def _discover_unique_images(folder):
    patterns = [
        '*.jpg','*.jpeg','*.png',
        '*.JPG','*.JPEG','*.PNG',
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

    # _nobg.png 있으면 → 그것만 사용 (API/원본 불필요)
    # 없으면 → 원본 이미지 사용
    nobg_files = [p for p in unique_paths
                  if '_nobg' in os.path.splitext(os.path.basename(p))[0]]
    orig_files = [p for p in unique_paths
                  if '_nobg' not in os.path.splitext(os.path.basename(p))[0]]

    if nobg_files:
        print(f'Found {len(nobg_files)} _nobg image(s) → API skipped.')
        return nobg_files, True   # True = nobg 모드
    print(f'Found {len(orig_files)} original image(s).')
    return orig_files, False      # False = 원본 모드

def _open_image(path):
    return Image.open(path).convert('RGB')

def _removebg(image_path, api_key, nobg_save_path):
    ext = os.path.splitext(image_path)[1].lower()
    if ext in ('.heic', '.heif'):
        print(f'  Converting HEIC to JPEG for API...')
        img = _open_image(image_path)
        orig_w, orig_h = img.size
        max_px = 4000
        if max(orig_w, orig_h) > max_px:
            scale = max_px / max(orig_w, orig_h)
            img = img.resize((int(orig_w*scale), int(orig_h*scale)), Image.LANCZOS)
        buf = io.BytesIO()
        img.save(buf, format='JPEG', quality=92)
        buf.seek(0)
        files = {'image_file': ('image.jpg', buf, 'image/jpeg')}
    else:
        with open(image_path, 'rb') as f:
            img_bytes = f.read()
        files = {'image_file': (os.path.basename(image_path), img_bytes)}

    response = requests.post(
        'https://api.remove.bg/v1.0/removebg',
        files=files,
        data={'size': 'auto'},
        headers={'X-Api-Key': api_key}
    )
    if response.status_code != 200:
        raise RuntimeError(f'remove.bg API error {response.status_code}: {response.text}')
    remaining = response.headers.get('X-Free-Calls', '?')
    print(f'  remove.bg OK. Free calls remaining: {remaining}')
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
                   removebg_api_key,
                   min_component_area_px, k_open, k_close,
                   object_min_area_cm2, rel_min_frac_of_largest, max_keep,
                   label_uses_aabb, contour_color, alpha_threshold, fill_opacity):

    if not os.path.isdir(input_folder):
        print(f'Input folder does not exist: {input_folder}')
        return pd.DataFrame()

    os.makedirs(output_folder, exist_ok=True)
    image_paths, is_nobg_mode = _discover_unique_images(input_folder)
    if not image_paths:
        print('No images found.')
        return pd.DataFrame()

    col  = tuple(int(c) for c in contour_color)
    k3   = np.ones((3,3), np.uint8)
    rows = []

    for path in image_paths:
        filename = os.path.basename(path)
        base     = os.path.splitext(filename)[0]

        if is_nobg_mode:
            # _nobg.png 직접 입력 모드
            # C1_nobg.png → stem = C1
            stem      = base.replace('_nobg', '')
            nobg_path = path   # 이 파일 자체가 nobg
        else:
            # 원본 이미지 모드
            stem      = base
            nobg_path = os.path.join(output_folder, stem + '_nobg.png')

        proc_path = os.path.join(output_folder, stem + '_nobg_processed.jpg')
        print(f'Processing: {filename}')

        # ── 1. nobg.png 준비 ─────────────────────────────────────────
        if is_nobg_mode:
            # _nobg.png 직접 사용
            with open(nobg_path, 'rb') as f:
                png_bytes = f.read()
            print(f'  Using _nobg directly: {nobg_path}')
        elif os.path.exists(nobg_path):
            # Output에 기존 nobg 있으면 재사용
            print(f'  Using existing nobg: {nobg_path}')
            with open(nobg_path, 'rb') as f:
                png_bytes = f.read()
        elif removebg_api_key:
            # API로 생성
            try:
                png_bytes = _removebg(path, removebg_api_key, nobg_path)
            except Exception as e:
                print(f'  remove.bg failed: {e} — skipping')
                continue
        else:
            print(f'  No _nobg.png and no API key — skipping')
            continue

        # ── 2. 알파채널 → 마스크 ─────────────────────────────────────
        rgba  = np.array(Image.open(io.BytesIO(png_bytes)).convert('RGBA'))
        alpha = rgba[:, :, 3]
        rgb   = rgba[:, :, :3]

        mask = (alpha > int(alpha_threshold)).astype(np.uint8) * 255
        mask = cv2.morphologyEx(mask, cv2.MORPH_OPEN,
                                np.ones(tuple(k_open),  np.uint8), iterations=1)
        mask = cv2.morphologyEx(mask, cv2.MORPH_CLOSE,
                                np.ones(tuple(k_close), np.uint8), iterations=2)

        eroded   = cv2.erode(mask, k3, iterations=1)
        boundary = mask - eroded

        h, w = mask.shape[:2]
        area_per_pixel_cm2 = (image_real_cm_W / w) * (image_real_cm_H / h)
        sx = image_real_cm_W / w
        sy = image_real_cm_H / h

        font_scale = max(0.6, round(max(w, h) / 1600, 2))
        line_thick = max(1, int(max(w, h) / 1500))

        # ── 3. Contour 탐지 ──────────────────────────────────────────
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

        # ── 4. 배경 합성 ─────────────────────────────────────────────
        comp = np.ones_like(rgb, dtype=np.uint8) * 255
        comp[alpha > 0] = rgb[alpha > 0]
        annotated = cv2.cvtColor(comp, cv2.COLOR_RGB2BGR)
        annotated[boundary > 0] = col

        if contours:
            for idx, cnt in enumerate(contours, start=1):
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

                if fill_opacity > 0:
                    overlay = annotated.copy()
                    overlay[actual_mask > 0] = col
                    cv2.addWeighted(overlay, float(fill_opacity),
                                    annotated, 1.0 - float(fill_opacity),
                                    0, annotated)
                    annotated[boundary > 0] = col

                cv2.line(annotated, (x, cy), (x+bw, cy), col, line_thick)
                cv2.line(annotated, (cx, y), (cx, y+bh), col, line_thick)

                label = f'Obj {idx} ({area_cm2:.2f} cm2)'
                ty = max(int(font_scale*20)+4, y - 6)
                cv2.putText(annotated, label, (x, ty),
                            cv2.FONT_HERSHEY_SIMPLEX, font_scale,
                            (255,255,255), line_thick+3, cv2.LINE_AA)
                cv2.putText(annotated, label, (x, ty),
                            cv2.FONT_HERSHEY_SIMPLEX, font_scale,
                            (30,30,30), line_thick+1, cv2.LINE_AA)

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
    removebg_api_key        = as.character(removebg_api_key),
    min_component_area_px   = as.integer(min_component_area_px),
    k_open                  = as.integer(k_open),
    k_close                 = as.integer(k_close),
    object_min_area_cm2     = as.numeric(object_min_area_cm2),
    rel_min_frac_of_largest = as.numeric(rel_min_frac_of_largest),
    max_keep                = as.integer(max_keep),
    label_uses_aabb         = isTRUE(label_uses_aabb),
    contour_color           = as.integer(contour_color),
    alpha_threshold         = as.integer(alpha_threshold),
    fill_opacity            = as.numeric(fill_opacity)
  )

  out= reticulate::py_to_r(df_py)
  if (is.null(out)) out= data.frame()
  out
}
