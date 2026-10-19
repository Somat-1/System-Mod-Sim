#!/usr/bin/env python3
"""Read an analogue dial indicator from video: needle angle -> micrometres.

The measurement chain is the same for any dial-indicator video:

    frame  ->  background subtraction  ->  needle pixels  ->  bearing [deg]
           ->  unwrap across frames    ->  turns + fraction
           ->  micrometres

Nothing here is specific to one indicator. The geometry (where the dial is,
how big it is, how many micrometres per turn) lives entirely in a calibration
dict, which `calibrate_dial.py` writes and `track_dial.py` consumes.

Bearing convention used throughout: **degrees clockwise from 12 o'clock**, so
0 = straight up, 90 = right, 180 = down, 270 = left. Image coordinates have y
growing downward, hence the `-dy` in `bearing_map`.

Why background subtraction rather than looking for a dark line: the dial face
is covered in printed ticks, numerals and a maker's name that are just as dark
as the needle. The needle is the only thing that *moves*, so a per-pixel
"what does this pixel look like when the needle is elsewhere" image removes
every static mark at once and leaves the needle almost alone.
"""

from __future__ import annotations

import json
import os
from pathlib import Path

import cv2
import numpy as np

# --------------------------------------------------------------------------
# calibration
# --------------------------------------------------------------------------

#: Everything the tracker needs to know about one video. Lengths in pixels,
#: angles in degrees, readings in micrometres unless the key says otherwise.
DEFAULT_CALIBRATION = {
    "video": "",
    "frame_size": [0, 0],
    "fps": 0.0,
    "n_frames": 0,
    "main": {
        "center": [0.0, 0.0],      # needle pivot, (x, y) in pixels
        "r_inner": 0.0,            # annulus the needle is measured in; must clear
        "r_outer": 0.0,            # the hub, the sub-dial and the bezel
        "um_per_rev": 200.0,       # one full turn of the long pointer
        "ticks_per_rev": 100,      # printed minor ticks (2 um each here)
        "zero_bearing_deg": 0.0,   # bearing the needle shows when the dial reads 0
        "direction": "cw",         # which way the needle turns for an increasing reading
    },
    "sub": {
        # The revolution counter. Optional: set "enabled": false to skip it.
        "enabled": True,
        "center": [0.0, 0.0],
        "r_inner": 0.0,
        "r_outer": 0.0,
        "mm_per_step": 0.2,        # one printed step = one full turn of the main needle
        "deg_per_step": -90.0,     # signed: how far the small needle moves per step
        "zero_bearing_deg": 0.0,   # bearing of the printed "0"
    },
    "detect": {
        "bg_frames": 240,          # frames sampled to build the background
        "bg_percentile": 88,       # used only for the sub-dial background
        "main_thresh_frac": 0.30,  # threshold = max(min, frac * brightest needle pixel)
        "main_thresh_min": 12.0,
        "main_tol_deg": 12.0,      # outlier rejection half-width, main needle
        "sub_thresh": 26.0,
        "sub_tol_deg": 20.0,
        "wedge_exclude_deg": 6.0,  # main needle is masked out of the sub-dial by this much
        "min_pixels_main": 40,
        "min_pixels_sub": 20,
    },
    "export": {
        "start_frame": 0,          # first frame written to the CSV
        "end_frame": -1,           # -1 = to the end
        "reference_frame": 0,      # frame whose reading is known
        "reference_reading_um": 0.0,
        "bad_frames": [],          # frames the tracker could not trust
    },
}


def _merge(base, over):
    out = dict(base)
    for k, v in (over or {}).items():
        out[k] = _merge(base[k], v) if isinstance(v, dict) and isinstance(base.get(k), dict) else v
    return out


def load_calibration(path):
    """Read a calibration JSON, filling anything absent from the defaults."""
    with open(path, "r", encoding="utf-8") as fh:
        return _merge(DEFAULT_CALIBRATION, json.load(fh))


def save_calibration(path, cal):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(cal, fh, indent=2)
        fh.write("\n")


# --------------------------------------------------------------------------
# video helpers
# --------------------------------------------------------------------------

def open_video(path):
    """Open a video, raising something readable if it is not there.

    OpenCV happily returns a closed capture for a missing file, and on Windows
    it fails the same way for a UNC path typed with the wrong number of
    backslashes, so check existence separately to tell the two apart.
    """
    path = str(path)
    if not os.path.exists(path):
        raise FileNotFoundError(f"video not found: {path}")
    cap = cv2.VideoCapture(path)
    if not cap.isOpened():
        raise RuntimeError(f"OpenCV could not open (codec not supported?): {path}")
    return cap


def video_info(path):
    cap = open_video(path)
    info = {
        "n_frames": int(cap.get(cv2.CAP_PROP_FRAME_COUNT)),
        "fps": float(cap.get(cv2.CAP_PROP_FPS)),
        "width": int(cap.get(cv2.CAP_PROP_FRAME_WIDTH)),
        "height": int(cap.get(cv2.CAP_PROP_FRAME_HEIGHT)),
    }
    cap.release()
    return info


def read_gray(cap, index=None):
    if index is not None:
        cap.set(cv2.CAP_PROP_POS_FRAMES, int(index))
    ok, frame = cap.read()
    if not ok:
        return None
    return cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def bearing_map(shape, center):
    """Per-pixel bearing in degrees clockwise from 12 o'clock."""
    h, w = shape[:2]
    yy, xx = np.mgrid[0:h, 0:w]
    dx = (xx - center[0]).astype(np.float32)
    dy = (yy - center[1]).astype(np.float32)
    return np.degrees(np.arctan2(dx, -dy)) % 360.0


def radius_map(shape, center):
    h, w = shape[:2]
    yy, xx = np.mgrid[0:h, 0:w]
    return np.hypot(xx - center[0], yy - center[1]).astype(np.float32)


def angdiff(a, b):
    """Signed smallest difference a - b, in (-180, 180]."""
    return (np.asarray(a) - np.asarray(b) + 180.0) % 360.0 - 180.0


def circular_mean(deg, weights):
    a = np.radians(deg)
    return np.degrees(np.arctan2((weights * np.sin(a)).sum(),
                                 (weights * np.cos(a)).sum())) % 360.0


def fit_bearing(deg, weights, tol_deg, min_pixels):
    """Intensity-weighted mean bearing of a set of needle pixels.

    Starts from the circular mean of everything, then re-centres a few times
    while discarding pixels more than `tol_deg` away. That rejection step is
    what stops a stray highlight or a bit of uncancelled print from dragging
    the estimate off the needle.

    Returns (bearing_deg, n_pixels_used, angular_spread_deg). The spread is a
    useful quality flag: it is the needle's apparent half-width, so it jumps
    when the needle is motion-blurred or only partly visible.
    """
    if deg.size < min_pixels:
        return np.nan, 0, np.nan
    th = circular_mean(deg, weights)
    keep = np.ones(deg.shape, bool)
    for _ in range(4):
        dev = np.abs(angdiff(deg, th))
        keep = dev < tol_deg
        if keep.sum() < min_pixels:
            return np.nan, int(keep.sum()), np.nan
        th = circular_mean(deg[keep], weights[keep])
    dev = np.abs(angdiff(deg, th))
    keep = dev < tol_deg
    w = weights[keep]
    spread = float(np.sqrt((w * dev[keep] ** 2).sum() / w.sum())) if w.sum() > 0 else np.nan
    return float(th), int(keep.sum()), spread


# --------------------------------------------------------------------------
# background
# --------------------------------------------------------------------------

def build_background(path, center, bearings=None, n_frames=240, percentile=88,
                     wedge_deg=5.0, mask_disc=None, block=180, progress=None):
    """Per-pixel image of the dial face *without* the needle.

    Two ways to get there:

    * `bearings=None` - take a high percentile over time. The needle is dark,
      so the bright percentile is mostly needle-free. This is the bootstrap
      pass, and it degrades wherever the needle lingers (at the turning points
      of a sweep) because then the needle occupies more than (100-percentile)%
      of the frames at that angle.
    * `bearings` given - we already know where the needle was in each frame, so
      blank out a wedge around it and take the plain median of what is left.
      This is exact wherever the needle was not sitting for the whole clip.

    `mask_disc` is an optional (x, y, r): a region blanked in every frame, used
    to keep the sub-dial's own needle out of the main background. Those pixels
    come back as the global median and must not be relied on - build a separate
    background for the sub-dial.
    """
    cap = open_video(path)
    total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    sel = np.linspace(0, total - 1, min(n_frames, total)).astype(int)
    frames = []
    for i in sel:
        g = read_gray(cap, i)
        if g is not None:
            frames.append(g)
    cap.release()
    stack = np.asarray(frames, dtype=np.float32)
    h, w = stack.shape[1:]
    amap = bearing_map((h, w), center)
    disc = None
    if mask_disc is not None:
        rr = radius_map((h, w), mask_disc[:2])
        disc = rr < mask_disc[2]
    bg = np.zeros((h, w), np.float32)
    for y0 in range(0, h, block):
        y1 = min(y0 + block, h)
        blk = stack[:, y0:y1, :].copy()
        sub_amap = amap[y0:y1]
        sub_disc = None if disc is None else disc[y0:y1]
        for j, i in enumerate(sel):
            m = None
            if bearings is not None and np.isfinite(bearings[i]):
                m = np.abs(angdiff(sub_amap, bearings[i])) < wedge_deg
            if sub_disc is not None:
                m = sub_disc if m is None else (m | sub_disc)
            if m is not None:
                blk[j][m] = np.nan
        with np.errstate(invalid="ignore"):
            if bearings is None:
                bg[y0:y1] = np.nanpercentile(blk, percentile, axis=0)
            else:
                bg[y0:y1] = np.nanmedian(blk, axis=0)
        if progress:
            progress(y1, h)
    fill = np.nanmedian(bg)
    return np.nan_to_num(bg, nan=fill).astype(np.uint8)


# --------------------------------------------------------------------------
# centre fitting
# --------------------------------------------------------------------------

def fit_rotation_center(path, bg, approx_center, r_inner, r_outer, n_frames=80,
                        thresh=30.0, exclude_disc=None, iterations=4):
    """Locate the needle pivot without being told where it is.

    Every frame gives a straight line through the needle; all those lines meet
    at the pivot. Fit a line per frame by weighted PCA, then solve for the
    point closest to all of them, re-weighting (Cauchy) so that a handful of
    bad frames cannot move the answer.
    """
    cap = open_video(path)
    total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    sel = np.linspace(0, total - 1, min(n_frames, total)).astype(int)
    diffs = []
    for i in sel:
        g = read_gray(cap, i)
        if g is not None:
            diffs.append(cv2.subtract(bg, g))
    cap.release()
    h, w = bg.shape
    yy, xx = np.mgrid[0:h, 0:w]
    center = np.asarray(approx_center, float)
    excl = None
    if exclude_disc is not None:
        excl = np.hypot(xx - exclude_disc[0], yy - exclude_disc[1]) < exclude_disc[2]
    for _ in range(iterations):
        rr = np.hypot(xx - center[0], yy - center[1])
        ann = (rr > r_inner) & (rr < r_outer)
        if excl is not None:
            ann &= ~excl
        normals, offsets = [], []
        for d in diffs:
            m = ann & (d > thresh)
            if m.sum() < 120:
                continue
            pts = np.c_[xx[m], yy[m]].astype(float)
            wgt = d[m].astype(float)
            mu = (pts * wgt[:, None]).sum(0) / wgt.sum()
            X = (pts - mu) * np.sqrt(wgt)[:, None]
            _, s, vt = np.linalg.svd(X, full_matrices=False)
            if s[0] / max(s[1], 1e-9) < 8:      # not line-like enough to trust
                continue
            n = np.array([-vt[0][1], vt[0][0]])
            normals.append(n)
            offsets.append(n @ mu)
        if len(normals) < 5:
            break
        A = np.asarray(normals)
        b = np.asarray(offsets)
        wt = np.ones(len(b))
        for _ in range(15):
            sol, *_ = np.linalg.lstsq(A * wt[:, None], b * wt, rcond=None)
            res = np.abs(A @ sol - b)
            scale = max(np.median(res), 0.3)
            wt = 1.0 / np.sqrt(1.0 + (res / (2.5 * scale)) ** 2)
        center = sol
    residual = float(np.median(np.abs(A @ center - b))) if len(normals) >= 5 else float("nan")
    return center, residual


# --------------------------------------------------------------------------
# tracking
# --------------------------------------------------------------------------

def track(path, cal, bg_main, bg_sub=None, progress=None):
    """Measure both needles in every frame.

    Returns a dict of per-frame arrays, all length n_frames:
        main_bearing, main_pixels, main_spread, main_contrast
        sub_bearing,  sub_pixels
    with NaN wherever a needle could not be found.
    """
    det = cal["detect"]
    mainc = cal["main"]
    subc = cal["sub"]
    cap = open_video(path)
    n = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    shape = bg_main.shape
    rm = radius_map(shape, mainc["center"])
    am = bearing_map(shape, mainc["center"])
    ann_m = (rm > mainc["r_inner"]) & (rm < mainc["r_outer"])
    am_in = am[ann_m]
    use_sub = bool(subc.get("enabled", True)) and bg_sub is not None
    if use_sub:
        rs = radius_map(shape, subc["center"])
        asb = bearing_map(shape, subc["center"])
        ann_s = (rs > subc["r_inner"]) & (rs < subc["r_outer"])

    out = {k: np.full(n, np.nan) for k in
           ("main_bearing", "main_spread", "sub_bearing", "main_contrast")}
    out["main_pixels"] = np.zeros(n, int)
    out["sub_pixels"] = np.zeros(n, int)

    cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
    for i in range(n):
        ok, frame = cap.read()
        if not ok:
            break
        g = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)

        d = cv2.subtract(bg_main, g)
        dv = d[ann_m]
        contrast = float(dv.max())
        out["main_contrast"][i] = contrast
        thr = max(det["main_thresh_min"], det["main_thresh_frac"] * contrast)
        sel = dv > thr
        if sel.sum() >= det["min_pixels_main"]:
            th, npx, sp = fit_bearing(am_in[sel], dv[sel].astype(float),
                                      det["main_tol_deg"], det["min_pixels_main"])
            out["main_bearing"][i] = th
            out["main_pixels"][i] = npx
            out["main_spread"][i] = sp

        if use_sub:
            ds = cv2.subtract(bg_sub, g)
            keep = ann_s & (ds > det["sub_thresh"])
            if np.isfinite(out["main_bearing"][i]):
                # the long needle sweeps across the sub-dial; mask it out
                keep &= np.abs(angdiff(am, out["main_bearing"][i])) > det["wedge_exclude_deg"]
            if keep.sum() >= det["min_pixels_sub"]:
                th, npx, _ = fit_bearing(asb[keep], ds[keep].astype(float),
                                         det["sub_tol_deg"], det["min_pixels_sub"])
                out["sub_bearing"][i] = th
                out["sub_pixels"][i] = npx

        if progress and i % 100 == 0:
            progress(i, n)
    cap.release()
    return out


# --------------------------------------------------------------------------
# quality / disturbance detection
# --------------------------------------------------------------------------

def flag_bad_frames(tr, spread_factor=3.0, contrast_factor=0.55):
    """Frames whose needle estimate should not be believed.

    A frame is bad if the needle could not be found at all, or if it looks
    smeared rather than sharp. Both happen for the same physical reason: the
    needle moved far within one exposure, so instead of a crisp dark line the
    frame holds a faint smear spread over many degrees. That shows up as a
    large angular spread and a collapsed peak contrast.
    """
    spread = tr["main_spread"]
    contrast = tr["main_contrast"]
    ref_spread = np.nanmedian(spread)
    ref_contrast = np.nanmedian(contrast)
    bad = ~np.isfinite(tr["main_bearing"])
    bad |= ~np.isfinite(spread)
    bad |= spread > spread_factor * ref_spread
    bad |= contrast < contrast_factor * ref_contrast
    return bad


def first_clean_run(bad, min_run=5, after=0):
    """Index of the first frame starting a run of `min_run` good frames.

    Used to answer "when did the needle settle back down after the startup
    jolt" - scan past the disturbance and take the first frame from which the
    tracker stays happy for a while.
    """
    n = len(bad)
    run = 0
    for i in range(after, n):
        run = 0 if bad[i] else run + 1
        if run >= min_run:
            return i - min_run + 1
    return after


def find_disturbance(bad, head_fraction=0.25, gap=8):
    """Locate the first burst of bad frames near the start of the clip.

    Returns (first_bad, last_bad, resume_frame) or None. Only the first part of
    the clip is searched, because this is meant for a power-on / initialisation
    transient, not for an arbitrary dropout halfway through.

    Bad frames are grouped into bursts and only the *first* burst is reported.
    Without that, a single marginal frame seconds later would be absorbed into
    the event and throw away all the good frames in between.
    """
    n = len(bad)
    head = max(1, int(n * head_fraction))
    idx = np.flatnonzero(bad[:head])
    if idx.size == 0:
        return None
    burst = np.split(idx, np.flatnonzero(np.diff(idx) > gap) + 1)[0]
    first, last = int(burst[0]), int(burst[-1])
    return first, last, first_clean_run(bad, min_run=5, after=last + 1)


# --------------------------------------------------------------------------
# angles -> readings
# --------------------------------------------------------------------------

def flag_sub_occlusion(sub_pixels, fraction=0.45):
    """Frames where too little of the small needle survived to trust its angle.

    Cause: the long pointer sweeps over the revolution counter twice per turn -
    in this video it passes within 3 px of the counter's centre - and hides the
    small needle underneath it. The tracker masks the long needle out, so what
    is left is a partly eaten blob whose fitted angle is pulled off true by
    tens of degrees.

    Rather than predict that from geometry, measure it: count how many pixels
    actually went into the fit and reject the frame when that collapses. This
    catches the same frames as the geometric test on a quarter as many false
    positives, because it responds to the thing that actually went wrong.

    Only the cross-check is affected - the main reading does not use the
    counter at all.
    """
    sp = np.asarray(sub_pixels, float)
    valid = sp > 0
    if valid.sum() < 10:
        return np.zeros(sp.shape, bool)
    return valid & (sp < fraction * np.median(sp[valid]))


def unwrap_bearings(bearing, bad=None):
    """Turn a wrapped 0-360 sequence into a continuous angle.

    Frames marked bad are interpolated across *before* unwrapping, so a
    garbage reading in the middle of a fast move cannot inject a spurious
    full turn. This is only safe while the needle moves less than half a turn
    between two good frames - `max_step_deg` in the report lets you check that.
    """
    b = np.asarray(bearing, float).copy()
    if bad is not None:
        b[bad] = np.nan
    idx = np.arange(b.size)
    ok = np.isfinite(b)
    if ok.sum() < 2:
        raise ValueError("not enough valid bearings to unwrap")
    # interpolate in the complex plane so the gap-filling respects wrap-around
    z = np.exp(1j * np.radians(b[ok]))
    zi = np.interp(idx, idx[ok], z.real) + 1j * np.interp(idx, idx[ok], z.imag)
    filled = np.degrees(np.angle(zi)) % 360.0
    return np.degrees(np.unwrap(np.radians(filled)))


def to_microns(unwrapped, cal):
    """Continuous angle -> dial reading in micrometres."""
    m = cal["main"]
    sign = 1.0 if m["direction"] == "cw" else -1.0
    ref_i = int(cal["export"]["reference_frame"])
    ref_um = float(cal["export"]["reference_reading_um"])
    scale = m["um_per_rev"] / 360.0
    return ref_um + sign * (unwrapped - unwrapped[ref_i]) * scale


def main_face_microns(bearing, cal):
    """What the long pointer alone reads, 0..um_per_rev, ignoring which turn."""
    m = cal["main"]
    sign = 1.0 if m["direction"] == "cw" else -1.0
    rel = (sign * (np.asarray(bearing, float) - m["zero_bearing_deg"])) % 360.0
    return rel * m["um_per_rev"] / 360.0


def sub_millimetres(sub_bearing, cal):
    """Independent coarse reading from the revolution counter.

    The counter's scale wraps once per *indicator range*, not once per turn of
    the small needle, so the travel is taken as the one-way sweep from the
    printed zero rather than the shortest angular difference. On a dial whose
    marks span 270 deg that leaves the last quarter turn as a dead zone, which
    is where an out-of-range reading will land - obvious rather than silently
    aliased onto a valid number.
    """
    s = cal["sub"]
    if not s.get("enabled", True) or not s.get("deg_per_step"):
        return np.full(np.shape(sub_bearing), np.nan)
    step = s["deg_per_step"]
    sweep = (np.sign(step) * (np.asarray(sub_bearing, float) - s["zero_bearing_deg"])) % 360.0
    return sweep / abs(step) * s["mm_per_step"]
