#!/usr/bin/env python3
"""Render an explainer video: the dial being tracked, next to what it measures.

Three panels, 1920x1080 at the source frame rate:

    +----------------------+------------------------+
    |                      |  shaft outline, polar  |
    |   the dial, with     |  ideal circle + trace  |
    |   the tracking drawn +------------------------+
    |   on top             |  linear sweep:         |
    |                      |  displacement vs angle |
    +----------------------+------------------------+

Left shows what the tracker actually sees: the two fitted pivots, the annuli
the needles are measured in, the pixels picked out as needle in this frame, and
the fitted bearings. Frames rejected as untrustworthy are called out rather
than hidden.

Right shows what comes out, both panels from the dial alone and both referenced
to the reading on the first tracked frame. The polar chart magnifies the
deviation about the ideal circle (default 15x) because an eccentricity of a few
hundred micrometres on a shaft is otherwise invisible; the linear sweep needs no
gain, since its y-axis is micrometres directly.

The shaft encoder is deliberately absent. It and the dial were not sampled off a
common clock, so putting them on one time axis would imply a synchronisation
that does not exist - the encoder has its own plot in
analysis/plot_encoder_full_rotation_xy.py.

Usage:
    python make_explainer_video.py --calib calibration_video_project_5.json
    python make_explainer_video.py --calib cal.json --gain 25 --fps 30
    python make_explainer_video.py --calib cal.json --max-frames 60   # quick look
"""

import argparse
import csv
import sys
from pathlib import Path

import cv2
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dial_core as dc  # noqa: E402

PANEL_W, PANEL_H = 840, 540
VIDEO_W = 1080
INK = {
    "ideal": "#333333",
    "trace": "#0072B2",
    "live": "#D55E00",
    "grid": "#D7D7D7",
}
BGR = {
    "centre": (0, 255, 255),
    "annulus": (235, 180, 0),
    "needle": (40, 40, 255),
    "sub": (255, 0, 255),
    "zero": (0, 230, 0),
    "mask": (60, 220, 60),
    "bad": (0, 0, 255),
    "text": (255, 255, 255),
}


# ------------------------------------------------------------ plot panels
def _figure():
    fig = plt.figure(figsize=(PANEL_W / 100.0, PANEL_H / 100.0), dpi=100)
    return fig


def _render(fig):
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())[:, :, :3]
    return cv2.cvtColor(buf.copy(), cv2.COLOR_RGB2BGR)


def _to_px(ax, x, y):
    """Data coordinates -> pixel coordinates in the rendered panel."""
    pts = ax.transData.transform(np.column_stack([np.ravel(x), np.ravel(y)]))
    h = ax.figure.canvas.get_width_height()[1]
    pts[:, 1] = h - pts[:, 1]
    return pts.astype(np.int32)


def build_polar_panel(theta_deg, dev_um, gain, r0_um, ecc, phase, tir):
    """Static backdrop for the shaft-outline chart, plus a data->pixel mapper.

    Radius is r0 + gain * deviation, so a perfectly round, perfectly centred
    shaft is the dashed circle. r0 is a drawing reference only - the real shaft
    radius is not known from this measurement and does not matter, because what
    is being shown is the deviation.
    """
    fig = _figure()
    ax = fig.add_axes([0.02, 0.06, 0.96, 0.84])
    ax.set_aspect("equal")
    lim = r0_um + gain * max(np.abs(dev_um).max(), 1.0) * 1.5
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)

    th = np.linspace(0, 2 * np.pi, 721)
    # graduation circles, labelled in true micrometres of deviation. The labels
    # sit on a quiet diagonal: straight up collides with the title and legend.
    # 9 o'clock: the panel is wider than it is tall, so that side is dead space
    # whatever the trace does, and the labels never land on the data
    lab = np.radians(264.0)
    # Grid every 50 µm, but only label every 100 µm and only out as far as the
    # data actually goes - labelled rings all the way to the frame edge collide
    # with each other and with the angle labels.
    step = 50.0
    reach = float(np.abs(dev_um).max())
    k = 1
    while k * step <= reach + step:
        for sign in (+1, -1):
            r = r0_um + sign * k * step * gain
            if r <= 0:
                continue
            ax.plot(r * np.sin(th), r * np.cos(th), color=INK["grid"], lw=0.7, zorder=1)
            if (k * step) % 100 == 0:
                ax.text(r * np.sin(lab), r * np.cos(lab), f"{sign*k*step:+.0f}",
                        fontsize=7, color="#999", ha="center", va="center", zorder=4,
                        bbox=dict(fc="white", ec="none", pad=0.6, alpha=.85))
        k += 1
    for a in range(0, 360, 30):                      # radial spokes
        r = np.radians(a)
        ax.plot([0, lim * np.sin(r)], [0, lim * np.cos(r)], color=INK["grid"],
                lw=0.6, zorder=1)
        ax.text(lim * 0.97 * np.sin(r), lim * 0.97 * np.cos(r), f"{a}°",
                fontsize=7.5, color="#888", ha="center", va="center", zorder=4)
    ax.plot(r0_um * np.sin(th), r0_um * np.cos(th), "--", color=INK["ideal"],
            lw=1.6, zorder=3, label="ideal: perfect shaft, no displacement")
    ax.plot([], [], color=INK["trace"], lw=1.6, label=f"measured displacement (×{gain:g})")
    # The trace is referenced to the reading at 0 deg, so it leaves the ideal
    # circle from this point - mark it, or the reader has to infer it.
    ax.plot([0], [r0_um], "o", ms=6, mfc="white", mec=INK["ideal"], mew=1.4, zorder=5)
    ax.annotate("start / zero", xy=(0, r0_um), xytext=(0, r0_um * 0.74),
                fontsize=7.5, color="#555", ha="center", va="top", zorder=5,
                arrowprops=dict(arrowstyle="-", color="#999", lw=0.8))
    # A pure eccentricity displaces the whole traced circle; showing both
    # centres makes the distinction between "off-centre" and "out-of-round"
    # visible, which is the whole point of the chart.
    ax.plot([0], [0], "+", color="#000", ms=13, mew=2.0, zorder=6,
            label="true rotation centre")
    ex, ey = gain * ecc * np.sin(np.radians(phase)), gain * ecc * np.cos(np.radians(phase))
    ax.plot([ex], [ey], "o", color="#D62728", ms=7, zorder=6,
            label="centre of the measured circle")
    ax.annotate("", xy=(ex, ey), xytext=(0, 0), zorder=5,
                arrowprops=dict(arrowstyle="->", color="#D62728", lw=1.3))
    ax.legend(loc="upper left", fontsize=8, framealpha=.92, ncol=1,
              bbox_to_anchor=(0.0, 1.02))
    ax.text(0.995, 0.015,
            f"eccentricity  {ecc:.1f} µm  at {phase:.0f}°\nTIR  {tir:.1f} µm",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=8.5,
            bbox=dict(boxstyle="round,pad=0.4", fc="white", ec="#AAA", alpha=.92))
    ax.set_title("Shaft outline over one revolution", fontsize=11.5, pad=8)
    ax.text(0.0, -0.02, f"radial displacement from the reading at 0°, µm, magnified ×{gain:g}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#666")

    r = r0_um + gain * dev_um
    xs, ys = r * np.sin(np.radians(theta_deg)), r * np.cos(np.radians(theta_deg))
    base = _render(fig)
    px = _to_px(ax, xs, ys)
    centre_px = _to_px(ax, [0], [0])[0]
    plt.close(fig)
    return base, px, centre_px


def runout_terms(theta_deg, dev_um):
    """Eccentricity and total indicated runout from the measured sweep.

    A shaft mounted off-centre by e reads as a single sinusoid of amplitude e,
    so the once-per-revolution term *is* the eccentricity and its phase is the
    direction of the high spot. TIR is simply the full swing of the indicator.
    Both are independent of where the sweep was zeroed.
    """
    th = np.radians(theta_deg)
    M = np.c_[np.ones_like(th), np.cos(th), np.sin(th)]
    c, *_ = np.linalg.lstsq(M, dev_um, rcond=None)
    ecc = float(np.hypot(c[1], c[2]))
    phase = float(np.degrees(np.arctan2(c[2], c[1])) % 360.0)
    tir = float(dev_um.max() - dev_um.min())
    return ecc, phase, tir


def build_sweep_panel(theta_deg, dev_um):
    """Indicator displacement against shaft angle - the sweep, unmagnified.

    No gain here: the y-axis is micrometres directly, so the sinusoid is the
    measurement at true scale. Returns the backdrop, the data->pixel mapping,
    and the markers to reveal as the trace reaches them.
    """
    fig = _figure()
    ax = fig.add_axes([0.095, 0.155, 0.885, 0.73])
    xmax = max(360.0, float(theta_deg.max()))
    pad = 0.08 * (dev_um.max() - dev_um.min())
    ax.axhline(0, ls="--", color=INK["ideal"], lw=1.4, zorder=2,
               label="zero reference (reading at 0°)")
    ax.plot([], [], color=INK["trace"], lw=1.8, label="measured displacement")
    ax.set(xlim=(0, xmax), ylim=(dev_um.min() - pad, dev_um.max() + pad),
           xlabel="shaft rotation angle [deg]",
           ylabel="indicator displacement [µm]")
    ax.set_xticks(np.arange(0, 361, 45))
    ax.grid(True, color=INK["grid"], lw=0.7, alpha=0.8)
    ax.set_title("Dial indicator linear sweep", fontsize=11.5, pad=10)
    ax.legend(loc="upper right", fontsize=8, framealpha=.92)

    hi, lo = int(np.argmax(dev_um)), int(np.argmin(dev_um))
    marks = [(hi, f"max high spot  {dev_um[hi]:+.1f} um", -1),
             (lo, f"max low spot  {dev_um[lo]:+.1f} um", +1)]
    # where the sweep comes back through its own zero
    sign = np.sign(dev_um)
    for j in np.flatnonzero((sign[:-1] * sign[1:] < 0)):
        marks.append((int(j), f"{theta_deg[j]:.0f} deg crossing", +1))

    base = _render(fig)
    px = _to_px(ax, theta_deg, dev_um)
    plt.close(fig)
    return base, px, [(j, text, side, tuple(px[j])) for j, text, side in marks]


def draw_trace(base, px, upto, colour, live_colour, width=2):
    """Copy the backdrop and draw the trace up to `upto`, with a live dot."""
    img = base.copy()
    if upto >= 2:
        cv2.polylines(img, [px[:upto]], False, colour, width, cv2.LINE_AA)
    if upto >= 1:
        cv2.circle(img, tuple(px[upto - 1]), 6, live_colour, -1, cv2.LINE_AA)
        cv2.circle(img, tuple(px[upto - 1]), 6, (255, 255, 255), 1, cv2.LINE_AA)
    return img


# -------------------------------------------------------------- left panel
def draw_dial(frame, bg, cal, bearing, sub_bearing, info, mask_alpha=0.75):
    """The video frame with everything the tracker used drawn on it."""
    m, s, det = cal["main"], cal["sub"], cal["detect"]
    img = frame.copy()
    C = np.asarray(m["center"], float)
    S = np.asarray(s["center"], float)

    # --- the pixels this frame's bearing was actually fitted from ----------
    gray = cv2.cvtColor(frame, cv2.COLOR_BGR2GRAY)
    d = cv2.subtract(bg, gray)
    rm = dc.radius_map(gray.shape, C)
    ann = (rm > m["r_inner"]) & (rm < m["r_outer"])
    thr = max(det["main_thresh_min"], det["main_thresh_frac"] * float(d[ann].max()))
    sel = ann & (d > thr)
    rs = dc.radius_map(gray.shape, S)
    sub_sel = (rs > s["r_inner"]) & (rs < s["r_outer"]) & (d > det["sub_thresh"])
    tint = np.zeros_like(img)
    tint[sel] = BGR["mask"]
    tint[sub_sel] = BGR["sub"]
    hit = sel | sub_sel
    img[hit] = cv2.addWeighted(img, 1 - mask_alpha, tint, mask_alpha, 0)[hit]

    # --- the geometry that was calibrated once ----------------------------
    for r in (m["r_inner"], m["r_outer"]):
        cv2.circle(img, tuple(C.astype(int)), int(r), BGR["annulus"], 1, cv2.LINE_AA)
    for r in (s["r_inner"], s["r_outer"]):
        cv2.circle(img, tuple(S.astype(int)), int(r), BGR["sub"], 1, cv2.LINE_AA)
    _ray(img, C, m["zero_bearing_deg"], m["r_inner"], m["r_outer"], BGR["zero"], 1)
    for c, size in ((C, 13), (S, 9)):
        p = c.astype(int)
        cv2.line(img, (p[0] - size, p[1]), (p[0] + size, p[1]), BGR["centre"], 2, cv2.LINE_AA)
        cv2.line(img, (p[0], p[1] - size), (p[0], p[1] + size), BGR["centre"], 2, cv2.LINE_AA)

    # --- what was measured in this frame ----------------------------------
    # The fitted bearing is drawn hairline where it crosses the measuring ring
    # so the highlighted pixels underneath stay visible - those pixels are the
    # evidence, the line is only the conclusion - and boldly outside it, where
    # it has nothing to hide.
    if np.isfinite(bearing):
        _ray(img, C, bearing, m["r_inner"] * 0.25, m["r_outer"], BGR["needle"], 1)
        _ray(img, C, bearing, m["r_outer"], m["r_outer"] + 46, BGR["needle"], 3)
    if np.isfinite(sub_bearing):
        _ray(img, S, sub_bearing, s["r_inner"], s["r_outer"], BGR["sub"], 1)
        _ray(img, S, sub_bearing, s["r_outer"], s["r_outer"] + 26, BGR["sub"], 3)

    _legend(img, info)
    _key(img)
    return img


def _key(img):
    """Colour key, bottom left, so the overlay reads without the README."""
    items = [("yellow cross", BGR["centre"], "fitted pivot (calibration)"),
             ("ring", BGR["annulus"], "measuring annulus"),
             ("green", BGR["mask"], "pixels used for this frame's fit"),
             ("red", BGR["needle"], "fitted needle bearing"),
             ("magenta", BGR["sub"], "revolution counter"),
             ("green ray", BGR["zero"], "dial zero")]
    h = img.shape[0]
    y = h - 18 - 23 * len(items)
    _plate(img, 0, y - 20, 360, h)
    for _, colour, label in items:
        cv2.rectangle(img, (18, y - 9), (38, y + 3), colour, -1)
        _text(img, 46, y, label, BGR["text"], 0.46, 1)
        y += 23


def _ray(img, centre, bearing_deg, r0, r1, colour, w):
    b = np.radians(bearing_deg)
    p0 = (int(centre[0] + r0 * np.sin(b)), int(centre[1] - r0 * np.cos(b)))
    p1 = (int(centre[0] + r1 * np.sin(b)), int(centre[1] - r1 * np.cos(b)))
    cv2.line(img, p0, p1, colour, w, cv2.LINE_AA)


def _text(img, x, y, s, colour=BGR["text"], scale=0.62, w=1):
    cv2.putText(img, s, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, (0, 0, 0), w + 3, cv2.LINE_AA)
    cv2.putText(img, s, (x, y), cv2.FONT_HERSHEY_SIMPLEX, scale, colour, w, cv2.LINE_AA)


def _plate(img, x0, y0, x1, y1, alpha=0.55):
    """Darken a rectangle so overlay text stays legible over the bright face."""
    roi = img[y0:y1, x0:x1]
    img[y0:y1, x0:x1] = cv2.addWeighted(roi, 1 - alpha, np.zeros_like(roi), alpha, 0)


def _legend(img, info):
    _plate(img, 0, 0, 480, 48 + 27 * len(info["lines"]))
    y = 34
    _text(img, 18, y, info["title"], scale=0.74, w=2); y += 34
    for line, colour in info["lines"]:
        _text(img, 18, y, line, colour); y += 27
    if info.get("flagged"):
        h, w = img.shape[:2]
        cv2.rectangle(img, (6, 6), (w - 6, h - 6), BGR["bad"], 5)
        _text(img, 18, 196, "FRAME REJECTED", BGR["bad"], 0.86, 2)
        _text(img, 18, 224, "needle smeared by fast motion", BGR["bad"], 0.58, 1)


# -------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calib", required=True)
    ap.add_argument("--out", default="../recordings/dial_explainer.mp4",
                    help="default lands beside the source video, which git ignores")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="resize the finished canvas, e.g. 0.667 for 720p")
    ap.add_argument("--gain", type=float, default=15.0,
                    help="magnification of the deviation in the polar chart only")
    ap.add_argument("--r0-um", type=float, default=5000.0,
                    help="drawing radius of the ideal circle [µm]")
    ap.add_argument("--fps", type=float, default=None)
    ap.add_argument("--max-frames", type=int, default=None, help="render only the first N")
    ap.add_argument("--step", type=int, default=1, help="render every Nth frame")
    args = ap.parse_args()

    here = Path(__file__).resolve().parent
    cal_path = Path(args.calib)
    if not cal_path.is_absolute():
        cal_path = here / cal_path
    cal = dc.load_calibration(cal_path)
    video = cal["video"]
    if not Path(video).is_absolute() and not Path(video).exists():
        video = str((cal_path.parent / video).resolve())
    info = dc.video_info(video)
    fps = args.fps or info["fps"]
    print(f"video   : {video}  ({info['n_frames']} frames @ {info['fps']:.3f} fps)")

    cache = cal_path.parent / "background.npz"
    if not cache.exists():
        raise SystemExit("background.npz missing — run track_dial.py once first")
    bg = np.load(cache)["bg_main"]

    print("tracking the needles …")
    tr = dc.track(video, cal, bg, bg)
    bad = dc.flag_bad_frames(tr)
    unwrapped = dc.unwrap_bearings(tr["main_bearing"], bad)
    um = dc.to_microns(unwrapped, cal)
    n = len(um)

    # --- shaft angle: the dial's own period, fitted from the dial trace ----
    export_start = int(cal["export"]["start_frame"])
    t = np.arange(n) / info["fps"]
    good = slice(export_start, n)
    period = _fit_period(t[good], um[good])
    theta = (t - t[export_start]) / period * 360.0
    # The indicator is wound in before the run so the stylus has room to travel
    # in both directions, which makes the absolute reading an artefact of how
    # far it was preloaded. Zero on the first tracked frame instead: the
    # displacement from there is the measurement, and it puts the start of the
    # polar trace on the ideal circle rather than offset from it by the mean.
    zero_um = float(um[export_start])
    dev = um - zero_um
    print(f"  one revolution = {period:.3f} s  (shaft angle from the dial's own period)")
    print(f"  zero reference = {zero_um:.2f} um at frame {export_start} (preload, not a measurement)")
    print(f"  displacement from zero: {dev[good].min():+.1f} .. {dev[good].max():+.1f} um")

    ecc, phase, tir = runout_terms(theta[good], dev[good])
    print(f"  eccentricity {ecc:.1f} um at {phase:.0f} deg,  TIR {tir:.1f} um")

    # --- static backdrops --------------------------------------------------
    polar_base, polar_px, _ = build_polar_panel(theta[good], dev[good], args.gain,
                                                args.r0_um, ecc, phase, tir)
    sweep_base, sweep_px, sweep_marks = build_sweep_panel(theta[good], dev[good])

    out = Path(args.out)
    if not out.is_absolute():
        out = here / out
    out.parent.mkdir(parents=True, exist_ok=True)
    out_w = int(round((VIDEO_W + PANEL_W) * args.scale)) // 2 * 2
    out_h = int(round(2 * PANEL_H * args.scale)) // 2 * 2
    # Only MPEG-4 Part 2 is available here; OpenCV's H.264 needs the OpenH264
    # DLL, which is not installed. Use --scale to trade size for resolution.
    writer = cv2.VideoWriter(str(out), cv2.VideoWriter_fourcc(*"mp4v"), fps,
                             (out_w, out_h))
    cap = dc.open_video(video)
    frames = range(0, n if args.max_frames is None else min(n, args.max_frames), args.step)

    print(f"rendering {len(list(frames))} frames -> {out.name}")
    cap.set(cv2.CAP_PROP_POS_FRAMES, 0)
    written = 0
    for i in range(n):
        ok, frame = cap.read()
        if not ok:
            break
        if args.max_frames is not None and i >= args.max_frames:
            break
        if i % args.step:
            continue
        started = i >= export_start
        lines = [
            (f"frame {i:4d}   t = {i/info['fps']:6.2f} s", BGR["text"]),
            (f"needle bearing   {tr['main_bearing'][i]:8.3f} deg", BGR["needle"]),
            (f"displacement   {dev[i]:+9.2f} um", (0, 255, 255)),
            (f"shaft angle      {max(theta[i],0):8.2f} deg" if started
             else "shaft angle         -- (before start)", BGR["zero"]),
            (f"zero ref {zero_um:.1f} um @ frame {export_start} (preload)", (185, 185, 185)),
        ]
        left = draw_dial(frame, bg, cal, tr["main_bearing"][i], tr["sub_bearing"][i],
                         {"title": "What the tracker sees", "lines": lines,
                          "flagged": bool(bad[i])})
        k = min(max(0, i - export_start + 1) if started else 0, len(polar_px))
        polar = draw_trace(polar_base, polar_px, k,
                           _bgr(INK["trace"]), _bgr(INK["live"]))
        sweep = draw_trace(sweep_base, sweep_px, k,
                           _bgr(INK["trace"]), _bgr(INK["live"]), width=2)
        _reveal(sweep, sweep_marks, k)
        canvas = np.zeros((2 * PANEL_H, VIDEO_W + PANEL_W, 3), np.uint8)
        canvas[:, :VIDEO_W] = left
        canvas[0:PANEL_H, VIDEO_W:] = polar
        canvas[PANEL_H:, VIDEO_W:] = sweep
        if (out_w, out_h) != canvas.shape[1::-1]:
            canvas = cv2.resize(canvas, (out_w, out_h), interpolation=cv2.INTER_AREA)
        writer.write(canvas)
        written += 1
        if written % 200 == 0:
            print(f"    {written} frames")
    cap.release(); writer.release()
    size_mb = out.stat().st_size / 1e6
    print(f"wrote {out}  ({written} frames, {size_mb:.1f} MB, {written/fps:.1f} s)")


def _fit_period(t, um):
    """Period of the dominant once-per-revolution component."""
    best = None
    for T in np.arange(t[-1] * 0.6, t[-1] * 1.6, 0.01):
        w = 2 * np.pi / T
        M = np.c_[np.ones_like(t), np.cos(w * t), np.sin(w * t)]
        c, *_ = np.linalg.lstsq(M, um, rcond=None)
        r = float((um - M @ c).std())
        if best is None or r < best[0]:
            best = (r, T)
    return best[1]


def _bgr(hex_colour):
    h = hex_colour.lstrip("#")
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    return (b, g, r)


def _reveal(img, marks, upto):
    """Draw the sweep's annotated points, but only once the trace reaches them.

    Putting them in the static backdrop would give the answer away before the
    measurement gets there, which defeats the point of animating it.
    """
    for j, text, side, (x, y) in marks:
        if j >= upto:
            continue
        cv2.circle(img, (x, y), 5, (40, 40, 220), -1, cv2.LINE_AA)
        cv2.circle(img, (x, y), 5, (255, 255, 255), 1, cv2.LINE_AA)
        (tw, th), _ = cv2.getTextSize(text, cv2.FONT_HERSHEY_SIMPLEX, 0.42, 1)
        tx = min(max(x - tw // 2, 4), img.shape[1] - tw - 4)
        # flip the label to the inside when the point sits near a panel edge,
        # otherwise the text is drawn off the canvas and lost
        h = img.shape[0]
        if side > 0 and y > h * 0.78:
            side = -1
        elif side < 0 and y < h * 0.18:
            side = 1
        ty = y + (30 + th if side > 0 else -22)
        cv2.putText(img, text, (tx, ty), cv2.FONT_HERSHEY_SIMPLEX, 0.42,
                    (255, 255, 255), 3, cv2.LINE_AA)
        cv2.putText(img, text, (tx, ty), cv2.FONT_HERSHEY_SIMPLEX, 0.42,
                    (60, 60, 60), 1, cv2.LINE_AA)


if __name__ == "__main__":
    main()
