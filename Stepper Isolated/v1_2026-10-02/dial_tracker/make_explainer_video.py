#!/usr/bin/env python3
"""Render an explainer video: the dial being tracked, next to what it measures.

Three panels, 1920x1080 at the source frame rate:

    +----------------------+------------------------+
    |                      |  shaft outline (polar) |
    |   the dial, with     |  ideal circle + trace  |
    |   the tracking drawn +------------------------+
    |   on top             |  commanded vs measured |
    |                      |  angle, ideal + trace  |
    +----------------------+------------------------+

Left shows what the tracker actually sees: the two fitted pivots, the annuli
the needles are measured in, the pixels picked out as needle in this frame, and
the fitted bearings. Frames rejected as untrustworthy are called out rather
than hidden.

Right shows what comes out. Both deviations are small next to the ideal shape
they sit on, so both are magnified about the ideal - the same trick, and by
default the same 15x, as analysis/plot_encoder_full_rotation_xy.py.

The two panels come from two instruments recording the same single revolution:
the dial on the left, the shaft encoder below right. They are aligned by
fraction of a revolution, not by a hardware trigger.

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


# ---------------------------------------------------------------- encoder
def load_encoder(path, cpr, full_steps, microsteps):
    """Commanded and measured shaft angle from the encoder recording.

    Same reconstruction as analysis/plot_encoder_full_rotation_xy.py: the move
    is a uniform ramp of one microstep at a time across the recording, and the
    encoder angle is counts scaled by counts-per-revolution.
    """
    meta, lines = {}, []
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#"):
                k, sep, v = line[1:].partition(":")
                if sep:
                    meta[k.strip()] = v.strip()
            elif line.strip():
                lines.append(line)
    rows = list(csv.DictReader(lines))
    t_us = np.array([float(r["t_us"]) for r in rows])
    count = np.array([float(r["count"]) for r in rows])
    elapsed = (t_us - t_us[0]) * 1e-6
    duration = float(meta.get("clip_duration_s", elapsed[-1]))
    commands = full_steps * microsteps
    step_angle = 360.0 / commands
    commanded = np.rint(np.clip(elapsed / duration, 0, 1) * commands) * step_angle
    measured = (count - count[0]) * 360.0 / cpr
    return commanded, measured, measured - commanded, meta


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


def build_polar_panel(theta_deg, dev_um, gain, r0_um):
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
            lw=1.6, zorder=3, label="ideal: perfectly round and centred")
    ax.plot([], [], color=INK["trace"], lw=1.6, label=f"measured outline (×{gain:g})")
    ax.legend(loc="upper left", fontsize=8, framealpha=.92, ncol=1,
              bbox_to_anchor=(0.0, 1.02))
    ax.set_title("Shaft outline over one revolution", fontsize=11.5, pad=8)
    ax.text(0.0, -0.02, f"radial deviation about the ideal circle, µm, magnified ×{gain:g}",
            transform=ax.transAxes, ha="left", va="top", fontsize=8, color="#666")

    r = r0_um + gain * dev_um
    xs, ys = r * np.sin(np.radians(theta_deg)), r * np.cos(np.radians(theta_deg))
    base = _render(fig)
    px = _to_px(ax, xs, ys)
    centre_px = _to_px(ax, [0], [0])[0]
    plt.close(fig)
    return base, px, centre_px


def build_xy_panel(commanded, displayed, gain):
    """Static backdrop for commanded-vs-measured angle, plus a pixel mapper."""
    fig = _figure()
    ax = fig.add_axes([0.13, 0.13, 0.80, 0.76])
    lo = min(0.0, float(displayed.min()))
    hi = max(360.0, float(displayed.max()))
    m = max(4.0, 0.02 * (hi - lo))
    lim = (lo - m, hi + m)
    ax.plot(lim, lim, "--", color=INK["ideal"], lw=1.5,
            label="ideal: measured = commanded", zorder=2)
    ax.plot([], [], color=INK["trace"], lw=1.4,
            label=f"encoder measurement (×{gain:g})")
    ax.set(xlim=lim, ylim=lim, xlabel="commanded angle [deg]",
           ylabel="measured angle [deg]")
    ax.set_aspect("equal", adjustable="box")
    ax.grid(True, color=INK["grid"], lw=0.7, alpha=0.8)
    ax.set_title("Shaft encoder: commanded vs measured", fontsize=11.5, pad=10)
    ax.legend(loc="upper left", fontsize=8, framealpha=.92)
    base = _render(fig)
    px = _to_px(ax, commanded, displayed)
    plt.close(fig)
    return base, px


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
    _plate(img, 0, 0, 430, 48 + 27 * len(info["lines"]))
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
    ap.add_argument("--encoder", default="../recordings/encoder_full_rotation_20261009_162013.csv")
    ap.add_argument("--out", default="../recordings/dial_explainer.mp4",
                    help="default lands beside the source video, which git ignores")
    ap.add_argument("--scale", type=float, default=1.0,
                    help="resize the finished canvas, e.g. 0.667 for 720p")
    ap.add_argument("--gain", type=float, default=15.0, help="deviation magnification")
    ap.add_argument("--r0-um", type=float, default=5000.0,
                    help="drawing radius of the ideal circle [µm]")
    ap.add_argument("--encoder-cpr", type=float, default=20000.0)
    ap.add_argument("--full-steps", type=int, default=200)
    ap.add_argument("--microsteps", type=int, default=16)
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
    enc_path = Path(args.encoder)
    if not enc_path.is_absolute():
        enc_path = (cal_path.parent / enc_path).resolve()

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
    dev = um - np.mean(um[good])
    print(f"  one revolution = {period:.3f} s  (shaft angle from the dial's own period)")
    print(f"  deviation about the mean circle: {dev[good].min():+.1f} .. {dev[good].max():+.1f} µm")

    # --- encoder, aligned by fraction of a revolution ---------------------
    commanded, measured, error, meta = load_encoder(
        enc_path, args.encoder_cpr, args.full_steps, args.microsteps)
    displayed = commanded + args.gain * error
    print(f"encoder : {enc_path.name}  ({len(commanded)} samples, "
          f"endpoint error {error[-1]:+.3f} deg)")

    # map each video frame onto an encoder sample by revolution fraction
    frac = np.clip((t - t[export_start]) / period, 0, 1)
    enc_idx = np.rint(frac * (len(commanded) - 1)).astype(int)

    # --- static backdrops --------------------------------------------------
    polar_base, polar_px, _ = build_polar_panel(theta[good], dev[good], args.gain, args.r0_um)
    xy_base, xy_px = build_xy_panel(commanded, displayed, args.gain)

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
            (f"dial reading     {um[i]:8.2f} um", (0, 255, 255)),
            (f"shaft angle      {max(theta[i],0):8.2f} deg" if started
             else "shaft angle         -- (before start)", BGR["zero"]),
        ]
        left = draw_dial(frame, bg, cal, tr["main_bearing"][i], tr["sub_bearing"][i],
                         {"title": "What the tracker sees", "lines": lines,
                          "flagged": bool(bad[i])})
        k = max(0, i - export_start + 1) if started else 0
        polar = draw_trace(polar_base, polar_px, min(k, len(polar_px)),
                           _bgr(INK["trace"]), _bgr(INK["live"]))
        xy = draw_trace(xy_base, xy_px, enc_idx[i] + 1 if started else 0,
                        _bgr(INK["trace"]), _bgr(INK["live"]), width=2)
        canvas = np.zeros((2 * PANEL_H, VIDEO_W + PANEL_W, 3), np.uint8)
        canvas[:, :VIDEO_W] = left
        canvas[0:PANEL_H, VIDEO_W:] = polar
        canvas[PANEL_H:, VIDEO_W:] = xy
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


if __name__ == "__main__":
    main()
