#!/usr/bin/env python3
"""Track a dial-indicator needle through a video and export the reading as CSV.

Reads a calibration written by `calibrate_dial.py`, measures both needles in
every frame, converts the long pointer's angle to micrometres and writes one
row per frame. Also writes a QC plot and, optionally, an overlay video so the
numbers can be checked against the picture they came from.

Usage:
    python track_dial.py --calib calibration_video_project_5.json
    python track_dial.py --calib cal.json --out readings.csv --qc qc.png
    python track_dial.py --calib cal.json --start auto      # skip a startup jolt
    python track_dial.py --calib cal.json --overlay check.mp4
    python track_dial.py --calib cal.json --rebuild-bg      # ignore cached background
"""

import argparse
import csv
import json
import sys
from pathlib import Path

import numpy as np

import dial_core as dc


def log(msg):
    print(msg, flush=True)


def cached_background(video, cal, cache_dir, rebuild=False):
    """Build (or reload) the two background images this video needs.

    They depend only on the video and the geometry, and take ~40 s to compute,
    so they are cached next to the calibration.
    """
    cache = Path(cache_dir) / "background.npz"
    key = json.dumps([cal["main"]["center"], cal["sub"]["center"],
                      cal["detect"]["bg_frames"], cal["detect"]["bg_percentile"]])
    if cache.exists() and not rebuild:
        z = np.load(cache, allow_pickle=True)
        if str(z.get("key", "")) == key:
            log(f"  background: reusing {cache.name}")
            return z["bg_main"], z["bg_sub"]

    det = cal["detect"]

    # Pass 1: a bright percentile is needle-free almost everywhere, which is
    # enough to find the needle. It degrades where the needle lingers, so it is
    # only ever used to bootstrap pass 2.
    log("  background pass 1/2 (bright percentile, bootstrap)")
    bg0 = dc.build_background(video, cal["main"]["center"], bearings=None,
                              n_frames=det["bg_frames"], percentile=det["bg_percentile"])

    log("            locating the needle with it")
    tmp = dict(cal)
    tmp["sub"] = dict(cal["sub"], enabled=False)
    bearings = dc.track(video, tmp, bg0, None)["main_bearing"]

    # Pass 2: now that the long needle's position is known in every frame, blank
    # a wedge around it and take a plain median. One image serves both dials:
    # the main annulus clears the sub-dial entirely, and the small needle sits
    # at any given angle in far fewer than half the frames, so the median steps
    # over it too.
    log("  background pass 2/2 (needle masked, median)")
    bg = dc.build_background(video, cal["main"]["center"], bearings=bearings,
                             n_frames=det["bg_frames"], wedge_deg=5.0)

    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, bg_main=bg, bg_sub=bg, key=key)
    log(f"  background: cached to {cache.name}")
    return bg, bg


def write_csv(path, rows, header, meta):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w", newline="", encoding="utf-8") as fh:
        for k, v in meta.items():
            fh.write(f"# {k}: {v}\n")
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def qc_plot(path, t, um, sub_mm, spread, contrast, bad, start, cal):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True,
                           gridspec_kw={"height_ratios": [3, 2, 2]})
    ax[0].plot(t, um, lw=1.0, color="#1f4e79")
    if bad.any():
        ax[0].plot(t[bad], um[bad], ".", ms=4, color="#c00", label="flagged frames")
        ax[0].legend(loc="best", fontsize=8)
    ax[0].set_ylabel("dial reading  [µm]")
    ax[0].set_title("Dial indicator read from video — %s" % Path(cal["video"]).name)
    ax[0].grid(alpha=.3)

    good = np.isfinite(sub_mm)
    if good.any():
        ax[1].plot(t[good], sub_mm[good] * 1000.0, lw=.8, color="#999",
                   label="revolution counter (independent)")
        ax[1].plot(t, um, lw=1.0, color="#1f4e79", alpha=.7, label="main needle")
        ax[1].legend(loc="best", fontsize=8)
    ax[1].set_ylabel("cross-check  [µm]")
    ax[1].grid(alpha=.3)

    ax[2].plot(t, spread, lw=.8, color="#b35", label="needle angular spread [°]")
    ax2 = ax[2].twinx()
    ax2.plot(t, contrast, lw=.8, color="#395", label="peak needle contrast [DN]")
    ax[2].set_ylabel("spread [°]", color="#b35")
    ax2.set_ylabel("contrast [DN]", color="#395")
    ax[2].set_xlabel("time  [s]")
    ax[2].grid(alpha=.3)
    if start > 0:
        for a in ax:
            a.axvline(t[0], color="#c00", ls="--", lw=.8)
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def make_overlay(video, out_path, cal, bearing, um, frames):
    """Burn the measured needle and reading into a short clip for eyeballing."""
    import cv2
    cap = dc.open_video(video)
    fps = cal["fps"] or 30.0
    h = w = None
    writer = None
    c = np.asarray(cal["main"]["center"], float)
    for i in frames:
        g = dc.read_gray(cap, i)
        if g is None:
            continue
        img = cv2.cvtColor(g, cv2.COLOR_GRAY2BGR)
        if writer is None:
            h, w = img.shape[:2]
            writer = cv2.VideoWriter(str(out_path), cv2.VideoWriter_fourcc(*"mp4v"),
                                     fps, (w, h))
        if np.isfinite(bearing[i]):
            th = np.radians(bearing[i])
            for r in (cal["main"]["r_inner"], cal["main"]["r_outer"]):
                p = (int(c[0] + r * np.sin(th)), int(c[1] - r * np.cos(th)))
                cv2.circle(img, p, 5, (0, 0, 255), -1)
            p0 = (int(c[0] + cal["main"]["r_inner"] * np.sin(th)),
                  int(c[1] - cal["main"]["r_inner"] * np.cos(th)))
            p1 = (int(c[0] + cal["main"]["r_outer"] * np.sin(th)),
                  int(c[1] - cal["main"]["r_outer"] * np.cos(th)))
            cv2.line(img, p0, p1, (0, 0, 255), 2)
        cv2.circle(img, (int(c[0]), int(c[1])), 4, (0, 255, 255), -1)
        cv2.putText(img, "f%d  %.2f um" % (i, um[i]), (20, 44),
                    cv2.FONT_HERSHEY_SIMPLEX, 1.1, (0, 255, 255), 2, cv2.LINE_AA)
        writer.write(img)
    cap.release()
    if writer:
        writer.release()


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--calib", required=True, help="calibration JSON from calibrate_dial.py")
    ap.add_argument("--video", help="override the video path stored in the calibration")
    ap.add_argument("--out", help="output CSV (default: output/<video stem>_dial.csv)")
    ap.add_argument("--qc", help="QC plot PNG (default: alongside the CSV)")
    ap.add_argument("--overlay", help="write an overlay video here for visual checking")
    ap.add_argument("--overlay-step", type=int, default=1, help="use every Nth frame in the overlay")
    ap.add_argument("--start", default=None,
                    help="'auto' to skip a startup disturbance, or a frame number")
    ap.add_argument("--end", type=int, default=None, help="last frame to export")
    ap.add_argument("--rebuild-bg", action="store_true", help="ignore the cached background")
    ap.add_argument("--drop-bad", action="store_true",
                    help="remove flagged frames from the CSV instead of marking them")
    args = ap.parse_args()

    cal_path = Path(args.calib)
    cal = dc.load_calibration(cal_path)
    video = args.video or cal["video"]
    if not video:
        ap.error("no video in the calibration and none given with --video")
    # a path stored in the calibration is relative to the calibration itself,
    # so the pair stays portable when the folder is moved or copied
    if not Path(video).is_absolute() and not Path(video).exists():
        video = str((cal_path.parent / video).resolve())
    info = dc.video_info(video)
    cal["fps"] = cal["fps"] or info["fps"]
    log(f"video   : {video}")
    log(f"          {info['n_frames']} frames, {info['fps']:.3f} fps, "
        f"{info['width']}x{info['height']}, {info['n_frames']/info['fps']:.2f} s")

    bg_main, bg_sub = cached_background(video, cal, cal_path.parent, args.rebuild_bg)

    log("tracking needles")
    tr = dc.track(video, cal, bg_main, bg_sub,
                  progress=lambda i, n: log(f"    frame {i}/{n}") if i % 400 == 0 else None)

    bad = dc.flag_bad_frames(tr)
    log(f"  flagged {int(bad.sum())} / {len(bad)} frames as untrustworthy")

    dist = dc.find_disturbance(bad)
    if dist:
        f0, f1, resume = dist
        log(f"  startup disturbance: frames {f0}-{f1} "
            f"(t {f0/info['fps']:.2f}-{f1/info['fps']:.2f} s), clean again from frame {resume}")

    # --- choose the exported window ------------------------------------
    start = cal["export"]["start_frame"]
    if args.start == "auto":
        start = dist[2] if dist else 0
    elif args.start is not None:
        start = int(args.start)
    end = args.end if args.end is not None else cal["export"]["end_frame"]
    if end is None or end < 0:
        end = len(bad) - 1
    log(f"  exporting frames {start}..{end}")

    unwrapped = dc.unwrap_bearings(tr["main_bearing"], bad)
    step = np.abs(np.diff(unwrapped))
    log(f"  max frame-to-frame needle step: {step.max():.2f} deg "
        f"({step.max()/360*cal['main']['um_per_rev']:.2f} um) - unwrap is safe below 180 deg")

    # The reference frame anchors the *absolute* reading, so it is deliberately
    # independent of the exported window: for this video it is a pre-jolt frame
    # where the dial is known to read a whole number of turns, which is outside
    # the window once a startup disturbance is skipped. Only move it if it does
    # not point at a usable frame.
    cal = dict(cal)
    cal["export"] = dict(cal["export"])
    ref = int(cal["export"]["reference_frame"])
    if not (0 <= ref < len(bad)) or bad[ref]:
        log(f"  reference frame {ref} is unusable; re-referencing to frame {start}")
        cal["export"]["reference_frame"] = start
        cal["export"]["reference_reading_um"] = 0.0
    um = dc.to_microns(unwrapped, cal)
    face = dc.main_face_microns(tr["main_bearing"], cal)
    occluded = dc.flag_sub_occlusion(tr["sub_pixels"])
    sub_b = np.where(occluded, np.nan, tr["sub_bearing"])
    sub_mm = dc.sub_millimetres(sub_b, cal)
    log(f"  revolution counter hidden by the long needle in {int(occluded.sum())} frames "
        f"(dropped from the cross-check only)")

    sel = np.arange(start, end + 1)
    if args.drop_bad:
        sel = sel[~bad[sel]]
    t = (sel - start) / info["fps"]

    rows = []
    for k, i in enumerate(sel):
        rows.append([i, f"{t[k]:.6f}", f"{um[i]:.4f}", f"{um[i]-um[sel[0]]:.4f}",
                     f"{tr['main_bearing'][i]:.4f}" if np.isfinite(tr["main_bearing"][i]) else "",
                     f"{face[i]:.3f}" if np.isfinite(face[i]) else "",
                     f"{sub_b[i]:.3f}" if np.isfinite(sub_b[i]) else "",
                     f"{sub_mm[i]:.5f}" if np.isfinite(sub_mm[i]) else "",
                     f"{tr['main_spread'][i]:.3f}" if np.isfinite(tr["main_spread"][i]) else "",
                     int(tr["main_pixels"][i]), int(tr["sub_pixels"][i]), int(bad[i])])

    out = Path(args.out) if args.out else \
        cal_path.parent / "output" / (Path(video).stem.replace(" ", "_") + "_dial.csv")
    meta = {
        "source_video": video,
        "calibration": cal_path.name,
        "fps": f"{info['fps']:.6f}",
        "frames_in_video": info["n_frames"],
        "exported_frames": len(sel),
        "start_frame": start,
        "end_frame": end,
        "um_per_rev": cal["main"]["um_per_rev"],
        "um_per_tick": cal["main"]["um_per_rev"] / cal["main"]["ticks_per_rev"],
        "reference_frame": cal["export"]["reference_frame"],
        "reference_reading_um": cal["export"]["reference_reading_um"],
        "flagged_frames_total": int(bad.sum()),
        "flagged_frames_dropped": "yes" if args.drop_bad else "no (marked in the flagged column)",
        "startup_disturbance": f"frames {dist[0]}-{dist[1]}, resumed {dist[2]}" if dist else "none",
        "sub_occluded_frames": int(occluded.sum()),
        "reading_um": "absolute dial reading (reference_reading_um at reference_frame)",
        "delta_um": "change since the first exported row",
    }
    header = ["frame", "t_s", "reading_um", "delta_um", "main_bearing_deg",
              "main_face_um", "sub_bearing_deg", "sub_reading_mm",
              "needle_spread_deg", "needle_pixels", "sub_pixels", "flagged"]
    write_csv(out, rows, header, meta)
    log(f"wrote {out}  ({len(rows)} rows)")

    qc = Path(args.qc) if args.qc else out.with_suffix(".png")
    qc_plot(qc, t, um[sel], sub_mm[sel], tr["main_spread"][sel],
            tr["main_contrast"][sel], bad[sel], start, cal)
    log(f"wrote {qc}")

    if args.overlay:
        frames = sel[::max(1, args.overlay_step)]
        log(f"writing overlay ({len(frames)} frames) ...")
        make_overlay(video, args.overlay, cal, tr["main_bearing"], um, frames)
        log(f"wrote {args.overlay}")

    # --- summary -------------------------------------------------------
    v = um[sel]
    log("")
    log("summary over the exported window")
    log(f"  duration        : {t[-1]:.3f} s   ({len(sel)} frames)")
    log(f"  reading         : {v.min():.2f} .. {v.max():.2f} um   (span {v.max()-v.min():.2f} um)")
    log(f"  start / end     : {v[0]:.2f} / {v[-1]:.2f} um")
    good = np.isfinite(sub_mm[sel])
    if good.sum() > 10:
        resid = v[good] - sub_mm[sel][good] * 1000.0
        log(f"  vs revolution counter: mean offset {resid.mean():.2f} um, "
            f"scatter {resid.std():.2f} um  (counter resolution ~{1/abs(cal['sub']['deg_per_step'])*cal['sub']['mm_per_step']*1000:.1f} um per deg)")


if __name__ == "__main__":
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    main()
