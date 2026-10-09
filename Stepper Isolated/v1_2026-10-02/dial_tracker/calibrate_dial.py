#!/usr/bin/env python3
"""Interactive calibration for `track_dial.py`: point at the dial, save a JSON.

Opens the video in a window with the current geometry drawn on top. You scrub
to a frame, pick what a click should set from the list on the left, and click
on the picture. Everything the tracker needs ends up in one JSON file.

What has to be set, roughly in this order:

  1. main centre     the pivot the long needle turns about
  2. main r_inner    inner edge of the ring the needle is measured in - must
                     clear the hub and the revolution counter
  3. main r_outer    outer edge - must stay inside the printed tick ring
  4. zero bearing    click along the direction of the printed 0
  5. sub centre      pivot of the revolution counter (skip with --no-sub)
  6. sub r_inner/outer, sub zero   same idea, for the small needle

`Auto centre` does step 1 properly once 2 and 3 are roughly right: it fits a
line to the needle in 80 frames and intersects them, which is far more precise
than clicking. `Fit sub scale` measures how far the small needle moves per turn
of the big one, so you do not have to read that off the dial.

Usage:
    python calibrate_dial.py --video "../recordings/Video Project 5.mp4"
    python calibrate_dial.py --video clip.mp4 --calib existing.json
    python calibrate_dial.py --video clip.mp4 --um-per-rev 100 --no-sub
"""

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import dial_core as dc  # noqa: E402

import matplotlib
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.widgets import Button, RadioButtons, Slider, TextBox  # noqa: E402

CLICK_TARGETS = [
    "main centre",
    "main r_inner",
    "main r_outer",
    "main zero bearing",
    "sub centre",
    "sub r_inner",
    "sub r_outer",
    "sub zero bearing",
]


class Calibrator:
    def __init__(self, video, cal):
        self.video = video
        self.cal = cal
        self.cap = dc.open_video(video)
        self.n = dc.video_info(video)["n_frames"]
        self.frame_i = 0
        self.bg = None
        self.detected = None       # (main_bearing, sub_bearing) for the shown frame
        self._build_ui()

    # ---------------------------------------------------------------- UI
    def _build_ui(self):
        self.fig = plt.figure(figsize=(12.5, 8.6))
        self.fig.canvas.manager.set_window_title("dial calibration — %s" % Path(self.video).name)
        self.ax = self.fig.add_axes([0.26, 0.14, 0.72, 0.82])
        self.ax.set_xticks([]); self.ax.set_yticks([])

        img = dc.read_gray(self.cap, 0)
        self.im = self.ax.imshow(img, cmap="gray", vmin=0, vmax=255)

        # overlay artists
        self.art = {
            "mc": self.ax.plot([], [], "+", color="#ff0", ms=16, mew=2)[0],
            "ri": self.ax.plot([], [], "-", color="#0cf", lw=1.1)[0],
            "ro": self.ax.plot([], [], "-", color="#0cf", lw=1.1)[0],
            "mz": self.ax.plot([], [], "-", color="#0f0", lw=1.6)[0],
            "mn": self.ax.plot([], [], "-", color="#f33", lw=2.0)[0],
            "sc": self.ax.plot([], [], "+", color="#ff0", ms=11, mew=2)[0],
            "sri": self.ax.plot([], [], "-", color="#f8f", lw=1.0)[0],
            "sro": self.ax.plot([], [], "-", color="#f8f", lw=1.0)[0],
            "sz": self.ax.plot([], [], "-", color="#0f0", lw=1.2)[0],
            "sn": self.ax.plot([], [], "-", color="#f0f", lw=2.0)[0],
        }
        self.readout = self.ax.text(0.015, 0.975, "", transform=self.ax.transAxes,
                                    va="top", ha="left", color="#ff0", fontsize=10,
                                    family="monospace",
                                    bbox=dict(fc="black", alpha=.55, ec="none"))

        ax_sl = self.fig.add_axes([0.26, 0.065, 0.72, 0.03])
        self.slider = Slider(ax_sl, "frame", 0, max(self.n - 1, 1), valinit=0, valstep=1)
        self.slider.on_changed(self.on_frame)

        ax_radio = self.fig.add_axes([0.015, 0.47, 0.225, 0.42])
        ax_radio.set_title("a click sets…", fontsize=9)
        self.radio = RadioButtons(ax_radio, CLICK_TARGETS)

        def mkbtn(x, y, w, h, label, cb):
            b = Button(self.fig.add_axes([x, y, w, h]), label)
            b.on_clicked(cb)
            return b

        self.b_auto = mkbtn(0.015, 0.385, 0.108, 0.045, "Auto centre", self.on_auto_centre)
        self.b_bg = mkbtn(0.132, 0.385, 0.108, 0.045, "Background", self.on_background)
        self.b_det = mkbtn(0.015, 0.330, 0.108, 0.045, "Detect", self.on_detect)
        self.b_sub = mkbtn(0.132, 0.330, 0.108, 0.045, "Fit sub scale", self.on_fit_sub)
        self.b_dir = mkbtn(0.015, 0.275, 0.108, 0.045, "Flip direction", self.on_flip)
        self.b_save = mkbtn(0.132, 0.275, 0.108, 0.045, "Save", self.on_save)

        ax_um = self.fig.add_axes([0.115, 0.205, 0.125, 0.042])
        self.tb_um = TextBox(ax_um, "µm / rev  ", initial=str(self.cal["main"]["um_per_rev"]))
        self.tb_um.on_submit(self.on_um)
        ax_mm = self.fig.add_axes([0.115, 0.150, 0.125, 0.042])
        self.tb_mm = TextBox(ax_mm, "mm / step  ", initial=str(self.cal["sub"]["mm_per_step"]))
        self.tb_mm.on_submit(self.on_mm)

        self.status = self.fig.text(0.015, 0.10, "", fontsize=8.5, family="monospace",
                                    va="top", wrap=True)
        self.fig.canvas.mpl_connect("button_press_event", self.on_click)
        self.set_status("click on the dial to place the geometry; 'Auto centre' once the "
                        "ring roughly covers the needle")
        self.refresh()

    def set_status(self, msg):
        self.status.set_text("\n".join(_wrap(msg, 46)))
        self.fig.canvas.draw_idle()

    # ------------------------------------------------------------ events
    def on_frame(self, val):
        self.frame_i = int(val)
        img = dc.read_gray(self.cap, self.frame_i)
        if img is not None:
            self.im.set_data(img)
        self.detected = None
        self.refresh()

    def on_click(self, ev):
        if ev.inaxes is not self.ax or ev.xdata is None:
            return
        x, y = float(ev.xdata), float(ev.ydata)
        t = self.radio.value_selected
        m, s = self.cal["main"], self.cal["sub"]
        if t == "main centre":
            m["center"] = [x, y]
        elif t == "main r_inner":
            m["r_inner"] = float(np.hypot(x - m["center"][0], y - m["center"][1]))
        elif t == "main r_outer":
            m["r_outer"] = float(np.hypot(x - m["center"][0], y - m["center"][1]))
        elif t == "main zero bearing":
            m["zero_bearing_deg"] = _bear(m["center"], x, y)
        elif t == "sub centre":
            s["center"] = [x, y]
        elif t == "sub r_inner":
            s["r_inner"] = float(np.hypot(x - s["center"][0], y - s["center"][1]))
        elif t == "sub r_outer":
            s["r_outer"] = float(np.hypot(x - s["center"][0], y - s["center"][1]))
        elif t == "sub zero bearing":
            s["zero_bearing_deg"] = _bear(s["center"], x, y)
        self.set_status("set %s" % t)
        self.refresh()

    def on_um(self, text):
        try:
            self.cal["main"]["um_per_rev"] = float(text)
            self.set_status("one turn of the long needle = %s µm" % text)
        except ValueError:
            self.set_status("not a number: %r" % text)

    def on_mm(self, text):
        try:
            self.cal["sub"]["mm_per_step"] = float(text)
        except ValueError:
            self.set_status("not a number: %r" % text)

    def on_flip(self, _):
        m = self.cal["main"]
        m["direction"] = "ccw" if m["direction"] == "cw" else "cw"
        self.set_status("increasing reading = %s" % m["direction"])
        self.refresh()

    def on_background(self, _):
        self.set_status("building background (two passes, ~1 min)…")
        det, m = self.cal["detect"], self.cal["main"]
        bg0 = dc.build_background(self.video, m["center"], None,
                                  det["bg_frames"], det["bg_percentile"])
        tmp = dict(self.cal, sub=dict(self.cal["sub"], enabled=False))
        bearings = dc.track(self.video, tmp, bg0, None)["main_bearing"]
        self.bg = dc.build_background(self.video, m["center"], bearings,
                                      det["bg_frames"], wedge_deg=5.0)
        self.set_status("background ready — 'Detect' now works on any frame")
        self.on_detect(None)

    def on_auto_centre(self, _):
        m = self.cal["main"]
        if self.bg is None:
            self.set_status("building a bootstrap background first…")
            self.bg = dc.build_background(self.video, m["center"], None,
                                          self.cal["detect"]["bg_frames"],
                                          self.cal["detect"]["bg_percentile"])
        c, resid = dc.fit_rotation_center(self.video, self.bg, m["center"],
                                          m["r_inner"], m["r_outer"])
        m["center"] = [float(c[0]), float(c[1])]
        self.set_status("centre fitted to (%.2f, %.2f), median residual %.2f px — "
                        "rebuild the background now" % (c[0], c[1], resid))
        self.refresh()

    def on_detect(self, _):
        if self.bg is None:
            self.set_status("press 'Background' first")
            return
        g = dc.read_gray(self.cap, self.frame_i)
        tr = _detect_one(g, self.bg, self.cal)
        self.detected = tr
        self.refresh()

    def on_fit_sub(self, _):
        """Measure the small needle's travel per turn of the big one."""
        if self.bg is None:
            self.set_status("press 'Background' first")
            return
        self.set_status("tracking both needles to fit the sub-dial scale…")
        tr = dc.track(self.video, self.cal, self.bg, self.bg)
        bad = dc.flag_bad_frames(tr)
        mu = dc.unwrap_bearings(tr["main_bearing"], bad)
        ok = np.isfinite(tr["sub_bearing"])
        if ok.sum() < 50:
            self.set_status("could not see the small needle — check its centre and radii")
            return
        ix = np.arange(len(mu))
        su = np.degrees(np.unwrap(np.radians(np.interp(ix, ix[ok], tr["sub_bearing"][ok]))))
        slope = np.polyfit(mu, su, 1)[0]
        deg_per_step = float(slope * 360.0)
        self.cal["sub"]["deg_per_step"] = round(deg_per_step, 3)
        resid = su - np.polyval(np.polyfit(mu, su, 1), mu)
        self.set_status("sub-dial moves %.3f deg per turn of the main needle "
                        "(scatter %.2f deg). Now set 'sub zero bearing' by clicking the "
                        "printed 0 on the counter." % (deg_per_step, resid.std()))

    def on_save(self, _):
        # keep the video path relative to the calibration when they share a
        # root, so the folder can be copied or moved as a unit
        try:
            rel = Path(self.video).resolve().relative_to(Path(self.out_path).resolve().parent)
            self.cal["video"] = rel.as_posix()
        except ValueError:
            try:
                import os
                self.cal["video"] = Path(os.path.relpath(
                    Path(self.video).resolve(),
                    Path(self.out_path).resolve().parent)).as_posix()
            except ValueError:
                self.cal["video"] = self.video
        info = dc.video_info(self.video)
        self.cal["fps"] = info["fps"]
        self.cal["n_frames"] = info["n_frames"]
        self.cal["frame_size"] = [info["width"], info["height"]]
        dc.save_calibration(self.out_path, self.cal)
        self.set_status("saved %s" % Path(self.out_path).name)

    # ----------------------------------------------------------- drawing
    def refresh(self):
        m, s = self.cal["main"], self.cal["sub"]
        th = np.linspace(0, 2 * np.pi, 400)
        self.art["mc"].set_data([m["center"][0]], [m["center"][1]])
        for key, r in (("ri", m["r_inner"]), ("ro", m["r_outer"])):
            self.art[key].set_data(m["center"][0] + r * np.cos(th),
                                   m["center"][1] + r * np.sin(th))
        self.art["mz"].set_data(*_ray(m["center"], m["zero_bearing_deg"], m["r_inner"], m["r_outer"]))
        self.art["sc"].set_data([s["center"][0]], [s["center"][1]])
        for key, r in (("sri", s["r_inner"]), ("sro", s["r_outer"])):
            self.art[key].set_data(s["center"][0] + r * np.cos(th),
                                   s["center"][1] + r * np.sin(th))
        self.art["sz"].set_data(*_ray(s["center"], s["zero_bearing_deg"], s["r_inner"], s["r_outer"]))

        txt = ["frame %d / %d" % (self.frame_i, self.n - 1),
               "main  c=(%.1f, %.1f)  r=%.0f..%.0f" % (m["center"][0], m["center"][1],
                                                       m["r_inner"], m["r_outer"]),
               "      zero %.2f deg   %s   %.1f um/rev" % (m["zero_bearing_deg"],
                                                           m["direction"], m["um_per_rev"]),
               "sub   c=(%.1f, %.1f)  r=%.0f..%.0f" % (s["center"][0], s["center"][1],
                                                       s["r_inner"], s["r_outer"]),
               "      zero %.2f deg   %.3f deg/step" % (s["zero_bearing_deg"], s["deg_per_step"])]
        if self.detected:
            mb, sb = self.detected
            self.art["mn"].set_data(*_ray(m["center"], mb, m["r_inner"] * .3, m["r_outer"]))
            if np.isfinite(sb):
                self.art["sn"].set_data(*_ray(s["center"], sb, s["r_inner"], s["r_outer"]))
            else:
                self.art["sn"].set_data([], [])
            face = float(dc.main_face_microns(np.array([mb]), self.cal)[0])
            smm = float(dc.sub_millimetres(np.array([sb]), self.cal)[0])
            txt.append("")
            txt.append("needle %.3f deg  ->  face %.2f um" % (mb, face))
            txt.append("counter %.2f deg  ->  %.4f mm" % (sb, smm) if np.isfinite(sb)
                       else "counter not found")
        else:
            self.art["mn"].set_data([], [])
            self.art["sn"].set_data([], [])
        self.readout.set_text("\n".join(txt))
        self.fig.canvas.draw_idle()


def _wrap(text, width):
    out, line = [], ""
    for word in text.split():
        if len(line) + len(word) + 1 > width:
            out.append(line); line = word
        else:
            line = (line + " " + word).strip()
    out.append(line)
    return out


def _bear(center, x, y):
    return float(np.degrees(np.arctan2(x - center[0], -(y - center[1]))) % 360.0)


def _ray(center, bearing_deg, r0, r1):
    b = np.radians(bearing_deg)
    return ([center[0] + r0 * np.sin(b), center[0] + r1 * np.sin(b)],
            [center[1] - r0 * np.cos(b), center[1] - r1 * np.cos(b)])


def _detect_one(gray, bg, cal):
    import cv2
    det, m, s = cal["detect"], cal["main"], cal["sub"]
    rm = dc.radius_map(gray.shape, m["center"]); am = dc.bearing_map(gray.shape, m["center"])
    ann = (rm > m["r_inner"]) & (rm < m["r_outer"])
    d = cv2.subtract(bg, gray); dv = d[ann]
    thr = max(det["main_thresh_min"], det["main_thresh_frac"] * float(dv.max()))
    sel = dv > thr
    mb = np.nan
    if sel.sum() >= det["min_pixels_main"]:
        mb, _, _ = dc.fit_bearing(am[ann][sel], dv[sel].astype(float),
                                  det["main_tol_deg"], det["min_pixels_main"])
    sb = np.nan
    if s.get("enabled", True):
        rs = dc.radius_map(gray.shape, s["center"]); asb = dc.bearing_map(gray.shape, s["center"])
        keep = (rs > s["r_inner"]) & (rs < s["r_outer"]) & (d > det["sub_thresh"])
        if np.isfinite(mb):
            keep &= np.abs(dc.angdiff(am, mb)) > det["wedge_exclude_deg"]
        if keep.sum() >= det["min_pixels_sub"]:
            sb, _, _ = dc.fit_bearing(asb[keep], d[keep].astype(float),
                                      det["sub_tol_deg"], det["min_pixels_sub"])
    return mb, sb


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--video", required=True)
    ap.add_argument("--calib", help="start from this calibration instead of defaults")
    ap.add_argument("--out", help="where to save (default: calibration_<video stem>.json)")
    ap.add_argument("--um-per-rev", type=float, help="micrometres per turn of the long needle")
    ap.add_argument("--no-sub", action="store_true", help="dial has no revolution counter")
    args = ap.parse_args()

    video = args.video
    if not Path(video).is_absolute() and not Path(video).exists() and args.calib:
        video = str((Path(args.calib).parent / video).resolve())

    cal = dc.load_calibration(args.calib) if args.calib else dict(dc.DEFAULT_CALIBRATION)
    cal = dc._merge(dc.DEFAULT_CALIBRATION, cal)
    info = dc.video_info(video)
    if not cal["main"]["r_outer"]:
        # a first guess that puts the ring somewhere sane on a square-ish frame
        h, w = info["height"], info["width"]
        cal["main"]["center"] = [w / 2.0, h / 2.0]
        cal["main"]["r_inner"] = min(h, w) * 0.23
        cal["main"]["r_outer"] = min(h, w) * 0.39
        cal["sub"]["center"] = [w / 2.0, h * 0.64]
        cal["sub"]["r_inner"] = min(h, w) * 0.028
        cal["sub"]["r_outer"] = min(h, w) * 0.070
    if args.um_per_rev:
        cal["main"]["um_per_rev"] = args.um_per_rev
    if args.no_sub:
        cal["sub"]["enabled"] = False

    c = Calibrator(video, cal)
    c.out_path = args.out or (Path(__file__).resolve().parent /
                              ("calibration_" + Path(video).stem.replace(" ", "_") + ".json"))
    print("calibrating %s  (%d frames)" % (video, info["n_frames"]))
    print("will save to %s" % c.out_path)
    plt.show()


if __name__ == "__main__":
    main()
