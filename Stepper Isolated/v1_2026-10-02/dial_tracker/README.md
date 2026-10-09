# dial_tracker — reading an analogue dial indicator out of a video

Turns a video of a dial indicator into a CSV of micrometres against time, by
measuring the needle's angle in every frame. Built for
[`../recordings/Video Project 5.mp4`](../recordings), a 61.5 s clip of a
**Mitutoyo 513-425** lever indicator, but the geometry lives in a JSON file so
any dial works.

This is the same idea as [Tracker](https://opensourcephysics.github.io/tracker-website/)
— mark a feature, follow it frame by frame — specialised to one job. Tracker
follows a *point* and would have you click the needle tip 1 845 times, or fit a
template that slides along a straight needle. A dial has one degree of freedom,
not two, so this measures the one thing that matters, the **angle**, from every
needle pixel at once. That is what buys the precision: roughly 1 070 pixels go
into each estimate, and over the 13 frames before the stage moves the measured
bearing repeats to **0.0049° — 0.0027 µm**, about 700× finer than the dial's own
2 µm graduation. The instrument, not the tracking, is the limit (see §7).

```
dial_tracker/
├── README.md                          ← this document
├── dial_core.py                       ← the measurement: background, centre, angle, unwrap
├── calibrate_dial.py                  ← interactive window: point at the dial, save a JSON
├── track_dial.py                      ← batch run: JSON + video → CSV + QC plot
├── calibration_video_project_5.json   ← the calibration for the clip above
└── output/
    ├── Video_Project_5_dial.csv       ← 1 827 rows, the extracted reading
    └── Video_Project_5_dial.png       ← QC plot (reading, cross-check, quality)
```

`background.npz` is a cache built on first run (~1 min) and is not in git.

---

## 1. Reading this dial

| Part | What it means |
|---|---|
| **Main pointer** | One full turn = **200 µm**. 100 printed ticks, so **2 µm per tick**. |
| **Numbering** | 0–100–0: counts up from 0 in *both* directions and meets at 100 opposite the zero, so a deviation either side of zero reads directly. |
| **Revolution counter** | Small dial marked 0, .2, .4, .6 mm. One step = one full turn of the main pointer. Tells you which turn you are on, out of ~0.6 mm of range. |

Reading = (revolution counter, in whole 0.2 mm turns) + (main pointer, in µm).

Because the main scale is symmetric, the printed number alone is ambiguous: 60
could be 60 µm or 140 µm. The tracker never reads printed numbers — it works in
continuous angle and unwraps across frames, so it does not have that problem.
The CSV reports both: `reading_um` is the unambiguous continuous value, and
`main_face_um` is what the long needle alone points at (0–200 µm), which is what
you would compare against a freeze-frame.

---

## 2. Quickstart

```bash
pip install opencv-python-headless numpy matplotlib

# A) use the calibration that is already here
python track_dial.py --calib calibration_video_project_5.json --start auto

# B) calibrate a new video, then track it
python calibrate_dial.py --video "../recordings/My Clip.mp4" --um-per-rev 200
python track_dial.py --calib calibration_My_Clip.json --start auto
```

Useful flags on `track_dial.py`:

| Flag | Effect |
|---|---|
| `--start auto` | Skip a start-of-clip disturbance and begin where the needle is trackable again |
| `--start 120` / `--end 900` | Pick the window by hand |
| `--drop-bad` | Remove flagged frames instead of marking them |
| `--overlay check.mp4` | Burn the measured needle into a video so you can watch it track |
| `--overlay-step 6` | Use every 6th frame in that overlay |
| `--rebuild-bg` | Rebuild the background cache (needed after moving the centre) |

### Calibrating

`calibrate_dial.py` opens a window. Pick what a click should set from the list
on the left, then click the picture. Order that works:

1. **main centre** — roughly, on the hub
2. **main r_inner / r_outer** — a ring that covers the needle but clears the hub,
   the revolution counter and the printed tick ring
3. **Auto centre** — fits the pivot from 80 frames of needle lines; far better
   than clicking (here: ±1.3 px)
4. **Background** — two passes, ~1 min; `Detect` works afterwards
5. **main zero bearing** — click along the direction of the printed 0
6. **sub centre / radii**, then **Fit sub scale**, then **sub zero bearing**
7. **Save**

`Detect` draws the measured needle on the current frame — scrub around and
check it sits on the needle before trusting a whole run.

---

## 3. How it works

**Background subtraction, not edge finding.** The face is covered in ticks,
numerals and a maker's name that are exactly as dark as the needle. The needle
is the only thing that *moves*, so a per-pixel image of "what this pixel looks
like when the needle is elsewhere" erases every printed mark at once.

It is built in two passes. The first takes a bright percentile over 240 frames,
which is needle-free almost everywhere — but it fails where the needle lingers,
and a needle at rest is exactly where a measurement starts. So the second pass
uses the first pass's answer to blank a ±5° wedge around the needle in each
frame and takes a plain median of what is left. That halved the apparent needle
width (p95 spread 0.84° → 0.56°).

**Angle from every needle pixel.** In the annulus, pixels darker than the frame's
own adaptive threshold are taken as needle, and their intensity-weighted
circular mean is the bearing, re-centred a few times while discarding anything
more than 12° out. The threshold is a fraction of each frame's peak needle
contrast, so a frame where the needle is faint still yields an answer.

**Unwrapping.** Per-frame bearings are 0–360; `np.unwrap` turns them into a
continuous angle, which is what makes multi-turn readings possible. This is
only valid while the needle moves less than half a turn between frames — the
run prints the largest step so you can check. Here it is **2.47° (1.37 µm)**
against a 180° limit, so there is a 70× margin.

**Quality flags.** A frame is flagged when the needle's angular spread jumps or
its peak contrast collapses. Both mean the same thing physically: the needle
moved far within one exposure, so the frame holds a faint smear instead of a
crisp line.

---

## 4. The calibration for this clip, and how it was obtained

Nothing below was typed in from the dial's datasheet — all of it was fitted
from the video, which is why it can be trusted.

| Quantity | Value | How |
|---|---|---|
| Camera motion | < 0.25 px, < 0.005°, scale within 0.02% | ORB + RANSAC against frame 0; the camera is fixed, so no stabilisation is applied |
| Main pivot | **(536.28, 605.36)** px | 75 needle lines intersected, Cauchy-weighted; median residual **1.27 px** |
| Scale | **100 ticks/rev = 3.600°** exactly | FFT of the angular darkness profile over the tick ring: the 100-cycle harmonic dominates. 100 × 2 µm = 200 µm/turn, confirming the printed scale |
| Measuring ring | r = 250–420 px | The counter lies entirely inside r = 223 px, so r ≥ 250 isolates the main needle |
| Counter pivot | **(628.41, 696.46)** px | Same line-intersection fit, main needle masked; residual **0.72 px** |
| Counter scale | **−90.007° per 0.2 mm step** | Regression of counter angle on main angle over all 1 845 frames. Essentially exactly 90°, confirming the printed 0/.2/.4/.6 layout at the quarter points |
| Dial zero | bearing **308.925°** | The needle's pre-jolt rest position, which sits **0.54 µm** from a printed tick |
| Absolute reading at rest | **200.0 µm** | Main needle on the printed 0 *and* the counter one full step along — two independent features agreeing on 1 turn |

An earlier value of −86.26°/step came from a sub-dial centre that was 4 px off,
pulled there by the main needle sweeping through the fit. Masking the long
needle first moved the centre and the scale snapped to 90.007°. Worth
remembering: on a small dial, a few pixels of centre error shows up as several
degrees of scale error.

---

## 5. Results for this clip

1 827 rows, 60.867 s, from frame 18. The stage sweeps smoothly out and back:

| | Value |
|---|---|
| Range | **106.70 → 425.22 µm** (span **318.52 µm**, about 1.6 turns) |
| Start / end of window | 208.18 µm → 203.62 µm |
| Frames flagged | 7 of 1 845 |
| Agreement with the revolution counter | mean offset **−0.23 µm**, scatter **2.40 µm** |

That last line is the real check. The counter is a physically separate needle on
separate gearing, measured through a separate fit, and it tracks the main needle
to within 2.4 µm across the whole run — which is the counter's own readability
(≈2.2 µm per degree), not a tracking error. The main reading never uses it.

### The initialisation jolt

Frames **13–17** (t = 0.43–0.57 s) are corrupt: when the motor energises the
needle snaps across so fast that it smears, and peak contrast collapses from
~140 to 21 while the apparent needle width grows 10×. `--start auto` finds that
burst, drops it, and begins at **frame 18**, the first of a sustained clean run.

Two things worth knowing about it:

- **No revolution was lost.** This was the real risk — a needle that moves more
  than half a turn between usable frames cannot be unwrapped. The counter settles
  it: it moved **−3.99°** across the jolt, which is **+9.3 µm**, while the main
  needle's shortest-path change is **+14.73° = +8.18 µm**. A lost turn would have
  shown ~86° on the counter. The two agree to ~1 µm, which is the counter's
  resolution.
- **The needle does not come back to where it started.** It resumes **8.2 µm
  above** its pre-jolt rest and climbs from there — the jolt leaves a permanent
  offset rather than settling back. So the export begins at the first clean
  frame, not at a return to the original reading, and `reading_um` is anchored
  to the pre-jolt rest via `reference_frame: 5` so the absolute scale survives
  the gap. If you would rather have the window start at zero, use `delta_um`.

The clean frames 0–12 before the jolt are still in the video and are a usable
at-rest baseline; set `--start 0` to include them.

---

## 6. CSV columns

| Column | Meaning |
|---|---|
| `frame` | Frame index in the source video |
| `t_s` | Seconds from the first exported frame |
| `reading_um` | **The measurement.** Absolute dial reading, continuous across turns |
| `delta_um` | Change since the first exported row |
| `main_bearing_deg` | Raw needle bearing, 0 = 12 o'clock, clockwise positive |
| `main_face_um` | What the long needle alone points at, 0–200 µm |
| `sub_bearing_deg` | Revolution counter bearing; blank where the long needle hides it |
| `sub_reading_mm` | Independent coarse reading from the counter |
| `needle_spread_deg` | Apparent needle half-width — the main quality indicator |
| `needle_pixels`, `sub_pixels` | Pixels that went into each fit |
| `flagged` | 1 = not trustworthy (value is interpolated) |

Header lines beginning `#` record the source video, calibration, scale and the
disturbance that was skipped.

---

## 7. Limits

- **The camera must not move.** It does not here, and that was checked rather
  than assumed, but a handheld clip needs per-frame registration first.
- **Between two usable frames the needle must turn less than 180°.** Fast motion
  is not merely imprecise, it is ambiguous. The printed max step is the check,
  and the revolution counter is the backstop.
- **Accuracy is limited by the instrument, not the tracking.** With the stage at
  rest the reading repeats to **0.0027 µm** (0.0049°). While the needle is
  moving, the scatter about a smooth local fit is **0.11 µm** — and most of that
  is real needle dynamics and motion blur, not measurement error. Either way the
  dial's own graduation is 2 µm and its mechanism has hysteresis and friction
  well above that, so treat the extra decimal places as *resolution*, not
  accuracy: they are good for seeing how something changes, not for claiming an
  absolute position to a nanometre.
- **The counter is hidden twice per turn.** The long needle passes within 3 px of
  its centre; 85 frames are dropped from the cross-check for that reason. The
  main reading is unaffected.
- **`--start auto` only looks in the first quarter of the clip,** since it is
  meant for a power-on transient. Use `--start N` for anything else.
