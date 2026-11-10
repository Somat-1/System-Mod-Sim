#!/usr/bin/env python3
"""Render and interactively inspect an ISO 230-2 TwinCAT scope export.

The measurement channel is kept raw.  The only correction is exact modulo
2**32 unwrapping of the EL5101 UINT32 counter.  Displacement is shown relative
to the first scope sample with the nominal IDS scale of 1 count = 1 nm.
"""

from __future__ import annotations

import argparse
import re
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
ISO_DIR = HERE.parent
DEFAULT_SOURCE = ISO_DIR.parent / "ISO230-2_v2.csv"
DEFAULT_LOG = ISO_DIR / "data" / "terminal_20261021_101911.log"
DEFAULT_OUTPUT = HERE / "plots" / "iso230_2_v2_raw_measurement.png"

COUNTS_PER_MM = 1_000_000.0  # nominal IDS/EL5101 scale: 1 count = 1 nm
U32 = 2**32
U31 = 2**31
U256_PER_MM = 25_600.0       # 200 full steps/rev, 2 mm lead, u256 units


@dataclass
class ScopeData:
    index: np.ndarray
    counter_raw: np.ndarray
    counter_unwrapped: np.ndarray
    time_s: np.ndarray
    displacement_mm: np.ndarray
    sample_ms: float
    start_epoch_s: float
    start_local_text: str
    end_local_text: str
    channel: str
    rollover_count: int


@dataclass
class ControllerData:
    sample_time_s: np.ndarray
    command_mm: np.ndarray
    sample_state: list[str]
    state_changes: list[tuple[float, str]]
    stops: list[tuple[float, str, int, int, int, float]]
    alignment_hint_s: float
    alignment_offset_s: float = 0.0
    alignment_mad_mm: float = float("nan")


def _filetime_to_unix(ticks: int) -> float:
    return (ticks - 116_444_736_000_000_000) / 10_000_000.0


def read_scope(path: Path) -> ScopeData:
    raw = path.read_bytes()
    match = re.search(rb"(?m)^0\t-?\d+\s*$", raw)
    if not match:
        raise ValueError(f"No TwinCAT sample rows found in {path}")

    header = raw[: match.start()].decode("latin-1", errors="replace")
    metadata: dict[str, list[str]] = {}
    for line in header.splitlines():
        parts = line.split("\t")
        if parts and parts[0].strip():
            metadata[parts[0].strip()] = [p.strip() for p in parts[1:]]

    body = raw[match.start() :].replace(b"EOF", b"")
    values = np.fromstring(body.replace(b"\t", b" "), sep=" ", dtype=np.int64)
    if values.size % 2:
        raise ValueError(f"Odd number of values in sample table: {values.size}")
    table = values.reshape(-1, 2)
    index = table[:, 0]
    counter = table[:, 1]
    index_delta = np.diff(index)
    if np.any(index_delta <= 0):
        raise ValueError("Scope sample index is not strictly increasing")

    sample_ms = float(metadata.get("SampleTime[ms]", ["1"])[0])
    delta_raw = np.diff(counter)
    rollover_count = int(np.count_nonzero(np.abs(delta_raw) > U31))
    delta = (delta_raw + U31) % U32 - U31
    unwrapped = np.empty(counter.size, dtype=np.int64)
    unwrapped[0] = counter[0]
    unwrapped[1:] = counter[0] + np.cumsum(delta, dtype=np.int64)

    start_fields = metadata.get("Starttime of export", [])
    if not start_fields:
        raise ValueError("Scope header has no export start time")
    start_epoch_s = _filetime_to_unix(int(start_fields[0]))
    start_local = " ".join(start_fields[1:])
    end_fields = metadata.get("Endtime of export", [])
    end_local = " ".join(end_fields[1:])
    channel = metadata.get("SymbolName", [""])[0]
    # Preserve the scope clock across occasional missing samples; do not interpolate.
    time_s = (index - index[0]).astype(np.float64) * sample_ms / 1000.0
    displacement_mm = (unwrapped - unwrapped[0]) / COUNTS_PER_MM

    return ScopeData(
        index=index,
        counter_raw=counter,
        counter_unwrapped=unwrapped,
        time_s=time_s,
        displacement_mm=displacement_mm,
        sample_ms=sample_ms,
        start_epoch_s=start_epoch_s,
        start_local_text=start_local,
        end_local_text=end_local,
        channel=channel,
        rollover_count=rollover_count,
    )


def read_controller_log(path: Path, scope_start_epoch_s: float) -> ControllerData:
    sample_time, command, sample_state = [], [], []
    state_changes: list[tuple[float, str]] = []
    stops: list[tuple[float, str, int, int, int, float]] = []
    clock_hints = []

    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = re.match(r"^([0-9]+(?:\.[0-9]+)?)\s+(.*)$", line)
        if not match:
            continue
        host_s = float(match.group(1))
        host_elapsed_s = host_s - scope_start_epoch_s
        payload = match.group(2)
        if payload.startswith("# S,"):
            fields = payload.split(",")
            if len(fields) >= 6:
                board_s = float(fields[1]) / 1000.0
                sample_time.append(board_s)
                command.append(float(fields[2]) / U256_PER_MM)
                sample_state.append(fields[5].strip())
                clock_hints.append(host_elapsed_s - board_s)
        elif payload.startswith("# EVENT,"):
            fields = payload.split(",")
            if len(fields) >= 5 and fields[3] == "state_change":
                state_changes.append((float(fields[1]) / 1000.0, fields[2]))
        elif payload.startswith("# STOP,"):
            fields = payload.split(",")
            if len(fields) >= 11:
                state = fields[1]
                run, direction, target_i = map(int, fields[2:5])
                x_mm = float(fields[6])
                # Host time belongs to record_end. Place the marker at the
                # middle of the board's recorded window.
                record_mid_s = (int(fields[8]) + int(fields[9])) / 2000.0
                stops.append((record_mid_s, state,
                              run, direction, target_i, x_mm))

    return ControllerData(
        sample_time_s=np.asarray(sample_time),
        command_mm=np.asarray(command),
        sample_state=sample_state,
        state_changes=state_changes,
        stops=stops,
        alignment_hint_s=float(np.median(clock_hints)) if clock_hints else 0.0,
    )


def align_controller_to_scope(scope: ScopeData,
                              controller: ControllerData) -> ControllerData:
    """Align the board millisecond clock to the scope from the trajectories.

    Serial host timestamps lag because the terminal consumes buffered output.
    The board clock itself is stable, so only one additive offset is fitted.
    """
    states = np.asarray(controller.sample_state)
    use = states == "ISO"
    board_t = controller.sample_time_s[use]
    command = controller.command_mm[use]
    if board_t.size < 100:
        offset = controller.alignment_hint_s
        mad = float("nan")
    else:
        stride = max(1, board_t.size // 5000)
        board_t = board_t[::stride]
        command = command[::stride]

        def score(candidate: float) -> float:
            sample_t = board_t + candidate
            valid = ((sample_t >= scope.time_s[0]) &
                     (sample_t <= scope.time_s[-1]))
            if np.count_nonzero(valid) < 100:
                return float("inf")
            residual = (np.interp(sample_t[valid], scope.time_s,
                                  scope.displacement_mm) - command[valid])
            residual -= np.median(residual)
            return float(np.median(np.abs(residual)))

        coarse = np.arange(controller.alignment_hint_s - 40.0,
                           controller.alignment_hint_s + 40.0001, 0.1)
        coarse_scores = np.asarray([score(x) for x in coarse])
        coarse_best = float(coarse[int(np.argmin(coarse_scores))])
        fine = np.arange(coarse_best - 0.15, coarse_best + 0.1501, 0.005)
        fine_scores = np.asarray([score(x) for x in fine])
        best_i = int(np.argmin(fine_scores))
        offset = float(fine[best_i])
        mad = float(fine_scores[best_i])

    controller.sample_time_s = controller.sample_time_s + offset
    controller.state_changes = [(t + offset, state)
                                for t, state in controller.state_changes]
    controller.stops = [(t + offset, state, run, direction, target_i, x_mm)
                        for t, state, run, direction, target_i, x_mm
                        in controller.stops]
    controller.alignment_offset_s = offset
    controller.alignment_mad_mm = mad
    return controller


def _relative_command(controller: ControllerData, scope: ScopeData) -> np.ndarray:
    if controller.command_mm.size == 0:
        return controller.command_mm
    # At scope start the stage is held at x_ref.  Use the first controller
    # sample inside the scope interval as the matching command reference.
    inside = np.flatnonzero(controller.sample_time_s >= 0)
    ref_i = int(inside[0]) if inside.size else 0
    return controller.command_mm - controller.command_mm[ref_i]


def _phase_intervals(controller: ControllerData, duration_s: float):
    changes = [(t, s) for t, s in controller.state_changes
               if -0.5 <= t <= duration_s + 0.5]
    intervals = []
    if not changes or changes[0][0] > 0:
        intervals.append((0.0, max(0.0, changes[0][0] if changes else duration_s),
                          "HOLD"))
    for i, (start, state) in enumerate(changes):
        if start > duration_s:
            break
        end = changes[i + 1][0] if i + 1 < len(changes) else duration_s
        intervals.append((max(0.0, start), min(duration_s, end), state))
    return [(a, b, s) for a, b, s in intervals if b > a]


def configure_style() -> None:
    mpl.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.linewidth": 0.8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "xtick.major.size": 4,
        "ytick.major.size": 4,
        "legend.frameon": False,
        "figure.facecolor": "white",
        "axes.facecolor": "white",
        "savefig.facecolor": "white",
    })


PHASE_COLORS = {
    "HOLD": "#d9d9d9",
    "DATUM": "#91bfdb",
    "ISO": "#fee090",
    "END_REF": "#fc8d59",
    "PARK": "#b3de69",
    "IDLE": "#d9d9d9",
}


def render_static(scope: ScopeData, controller: ControllerData, source: Path,
                  output: Path) -> None:
    configure_style()
    fig = plt.figure(figsize=(14.2, 7.6), layout="constrained")
    gs = fig.add_gridspec(2, 1, height_ratios=[12, 1.35], hspace=0.04)
    ax = fig.add_subplot(gs[0])
    strip = fig.add_subplot(gs[1], sharex=ax)

    x_min = scope.time_s / 60.0
    raw_line, = ax.plot(x_min, scope.displacement_mm, color="#155f8a",
                        lw=0.72,
                        label=f"Raw interferometer ({scope.index.size:,} samples)",
                        rasterized=True, zorder=3)
    command_rel = _relative_command(controller, scope)
    command_line = None
    if command_rel.size:
        mask = ((controller.sample_time_s >= 0) &
                (controller.sample_time_s <= scope.time_s[-1]))
        command_line, = ax.plot(controller.sample_time_s[mask] / 60,
                                command_rel[mask], color="#5c5c5c", lw=0.8,
                                alpha=0.72, linestyle=(0, (3, 2)),
                                label="Controller command (context)", zorder=2)

    ax.axhline(0, color="0.55", lw=0.65, zorder=1)
    ax.set_ylabel("Relative interferometer displacement [mm]")
    ax.set_title("ISO 230-2 sequence â€” raw interferometer measurement", loc="left",
                 fontweight="bold", pad=30)
    ax.grid(True, color="0.9", lw=0.55)
    ax.set_xlim(0, scope.time_s[-1] / 60)
    handles = [raw_line] + ([command_line] if command_line is not None else [])
    ax.legend(handles=handles, loc="upper right", ncol=2)

    subtitle = (
        f"{source.name}  |  {scope.index.size:,} samples at "
        f"{1000 / scope.sample_ms:.0f} Hz  |  duration {scope.time_s[-1] / 60:.3f} min  |  "
        f"export {scope.start_local_text}â€“{scope.end_local_text.split()[-1] if scope.end_local_text else ''}"
    )
    ax.text(0, 1.012, subtitle, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.7, color="0.27")
    note = (
        f"Display transform only: {scope.rollover_count} exact UINT32 rollover(s) unwrapped; "
        "first sample subtracted; nominal 1 count = 1 nm.  No filtering, averaging, "
        "detrending, calibration, or eccentricity compensation.\nController context "
        f"aligned by board-clock trajectory fit (offset {controller.alignment_offset_s:+.3f} s)."
    )
    ax.text(0.005, 0.018, note, transform=ax.transAxes, ha="left", va="bottom",
            fontsize=8.0, color="0.27",
            bbox=dict(boxstyle="round,pad=0.28", fc="white", ec="0.78", alpha=0.92))

    intervals = _phase_intervals(controller, scope.time_s[-1])
    for start, end, state in intervals:
        strip.axvspan(start / 60, end / 60,
                      color=PHASE_COLORS.get(state, "#eeeeee"), alpha=0.82,
                      ec="white", lw=0.5)
        if state == "ISO" or end - start > 35:
            strip.text((start + end) / 120, 0.72, state.replace("_", " "),
                       ha="center", va="center", fontsize=7.5, color="0.22")

    stop_styles = {1: ("^", "#0072b2"), -1: ("v", "#d55e00"), 0: ("o", "#555555")}
    for direction in (1, -1, 0):
        times = [s[0] / 60 for s in controller.stops
                 if s[3] == direction and 0 <= s[0] <= scope.time_s[-1]]
        if times:
            marker, color = stop_styles[direction]
            label = {1: "outward stop", -1: "return stop", 0: "reference stop"}[direction]
            strip.scatter(times, np.full(len(times), 0.22), s=7, marker=marker,
                          c=color, linewidths=0, label=label, zorder=4)

    strip.set_ylim(0, 1)
    strip.set_yticks([])
    strip.set_xlabel("Elapsed time from scope export start [min]")
    strip.spines[["left", "right", "top"]].set_visible(False)
    strip.legend(loc="upper center", bbox_to_anchor=(0.5, -0.78), ncol=3,
                 fontsize=8, handletextpad=0.4, columnspacing=1.2)

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


class InteractiveViewer:
    def __init__(self, scope: ScopeData, controller: ControllerData, source: Path):
        configure_style()
        self.scope = scope
        self.controller = controller
        self.source = source
        self.auto_y = True
        self.command_visible = True
        self.grid_visible = True

        self.fig, (self.ax_pos, self.ax_step) = plt.subplots(
            2, 1, figsize=(15.5, 8.8), sharex=True,
            gridspec_kw={"height_ratios": [3.4, 1.25]}, layout="constrained")
        self.fig.canvas.manager.set_window_title("ISO 230-2 raw measurement viewer")

        self.x = scope.time_s
        self.y = scope.displacement_mm
        self.dy_um = np.r_[np.nan, np.diff(scope.counter_unwrapped) / 1000.0]
        self.raw_line, = self.ax_pos.plot(self.x, self.y, color="#155f8a", lw=0.68,
                                          label="raw interferometer", rasterized=True)
        command_rel = _relative_command(controller, scope)
        self.command_line = None
        if command_rel.size:
            mask = ((controller.sample_time_s >= 0) &
                    (controller.sample_time_s <= scope.time_s[-1]))
            self.command_line, = self.ax_pos.plot(
                controller.sample_time_s[mask], command_rel[mask], color="#555555",
                lw=0.85, alpha=0.7, linestyle=(0, (3, 2)), label="controller command")
        self.step_line, = self.ax_step.plot(self.x, self.dy_um, color="#8a3b12",
                                            lw=0.5, rasterized=True)

        self.ax_pos.set_ylabel("Relative displacement [mm]")
        self.ax_step.set_ylabel("Raw increment\n[Âµm/sample]")
        self.ax_step.set_xlabel("Elapsed time from scope export start [s]")
        self.ax_pos.set_title("ISO 230-2 â€” interactive raw-data viewer", loc="left",
                              fontweight="bold")
        self.ax_pos.legend(loc="upper right")
        for ax in (self.ax_pos, self.ax_step):
            ax.grid(True, color="0.9", lw=0.55)
        self.ax_pos.axhline(0, color="0.55", lw=0.65)
        self.ax_step.axhline(0, color="0.55", lw=0.65)
        self.ax_pos.set_xlim(0, scope.time_s[-1])

        help_text = (
            "Toolbar: pan / box zoom / home   â€¢   scroll: horizontal zoom at cursor   â€¢   "
            "A: toggle visible-window auto-Y   â€¢   Y: rescale Y now   â€¢   R: reset   â€¢   "
            "C: command on/off   â€¢   G: grid on/off"
        )
        self.fig.text(0.5, 0.004, help_text, ha="center", va="bottom", fontsize=8.5,
                      color="0.25")
        self.status = self.ax_pos.text(
            0.005, 0.985, "", transform=self.ax_pos.transAxes, ha="left", va="top",
            fontsize=8.2, color="0.25",
            bbox=dict(boxstyle="round,pad=0.25", fc="white", ec="0.8", alpha=0.9))
        self._update_status()

        self.ax_pos.callbacks.connect("xlim_changed", self._xlim_changed)
        self.fig.canvas.mpl_connect("key_press_event", self._key_press)
        self.fig.canvas.mpl_connect("scroll_event", self._scroll)
        self._autoscale_visible()

    def _visible_slice(self) -> slice:
        left, right = sorted(self.ax_pos.get_xlim())
        i0 = max(0, int(np.searchsorted(self.x, left, side="left")) - 1)
        i1 = min(self.x.size, int(np.searchsorted(self.x, right, side="right")) + 1)
        return slice(i0, max(i0 + 1, i1))

    @staticmethod
    def _set_limits(ax, values: np.ndarray, minimum_pad: float) -> None:
        finite = values[np.isfinite(values)]
        if not finite.size:
            return
        low, high = float(np.min(finite)), float(np.max(finite))
        span = high - low
        pad = max(minimum_pad, span * 0.08)
        if span == 0:
            pad = max(minimum_pad, abs(low) * 0.02)
        ax.set_ylim(low - pad, high + pad)

    def _autoscale_visible(self) -> None:
        view = self._visible_slice()
        self._set_limits(self.ax_pos, self.y[view], 2e-6)
        self._set_limits(self.ax_step, self.dy_um[view], 0.002)

    def _xlim_changed(self, _ax) -> None:
        if self.auto_y:
            self._autoscale_visible()
        self._update_status()
        self.fig.canvas.draw_idle()

    def _update_status(self) -> None:
        left, right = sorted(self.ax_pos.get_xlim())
        view = self._visible_slice()
        count = max(0, view.stop - view.start)
        width = max(0.0, right - left)
        self.status.set_text(
            f"visible: {left:.3f}â€“{right:.3f} s  ({width:.3f} s, {count:,} samples)\n"
            f"visible-window auto-Y: {'ON' if self.auto_y else 'OFF'}  |  "
            f"UINT32 rollovers unwrapped: {self.scope.rollover_count}")

    def _key_press(self, event) -> None:
        key = (event.key or "").lower()
        if key == "a":
            self.auto_y = not self.auto_y
            if self.auto_y:
                self._autoscale_visible()
        elif key == "y":
            self._autoscale_visible()
        elif key == "r":
            self.ax_pos.set_xlim(0, self.scope.time_s[-1])
            self._autoscale_visible()
        elif key == "c" and self.command_line is not None:
            self.command_visible = not self.command_visible
            self.command_line.set_visible(self.command_visible)
        elif key == "g":
            self.grid_visible = not self.grid_visible
            for ax in (self.ax_pos, self.ax_step):
                ax.grid(self.grid_visible, color="0.9", lw=0.55)
        else:
            return
        self._update_status()
        self.fig.canvas.draw_idle()

    def _scroll(self, event) -> None:
        if event.xdata is None or event.inaxes not in (self.ax_pos, self.ax_step):
            return
        left, right = self.ax_pos.get_xlim()
        factor = 0.75 if event.button == "up" else 1 / 0.75
        centre = event.xdata
        new_left = centre - (centre - left) * factor
        new_right = centre + (right - centre) * factor
        duration = self.scope.time_s[-1]
        if new_right - new_left >= duration:
            new_left, new_right = 0, duration
        else:
            if new_left < 0:
                new_right -= new_left
                new_left = 0
            if new_right > duration:
                new_left -= new_right - duration
                new_right = duration
        self.ax_pos.set_xlim(new_left, new_right)

    def show(self) -> None:
        plt.show()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=DEFAULT_SOURCE)
    parser.add_argument("--log", type=Path, default=DEFAULT_LOG)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--interactive", action="store_true",
                        help="open the full-resolution interactive viewer")
    args = parser.parse_args()

    scope = read_scope(args.source)
    controller = read_controller_log(args.log, scope.start_epoch_s)
    controller = align_controller_to_scope(scope, controller)
    if args.interactive:
        InteractiveViewer(scope, controller, args.source).show()
    else:
        render_static(scope, controller, args.source, args.output)
        print(f"Wrote {args.output}")
        print(f"Samples: {scope.index.size:,}; rate: {1000 / scope.sample_ms:.0f} Hz; "
              f"duration: {scope.time_s[-1]:.3f} s; rollovers: {scope.rollover_count}")


if __name__ == "__main__":
    main()

