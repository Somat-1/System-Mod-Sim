#!/usr/bin/env python3
"""Align the pre-emptive IDS capture to the ESP command log and plot segments."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
V2 = HERE.parent
IDS_PATH = V2 / "data" / "PreemptiveRundata.csv"
ESP_PATH = V2 / "data" / "esp32_runs" / "full_campaign_marked_complete_20260828.txt"
PLOTS = HERE / "plots"
SEGMENTS = PLOTS / "segments"


def load_ids() -> tuple[np.ndarray, np.ndarray]:
    indices: list[int] = []
    values: list[int] = []
    with IDS_PATH.open(encoding="utf-8-sig", errors="replace") as handle:
        for line in handle:
            fields = line.rstrip().split("\t")
            if len(fields) < 2 or not fields[0].isdigit() or not fields[1].isdigit():
                continue
            indices.append(int(fields[0]))
            raw = int(fields[1])
            values.append(raw - 2**32 if raw >= 2**31 else raw)
    t = np.asarray(indices, dtype=float) * 0.001
    y = np.asarray(values, dtype=float)
    return t, y


def load_esp() -> list[dict[str, str]]:
    header = None
    rows: list[dict[str, str]] = []
    with ESP_PATH.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            text = line.strip().strip("\r")
            if text.startswith("timestamp_us,event,"):
                header = next(csv.reader([text]))
                continue
            if header is None or not text or not text[0].isdigit():
                continue
            fields = next(csv.reader([text]))
            if len(fields) == len(header):
                rows.append(dict(zip(header, fields)))
    return rows


def command_trace(rows: list[dict[str, str]], dt: float = 0.02):
    campaign = next(r for r in rows if r["event"] == "CAMPAIGN_START")
    stop = next(r for r in reversed(rows) if r["event"] == "SAFE_STOP")
    t0 = int(campaign["timestamp_us"]) / 1e6
    t1 = int(stop["timestamp_us"]) / 1e6
    grid = np.arange(0, t1 - t0 + dt, dt)
    position = np.zeros_like(grid)
    batches = [r for r in rows if r["event"] == "PULSE_BATCH" and int(r["timestamp_us"]) / 1e6 >= t0]
    knots_t = [0.0]
    knots_p = [0.0]
    prior = 0.0
    for row in batches:
        start = int(row["timestamp_us"]) / 1e6 - t0
        done = int(row["rmt_done_us"]) / 1e6 - t0
        after = float(row["position_u16"])
        knots_t.extend((max(start, knots_t[-1]), max(done, start + 1e-6)))
        knots_p.extend((prior, after))
        prior = after
    position[:] = np.interp(grid, np.asarray(knots_t), np.asarray(knots_p))
    return t0, t1, grid, position


def alignment_score(ids_t, ids_y, cmd_t, cmd_y, offset):
    sample_t = np.arange(0, cmd_t[-1], 0.2)
    x = np.interp(sample_t, cmd_t, cmd_y)
    y = np.interp(sample_t + offset, ids_t, ids_y)
    # First differences emphasize moves and suppress slow encoder drift.
    dx = np.diff(x)
    dy = np.diff(y)
    dx -= dx.mean()
    dy -= dy.mean()
    denom = np.linalg.norm(dx) * np.linalg.norm(dy)
    return 0.0 if denom == 0 else abs(float(np.dot(dx, dy) / denom))


def find_alignment(ids_t, ids_y, cmd_t, cmd_y):
    max_offset = max(0.0, ids_t[-1] - cmd_t[-1])
    coarse = np.arange(0.0, max_offset + 0.1001, 0.1)
    scores = np.asarray([alignment_score(ids_t, ids_y, cmd_t, cmd_y, x) for x in coarse])
    best = coarse[int(np.argmax(scores))]
    fine = np.arange(max(0, best - 0.2), min(max_offset, best + 0.2) + 0.0051, 0.005)
    fine_scores = np.asarray([alignment_score(ids_t, ids_y, cmd_t, cmd_y, x) for x in fine])
    offset = float(fine[int(np.argmax(fine_scores))])
    y = np.interp(cmd_t + offset, ids_t, ids_y)
    design = np.column_stack((cmd_y, np.ones_like(cmd_y)))
    scale, intercept = np.linalg.lstsq(design, y, rcond=None)[0]
    fitted = scale * cmd_y + intercept
    alignment_corr = float(np.max(fine_scores))
    return offset, float(scale), float(intercept), alignment_corr, coarse, scores


def row_time(row, esp_t0):
    return int(row["timestamp_us"]) / 1e6 - esp_t0


def build_segments(rows, esp_t0, stop_rel):
    segments = []
    measured_blocks = {"BLOCK_0_START", "A1", "A2", "B", "E", "BLOCK_0_END"}
    for run in sorted({int(r["run_index"]) for r in rows if r["run_index"].isdigit() and int(r["run_index"]) > 0}):
        run_rows = [r for r in rows if int(r["run_index"] or 0) == run]
        configs = [r for r in run_rows if r["event"] == "RUN_CONFIG"]
        if not configs:
            continue
        config = configs[0]
        for block in measured_blocks:
            block_rows = [r for r in run_rows if r["block"] == block]
            if not block_rows:
                continue
            if block.startswith("BLOCK_0"):
                groups = [(block, block_rows)]
            else:
                labels = []
                for r in block_rows:
                    if r["event"] == "PULSE_BATCH" and r["label"] not in labels:
                        labels.append(r["label"])
                groups = [(label, [r for r in block_rows if r["label"] == label]) for label in labels]
            for label, group in groups:
                times = [row_time(r, esp_t0) for r in group]
                start = max(0.0, min(times) - 0.15)
                end = min(stop_rel, max(times) + 0.55)
                complete = not (end >= stop_rel - 0.6)
                segments.append({
                    "run_index": run, "current": config["current"], "mres": int(config["mres"]),
                    "block": block, "condition": label, "esp_start_s": start,
                    "esp_end_s": end, "complete": complete,
                })
    return segments


def plot_segment(segment, ids_t, ids_y, cmd_t, cmd_y, offset, scale, intercept):
    a, b = segment["esp_start_s"], segment["esp_end_s"]
    pad = 0.25
    mask = (ids_t >= offset + a - pad) & (ids_t <= offset + b + pad)
    cmask = (cmd_t >= a - pad) & (cmd_t <= b + pad)
    fig, ax = plt.subplots(figsize=(11, 4.8), constrained_layout=True)
    ax.plot(ids_t[mask] - offset, ids_y[mask], color="#1769aa", lw=0.7, label="IDS encoder")
    command_ax = ax.twinx()
    command_ax.plot(cmd_t[cmask], cmd_y[cmask] / 16.0, color="#d94801", lw=1.2,
                    label="ESP commanded position")
    ax.axvspan(a, b, color="#808080", alpha=0.08)
    status = "complete" if segment["complete"] else "partial (aborted)"
    ax.set_title(f"Run {segment['run_index']:02d} | {segment['current']} | 1/{segment['mres']} | "
                 f"{segment['block']} | {segment['condition']} | {status}")
    ax.set_xlabel("Time since ESP CAMPAIGN_START [s]")
    ax.set_ylabel("Encoder counts")
    command_ax.set_ylabel("Commanded position [full steps]", color="#d94801")
    ax.grid(alpha=0.22)
    lines = ax.lines + command_ax.lines
    ax.legend(lines, [line.get_label() for line in lines], loc="best")
    name = (f"run_{segment['run_index']:02d}_{segment['current']}_mres_{segment['mres']}_"
            f"{segment['block']}_{segment['condition']}").replace("/", "_")
    path = SEGMENTS / f"{name}.png"
    fig.savefig(path, dpi=160)
    plt.close(fig)
    return path


def main():
    PLOTS.mkdir(parents=True, exist_ok=True)
    SEGMENTS.mkdir(parents=True, exist_ok=True)
    ids_t, ids_y = load_ids()
    rows = load_esp()
    esp_t0, esp_t1, cmd_t, cmd_y = command_trace(rows)
    offset, scale, intercept, corr, search_x, search_y = find_alignment(ids_t, ids_y, cmd_t, cmd_y)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 9), sharex=True, constrained_layout=True)
    stride = 20
    ax1.plot(ids_t[::stride], ids_y[::stride], color="#1769aa", lw=0.55)
    ax1.axvspan(offset, offset + cmd_t[-1], color="#fdae6b", alpha=0.18, label="aligned ESP run")
    ax1.set_ylabel("IDS encoder [signed counts]")
    ax1.set_title("Preemptive IDS capture — complete recording")
    ax1.grid(alpha=0.2)
    ax1.legend()
    ax2.plot(ids_t[::stride], ids_y[::stride], color="#1769aa", lw=0.5, label="IDS encoder")
    command_ax = ax2.twinx()
    command_ax.plot(cmd_t + offset, cmd_y / 16.0, color="#d94801", lw=0.8,
                    label="ESP commanded position")
    ax2.set_xlabel("Time since IDS recording start [s]")
    ax2.set_ylabel("Encoder counts")
    command_ax.set_ylabel("Commanded position [full steps]", color="#d94801")
    ax2.grid(alpha=0.2)
    lines = ax2.lines + command_ax.lines
    ax2.legend(lines, [line.get_label() for line in lines], loc="best")
    fig.savefig(PLOTS / "00_full_run_overview.png", dpi=170)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(10, 4), constrained_layout=True)
    ax.plot(search_x, search_y, color="#54278f")
    ax.axvline(offset, color="#d94801", ls="--", label=f"best offset {offset:.3f} s")
    ax.set(xlabel="IDS time of ESP CAMPAIGN_START [s]", ylabel="|derivative correlation|",
           title="Clock-alignment search")
    ax.grid(alpha=0.25)
    ax.legend()
    fig.savefig(PLOTS / "01_alignment_diagnostic.png", dpi=170)
    plt.close(fig)

    segments = build_segments(rows, esp_t0, cmd_t[-1])
    for segment in segments:
        segment["ids_start_s"] = offset + segment["esp_start_s"]
        segment["ids_end_s"] = offset + segment["esp_end_s"]
        segment["plot"] = str(plot_segment(segment, ids_t, ids_y, cmd_t, cmd_y,
                                            offset, scale, intercept).relative_to(HERE))
    with (HERE / "segment_manifest.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(segments[0]))
        writer.writeheader()
        writer.writerows(segments)

    complete = sum(bool(x["complete"]) for x in segments)
    summary = f"""# Preemptive run visualization

- IDS source: `{IDS_PATH.relative_to(V2)}` ({len(ids_y):,} samples at 1 kHz; {ids_t[-1]:.3f} s)
- ESP source: `{ESP_PATH.relative_to(V2)}`
- ESP campaign duration captured: {cmd_t[-1]:.3f} s
- Estimated ESP `CAMPAIGN_START` in IDS time: {offset:.3f} s
- Alignment derivative correlation: {corr:.4f}
- Fitted encoder counts per commanded 1/16-full-step unit: {scale:.6g}
- Extracted labeled segments: {len(segments)} ({complete} complete, {len(segments)-complete} partial)
- Run ended by commanded abort during run 3, A1, N2 negative; later tests are absent, not inferred.

The command overlays use a separate full-step axis; no encoder calibration is assumed.
Segment boundaries come from the ESP timestamps shifted by the marker-derived clock offset. `segment_manifest.csv`
contains both ESP-relative and IDS-relative boundaries and should be treated as the
machine-readable splice index.
"""
    (HERE / "README.md").write_text(summary, encoding="utf-8")
    print(summary)


if __name__ == "__main__":
    main()
