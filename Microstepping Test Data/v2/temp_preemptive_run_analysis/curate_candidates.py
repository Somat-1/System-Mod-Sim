#!/usr/bin/env python3
"""Score preemptive-run cuts and render contact sheets for visual review."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("preemptive", HERE / "analyze_preemptive_run.py")
preemptive = importlib.util.module_from_spec(spec)
spec.loader.exec_module(preemptive)


def moving_average(values, width=51):
    if values.size < width:
        return values.copy()
    kernel = np.ones(width) / width
    return np.convolve(values, kernel, mode="same")


def normalized(values):
    values = values - np.median(values)
    scale = np.percentile(values, 95) - np.percentile(values, 5)
    return values / scale if scale > 0 else values


def evaluate(segment, ids_t, ids_y, cmd_t, cmd_y, offset):
    a, b = segment["esp_start_s"], segment["esp_end_s"]
    grid = np.linspace(a, b, 1200)
    y = np.interp(grid + offset, ids_t, ids_y)
    x = np.interp(grid, cmd_t, cmd_y)
    ys = moving_average(y, 31)
    # Ignore the convolution margins and remove an affine drift from measurement.
    valid = slice(40, -40)
    tt = grid[valid]
    trend = np.polyval(np.polyfit(tt, ys[valid], 1), tt)
    yd = ys[valid] - trend
    xd = x[valid] - np.mean(x[valid])
    position_corr = abs(float(np.corrcoef(xd, yd)[0, 1])) if np.std(xd) and np.std(yd) else 0.0
    dx = np.diff(x[valid])
    dy = np.diff(ys[valid])
    derivative_corr = abs(float(np.corrcoef(dx, dy)[0, 1])) if np.std(dx) and np.std(dy) else 0.0
    measured_span = float(np.percentile(ys[valid], 95) - np.percentile(ys[valid], 5))
    noise = float(np.median(np.abs(np.diff(y)))) * 1.4826
    snr = measured_span / max(noise, 1e-9)
    jumps = np.abs(np.diff(y))
    artifact_ratio = float(np.max(jumps) / max(measured_span, 1e-9))
    return {
        "position_corr": position_corr,
        "derivative_corr": derivative_corr,
        "measured_span_counts": measured_span,
        "sample_noise_counts": noise,
        "span_to_noise": snr,
        "artifact_ratio": artifact_ratio,
        "grid": grid, "y": y, "x": x,
    }


def main():
    ids_t, ids_y = preemptive.load_ids()
    rows = preemptive.load_esp()
    esp_t0, _, cmd_t, cmd_y = preemptive.command_trace(rows)
    offset, *_ = preemptive.find_alignment(ids_t, ids_y, cmd_t, cmd_y)
    segments = preemptive.build_segments(rows, esp_t0, cmd_t[-1])
    records = []
    raw = {}
    for segment in segments:
        key = (segment["run_index"], segment["block"], segment["condition"])
        metrics = evaluate(segment, ids_t, ids_y, cmd_t, cmd_y, offset)
        raw[key] = metrics
        records.append({**segment, **{k: v for k, v in metrics.items() if not isinstance(v, np.ndarray)}})

    by_condition = {}
    for record in records:
        key = (record["block"], record["condition"])
        by_condition.setdefault(key, []).append(record)
    for record in records:
        peers = [p for p in by_condition[(record["block"], record["condition"])]
                 if p["run_index"] in (1, 2)]
        repeatability = np.nan
        if len(peers) == 2:
            a = normalized(raw[(peers[0]["run_index"], peers[0]["block"], peers[0]["condition"])]["y"])
            b = normalized(raw[(peers[1]["run_index"], peers[1]["block"], peers[1]["condition"])]["y"])
            repeatability = float(np.corrcoef(a, b)[0, 1])
        record["cross_current_repeatability"] = repeatability
        repeat_score = max(0.0, repeatability) if np.isfinite(repeatability) else 0.0
        artifact_penalty = max(0.0, record["artifact_ratio"] - 2.0) * 0.05
        record["quality_score"] = (
            0.42 * record["position_corr"] + 0.18 * record["derivative_corr"]
            + 0.40 * repeat_score - artifact_penalty
        )

    records.sort(key=lambda r: r["quality_score"], reverse=True)
    fields = [k for k in records[0] if k not in {"grid", "y", "x"}]
    with (HERE / "candidate_quality.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader(); writer.writerows(records)

    # One representative per condition, ranked; the contact sheets are for
    # explicit human screening before any model comparison is accepted.
    representatives = []
    used = set()
    for record in records:
        key = (record["block"], record["condition"])
        if key not in used:
            used.add(key); representatives.append(record)
    review_dir = HERE / "plots" / "curation_review"
    review_dir.mkdir(parents=True, exist_ok=True)
    for page, start in enumerate(range(0, len(representatives), 12), 1):
        subset = representatives[start:start + 12]
        fig, axes = plt.subplots(4, 3, figsize=(17, 13), constrained_layout=True)
        for ax, record in zip(axes.flat, subset):
            data = raw[(record["run_index"], record["block"], record["condition"])]
            t = data["grid"] - data["grid"][0]
            ax.plot(t, normalized(data["y"]), color="#1769aa", lw=.75, label="IDS normalized")
            ax.plot(t, normalized(data["x"]), color="#d94801", lw=1.0, label="command normalized")
            ax.set_title(f"R{record['run_index']} {record['block']} {record['condition']}\n"
                         f"Q={record['quality_score']:.2f}, repeat={record['cross_current_repeatability']:.2f}, "
                         f"artifact={record['artifact_ratio']:.1f}", fontsize=9)
            ax.grid(alpha=.2)
        for ax in axes.flat[len(subset):]: ax.axis("off")
        axes.flat[0].legend(fontsize=7)
        fig.suptitle("Candidate screening — normalized shapes (ranking is advisory)")
        fig.savefig(review_dir / f"candidate_review_{page:02d}.png", dpi=150)
        plt.close(fig)
    print("Top candidates:")
    for r in representatives[:15]:
        print(f"{r['block']:14s} {r['condition']:18s} run={r['run_index']} "
              f"Q={r['quality_score']:.3f} repeat={r['cross_current_repeatability']:.3f} "
              f"artifact={r['artifact_ratio']:.2f}")


if __name__ == "__main__":
    main()
