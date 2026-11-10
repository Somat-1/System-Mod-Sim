#!/usr/bin/env python3
'''Render the chirp accelerometer samples without signal processing.

The 20 kHz captures contain many more samples than a raster image has pixels.
For display only, each horizontal bin is therefore drawn as the exact minimum
to maximum of the raw samples that fall in that bin.  Nothing is averaged,
filtered, normalised, detrended, calibrated, or converted to the frequency
domain.
'''

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import vizstyle as vs


HERE = Path(__file__).resolve().parent
AXIS_NAMES = ('axis 1', 'axis 2', 'axis 3')


def minmax_envelope(values, first, last, bins=6000):
    '''Return bin centres and exact raw minima/maxima for a sample interval.'''
    values = values[first:last]
    count = values.shape[0]
    bins = min(bins, count)
    starts = np.linspace(0, count, bins + 1, dtype=np.int64)
    # linspace is strictly increasing because bins <= count.
    lo = np.minimum.reduceat(values, starts[:-1], axis=0)
    hi = np.maximum.reduceat(values, starts[:-1], axis=0)
    centres = (starts[:-1] + starts[1:] - 1) / 2.0 + first
    return centres, lo, hi


def draw_page(name, accel, dt, units, first, last, interval_label):
    '''Build one three-axis raw time-history page.'''
    indices, lo, hi = minmax_envelope(accel, first, last)
    time_s = indices * dt
    start_s, end_s = first * dt, (last - 1) * dt

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(11.69, 8.27))
    fig.subplots_adjust(left=0.09, right=0.985, top=0.87, bottom=0.105,
                        hspace=0.18)
    fig.suptitle(f'Raw chirp measurement — {name} — {interval_label}',
                 x=0.09, ha='left', fontsize=16, fontweight='bold')
    fig.text(0.09, 0.91,
             f'{last - first:,} original samples · {1 / dt:,.2f} samples/s · '
             f'{start_s:.3f}–{end_s:.3f} s in source recording',
             ha='left', color=vs.INK_SOFT, fontsize=9)
    fig.text(0.09, 0.035,
             'Unprocessed acceleration. Each vertical stroke is the exact raw '
             'minimum-to-maximum range in one display bin; no samples are '
             'averaged or filtered.',
             ha='left', color=vs.INK_MUTED, fontsize=8)

    colours = (vs.SERIES[0], vs.SERIES[1], vs.SERIES[2])
    for axis, ax in enumerate(axes):
        ax.vlines(time_s, lo[:, axis], hi[:, axis], color=colours[axis],
                  linewidth=0.45, alpha=0.92, rasterized=True)
        ax.axhline(0, color=vs.INK_MUTED, linewidth=0.65, zorder=0)
        ymin = float(np.min(lo[:, axis]))
        ymax = float(np.max(hi[:, axis]))
        pad = max((ymax - ymin) * 0.055, 1e-6)
        ax.set_ylim(ymin - pad, ymax + pad)
        ax.set_xlim(start_s, end_s)
        ax.set_ylabel(f'{AXIS_NAMES[axis]}\n[{units[axis]}]')
        ax.text(0.997, 0.94, f'raw range {ymin:.3g} to {ymax:.3g} {units[axis]}',
                transform=ax.transAxes, ha='right', va='top', fontsize=7.5,
                color=vs.INK_MUTED)
        vs.despine(ax)
    axes[-1].set_xlabel('time in source recording [s]')
    return fig


def recording_intervals(result, sample_count, dt):
    '''Yield the complete recording followed by each identified chirp run.'''
    yield 0, sample_count, 'complete recording', 'full'
    for run_number, run in enumerate(result.get('runs', []), 1):
        legs = run.get('legs', [])
        if not legs:
            continue
        first = max(0, int(np.floor(min(leg['t0'] for leg in legs) / dt)))
        last = min(sample_count,
                   int(np.ceil(max(leg['t1'] for leg in legs) / dt)) + 1)
        profile = str(run.get('profile', 'chirp')).upper()
        yield first, last, f'run {run_number} ({profile} sweep)', f'run_{run_number}'


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cache', type=Path, default=HERE / 'cache')
    parser.add_argument('--results', type=Path,
                        default=HERE / 'chirp_results.json')
    parser.add_argument('--out', type=Path,
                        default=HERE / 'plots' / 'chirp_raw')
    args = parser.parse_args()

    vs.use()
    results = json.loads(args.results.read_text(encoding='utf-8'))
    by_name = {item['recording']: item
               for item in results.get('recordings', [])}
    args.out.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out / 'chirp_raw_measurements.pdf'

    written = []
    with PdfPages(pdf_path) as pdf:
        for name, result in by_name.items():
            cache_path = args.cache / f'{name}.npz'
            if not cache_path.exists():
                print(f'{name}: cache missing, skipped')
                continue
            with np.load(cache_path) as cache:
                accel = cache['accel']
                dt = float(cache['dt'])
                units = [str(unit) for unit in cache['units']]
                for first, last, label, suffix in recording_intervals(
                        result, accel.shape[0], dt):
                    fig = draw_page(name, accel, dt, units, first, last, label)
                    path = args.out / f'chirp_raw_{name}_{suffix}.png'
                    fig.savefig(path, dpi=200)
                    pdf.savefig(fig, dpi=200)
                    plt.close(fig)
                    written.append(path)
                    print(f'wrote {path.name}')

    print(f'wrote {pdf_path.name} ({len(written)} pages)')


if __name__ == '__main__':
    main()
