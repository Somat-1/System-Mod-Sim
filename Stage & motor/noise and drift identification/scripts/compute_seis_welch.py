#!/usr/bin/env python3
'''Compute a gap-aware Welch PSD for a three-channel SEIS voltage recording.

Active intervals in analysis/exclusion_windows.csv are omitted before any
spectral estimate is formed. Periodic Hann windows overlap by 50 percent, are
mean-detrended individually, and never cross an excluded interval. Recordings
without an active exclusion are processed continuously in full.
'''

import argparse
import csv
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal.windows import hann

from render_seis_sequence import csv_blocks


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_SOURCE = ROOT / 'data' / '2026-10-15_seis' / 'TvalSeis.csv'
DEFAULT_SUMMARY = ROOT / 'analysis' / 'cache' / 'TvalSeis_sequence_50ms.npz'
DEFAULT_EXCLUSIONS = ROOT / 'analysis' / 'exclusion_windows.csv'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'
DEFAULT_RESULT = ROOT / 'analysis' / 'cache' / 'TvalSeis_welch_psd.npz'

COLOURS = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#161616'
MUTED = '#696762'
GRID = '#dfded9'
ALERT = '#bd3a32'


class WelchAccumulator:
    def __init__(self, fs, nperseg):
        self.fs = float(fs)
        self.nperseg = int(nperseg)
        self.hop = self.nperseg // 2
        self.window = hann(self.nperseg, sym=False).astype(np.float64)
        self.scale = self.fs * np.sum(self.window ** 2)
        self.total = np.zeros((self.nperseg // 2 + 1, 3), dtype=np.float64)
        self.count = 0
        self.buffer = np.empty((0, 3), dtype=np.float32)

    def reset_gap(self):
        '''Discard overlap state so a segment cannot bridge an excluded gap.'''
        self.buffer = np.empty((0, 3), dtype=np.float32)

    def feed(self, values):
        if not len(values):
            return
        self.buffer = (np.vstack((self.buffer, values))
                       if self.buffer.size else values.copy())
        while len(self.buffer) >= self.nperseg:
            segment = self.buffer[:self.nperseg].astype(np.float64)
            segment -= np.mean(segment, axis=0, keepdims=True)
            spectrum = np.fft.rfft(segment * self.window[:, None], axis=0)
            density = np.abs(spectrum) ** 2 / self.scale
            if self.nperseg % 2 == 0:
                density[1:-1] *= 2.0
            else:
                density[1:] *= 2.0
            self.total += density
            self.count += 1
            self.buffer = self.buffer[self.hop:]

    def finish(self):
        if self.count == 0:
            raise ValueError('no complete Welch windows remain')
        frequency = np.fft.rfftfreq(self.nperseg, d=1.0/self.fs)
        return frequency, self.total / self.count


def active_exclusions(path, source_name):
    windows = []
    with path.open(newline='', encoding='utf-8') as stream:
        for row in csv.DictReader(stream):
            if row['source'] == source_name and row['status'].lower() == 'active':
                windows.append((float(row['start_s']), float(row['end_s']),
                                row['reason']))
    return sorted(windows)


def kept_regions(duration_s, exclusions):
    regions = []
    cursor = 0.0
    for start, end, _ in exclusions:
        start = max(cursor, min(duration_s, start))
        end = max(start, min(duration_s, end))
        if start > cursor:
            regions.append((cursor, start))
        cursor = end
    if cursor < duration_s:
        regions.append((cursor, duration_s))
    return regions


def compute(source, summary_path, exclusion_path, nperseg):
    with np.load(summary_path) as summary:
        dt = float(summary['dt'])
        duration_s = float(summary['duration_s'])
        samples = int(summary['samples'])
    fs = 1.0 / dt
    exclusions = active_exclusions(exclusion_path, source.name)
    regions = kept_regions(duration_s, exclusions)
    index_regions = [(int(np.ceil(start * fs)), int(np.floor(end * fs)))
                     for start, end in regions]

    accumulator = WelchAccumulator(fs, nperseg)
    cursor = 0
    active_region = None
    for block_number, table in enumerate(csv_blocks(source), 1):
        table = np.atleast_2d(table)
        values = table[:, 1:4]
        block_start, block_end = cursor, cursor + len(values)
        for region_number, (start, end) in enumerate(index_regions):
            left, right = max(block_start, start), min(block_end, end)
            if left >= right:
                continue
            if active_region != region_number:
                accumulator.reset_gap()
                active_region = region_number
            accumulator.feed(values[left-block_start:right-block_start])
        cursor = block_end
        print(f'chunk {block_number}: {cursor:,}/{samples:,} samples; '
              f'{accumulator.count} Welch windows', flush=True)

    frequency, psd = accumulator.finish()
    kept_s = sum(end-start for start, end in regions)
    return dict(frequency_hz=frequency, psd_mV2_per_Hz=psd,
                fs_hz=np.array(fs), nperseg=np.array(nperseg),
                overlap_samples=np.array(nperseg//2),
                segment_count=np.array(accumulator.count),
                kept_duration_s=np.array(kept_s),
                excluded_duration_s=np.array(duration_s-kept_s),
                duration_s=np.array(duration_s)), exclusions, regions


def style():
    mpl.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'font.family': 'DejaVu Sans', 'font.size': 8.5,
        'axes.grid': True, 'axes.axisbelow': True,
        'grid.color': GRID, 'grid.linewidth': 0.55,
        'axes.edgecolor': GRID, 'axes.linewidth': 0.8,
        'xtick.color': MUTED, 'ytick.color': MUTED,
        'legend.frameon': False, 'pdf.fonttype': 42,
    })


def render(result, exclusions, regions, recording_label='SEIS voltage noise'):
    frequency = result['frequency_hz']
    psd = result['psd_mV2_per_Hz']
    fs = float(result['fs_hz'])
    nperseg = int(result['nperseg'])
    segments = int(result['segment_count'])
    resolution = fs / nperseg
    duration = float(result['duration_s'])

    fig = plt.figure(figsize=(15.5, 8.5))
    grid = fig.add_gridspec(2, 2, height_ratios=(0.20, 1.0),
                            width_ratios=(2.35, 1.0),
                            left=0.065, right=0.985, bottom=0.105, top=0.84,
                            wspace=0.20, hspace=0.34)
    timeline = fig.add_subplot(grid[0, :])
    focus = fig.add_subplot(grid[1, 0])
    full = fig.add_subplot(grid[1, 1])
    method = 'gap-aware' if exclusions else 'all samples admitted'
    fig.suptitle(f'{recording_label} — Welch power spectral density',
                 x=0.065, y=0.965, ha='left', fontsize=15,
                 fontweight='bold', color=INK)
    fig.text(0.065, 0.915,
             f'Periodic Hann window · {nperseg:,} samples ({nperseg/fs:.2f} s) '
             f'· 50% overlap · Δf = {resolution:.4f} Hz · {segments} windows '
             f'· {method}',
             ha='left', fontsize=8.5, color=MUTED)

    timeline.set_xlim(0, duration/60.0)
    timeline.set_ylim(0, 1)
    for start, end in regions:
        timeline.axvspan(start/60.0, end/60.0, color='#dce9f6', alpha=1.0,
                         linewidth=0)
    for start, end, _ in exclusions:
        timeline.axvspan(start/60.0, end/60.0, color=ALERT, alpha=0.22,
                         linewidth=0)
        timeline.text((start+end)/120.0, 0.5,
                      f'EXCLUDED\n{start/60:g}–{end/60:g} min',
                      ha='center', va='center', color=ALERT,
                      fontsize=8, fontweight='bold')
    if not exclusions:
        timeline.text(duration/120.0, 0.5, 'ALL SAMPLES ADMITTED\nno exclusions',
                      ha='center', va='center', color='#315f86',
                      fontsize=8, fontweight='bold')
    timeline.text(0.005, 0.82, '(a) Samples admitted to Welch estimate',
                  transform=timeline.transAxes, ha='left', va='top',
                  fontweight='bold')
    timeline.set_xlabel('elapsed recording time [min]')
    timeline.set_yticks([])

    keep = frequency > 0
    for channel, colour in enumerate(COLOURS):
        label = f'channel {channel+1}'
        focus.loglog(frequency[keep], psd[keep, channel], color=colour,
                     linewidth=0.85, label=label)
        full.loglog(frequency[keep], psd[keep, channel], color=colour,
                    linewidth=0.72, label=label)
    focus.set_xlim(max(resolution, 0.01), min(400.0, fs/2))
    focus.set_xlabel('frequency [Hz]')
    focus.set_ylabel('PSD [mV²/Hz]')
    focus.set_title('(b) Primary view: frequencies up to 400 Hz', loc='left',
                    fontweight='bold')
    focus.legend(loc='best', ncol=3)
    full.set_xlim(max(resolution, 0.01), fs/2)
    full.set_xlabel('frequency [Hz]')
    full.set_ylabel('PSD [mV²/Hz]')
    full.set_title('(c) Full measured band', loc='left', fontweight='bold')
    full.legend(loc='best')

    for axis in (timeline, focus, full):
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
    exclusion_note = (
        'No segment crosses an active exclusion.'
        if exclusions else 'No time interval is excluded from this recording.'
    )
    fig.text(0.065, 0.035,
             f'Each segment is mean-detrended. {exclusion_note} Voltage has not '
             'been converted to acceleration because no SEIS sensitivity '
             'calibration was supplied.',
             ha='left', fontsize=7.5, color=MUTED)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--summary', type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument('--exclusions', type=Path, default=DEFAULT_EXCLUSIONS)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--result', type=Path, default=DEFAULT_RESULT)
    parser.add_argument('--nperseg', type=int, default=2**20)
    parser.add_argument('--render-only', action='store_true',
                        help='render an existing result without rereading raw data')
    parser.add_argument('--no-save-result', action='store_true',
                        help='render the computed PSD without writing an NPZ')
    parser.add_argument('--output-stem', default='seis_03_welch_psd')
    parser.add_argument('--title', default='SEIS voltage noise')
    parser.add_argument('--dpi', type=int, default=180)
    parser.add_argument('--pdf', action='store_true',
                        help='also write a PDF copy')
    args = parser.parse_args()

    style()
    exclusions = active_exclusions(args.exclusions, args.source.name)
    if args.render_only:
        with np.load(args.result, allow_pickle=False) as saved:
            result = {key: saved[key] for key in saved.files}
        regions = kept_regions(float(result['duration_s']), exclusions)
    else:
        result, exclusions, regions = compute(
            args.source, args.summary, args.exclusions, args.nperseg)
        if not args.no_save_result:
            args.result.parent.mkdir(parents=True, exist_ok=True)
            np.savez_compressed(
                args.result, **result,
                exclusion_start_s=np.array([x[0] for x in exclusions]),
                exclusion_end_s=np.array([x[1] for x in exclusions]))

    args.out.mkdir(parents=True, exist_ok=True)
    fig = render(result, exclusions, regions, recording_label=args.title)
    png = args.out / f'{args.output_stem}.png'
    pdf = args.out / f'{args.output_stem}.pdf'
    fig.savefig(png, dpi=args.dpi)
    if args.pdf:
        fig.savefig(pdf)
    plt.close(fig)
    if not args.render_only and not args.no_save_result:
        print(f'wrote {args.result}')
    print(f'wrote {png}')
    if args.pdf:
        print(f'wrote {pdf}')


if __name__ == '__main__':
    main()
