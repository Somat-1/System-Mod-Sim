#!/usr/bin/env python3
'''Render direct FFT periodograms without Welch segmentation or averaging.

Each continuous admitted interval is mean-detrended and transformed once with
a rectangular window. Intervals separated by an active exclusion are kept as
independent traces so no artificial sample adjacency is introduced.
'''

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.fft import rfft, rfftfreq

from compute_seis_welch import active_exclusions, kept_regions
from render_seis_sequence import csv_blocks


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_SOURCE = ROOT / 'data' / '2026-10-15_seis' / 'TvalSeis.csv'
DEFAULT_SUMMARY = ROOT / 'analysis' / 'cache' / 'TvalSeis_sequence_50ms.npz'
DEFAULT_EXCLUSIONS = ROOT / 'analysis' / 'exclusion_windows.csv'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'

COLOURS = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#161616'
MUTED = '#696762'
GRID = '#dfded9'
ALERT = '#bd3a32'
REGION_STYLES = ('-', '--', ':', '-.')


def style():
    mpl.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'font.family': 'DejaVu Sans', 'font.size': 8.5,
        'axes.grid': True, 'axes.axisbelow': True,
        'grid.color': GRID, 'grid.linewidth': 0.55,
        'axes.edgecolor': GRID, 'axes.linewidth': 0.8,
        'xtick.color': MUTED, 'ytick.color': MUTED,
        'legend.frameon': False,
    })


def load_record_regions(source, fs, regions):
    '''Stream selected continuous intervals into float32 channel arrays.'''
    index_regions = [(int(np.ceil(start * fs)), int(np.floor(end * fs)))
                     for start, end in regions]
    arrays = [np.empty((end-start, 3), dtype=np.float32)
              for start, end in index_regions]
    filled = [0 for _ in arrays]
    cursor = 0

    for block_number, table in enumerate(csv_blocks(source), 1):
        table = np.atleast_2d(table)
        values = table[:, 1:4]
        block_start, block_end = cursor, cursor + len(values)
        for number, (start, end) in enumerate(index_regions):
            left, right = max(block_start, start), min(block_end, end)
            if left >= right:
                continue
            count = right-left
            arrays[number][filled[number]:filled[number]+count] = (
                values[left-block_start:right-block_start])
            filled[number] += count
        cursor = block_end
        print(f'load chunk {block_number}: {cursor:,} samples', flush=True)

    for number, (array, expected) in enumerate(zip(arrays, index_regions)):
        target = expected[1]-expected[0]
        if filled[number] != target:
            raise RuntimeError(
                f'region {number} incomplete: {filled[number]:,}/{target:,}')
    return arrays


def log_bin_mean(frequency, density, edges):
    '''Arithmetic mean PSD in display-only logarithmic frequency bins.'''
    keep = frequency > 0
    index = np.searchsorted(edges, frequency[keep], side='right') - 1
    valid = (index >= 0) & (index < len(edges)-1)
    index = index[valid]
    values = density[keep][valid]
    count = np.bincount(index, minlength=len(edges)-1)
    total = np.bincount(index, weights=values, minlength=len(edges)-1)
    result = np.full(len(edges)-1, np.nan, dtype=np.float64)
    occupied = count > 0
    result[occupied] = total[occupied] / count[occupied]
    return result


def periodogram_region(values, fs, focus_max_hz=400.0, full_bins=16000):
    '''Calculate one rectangular-window periodogram for a continuous record.'''
    samples = len(values)
    frequency = rfftfreq(samples, d=1.0/fs)
    focus = (frequency >= max(frequency[1], 0.01)) & (
        frequency <= min(focus_max_hz, fs/2))
    edges = np.geomspace(max(frequency[1], 0.01), fs/2, full_bins+1)
    full_frequency = np.sqrt(edges[:-1] * edges[1:]).astype(np.float32)
    full_density = np.empty((full_bins, 3), dtype=np.float32)
    dominant = []

    for channel in range(3):
        signal = values[:, channel]
        signal -= np.mean(signal, dtype=np.float64)
        spectrum = rfft(signal, workers=-1)
        density = ((spectrum.real*spectrum.real + spectrum.imag*spectrum.imag)
                   / (fs*samples))
        if samples % 2 == 0:
            density[1:-1] *= 2.0
        else:
            density[1:] *= 2.0
        full_density[:, channel] = log_bin_mean(frequency, density, edges)
        peak = int(np.argmax(density[focus]))
        dominant.append((float(frequency[focus][peak]),
                         float(density[focus][peak])))
        del spectrum, density
        print(f'  channel {channel+1}: direct FFT complete', flush=True)

    return dict(full_frequency_hz=full_frequency,
                full_psd_mV2_per_Hz=full_density,
                samples=samples,
                duration_s=samples/fs,
                dominant=dominant)


def render(results, exclusions, regions, duration_s, fs, title):
    fig = plt.figure(figsize=(15.5, 8.5))
    grid = fig.add_gridspec(2, 2, height_ratios=(0.20, 1.0),
                            width_ratios=(2.35, 1.0),
                            left=0.065, right=0.985, bottom=0.11, top=0.84,
                            wspace=0.20, hspace=0.34)
    timeline = fig.add_subplot(grid[0, :])
    focus = fig.add_subplot(grid[1, 0])
    full = fig.add_subplot(grid[1, 1])

    fig.suptitle(f'{title} — direct FFT periodogram',
                 x=0.065, y=0.965, ha='left', fontsize=15,
                 fontweight='bold', color=INK)
    region_text = (f'{len(regions)} independent continuous records'
                   if len(regions) > 1 else 'one continuous record')
    fig.text(
        0.065, 0.915,
        f'Rectangular window · {region_text} · no overlap · no Welch averaging '
        f'· sample rate {fs:,.1f} Hz',
        ha='left', fontsize=8.5, color=MUTED,
    )

    timeline.set_xlim(0, duration_s/60.0)
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
        timeline.text(duration_s/120.0, 0.5,
                      'ALL SAMPLES ADMITTED\nno exclusions',
                      ha='center', va='center', color='#315f86',
                      fontsize=8, fontweight='bold')
    timeline.text(0.005, 0.82, '(a) Continuous records transformed',
                  transform=timeline.transAxes, ha='left', va='top',
                  fontweight='bold')
    timeline.set_xlabel('elapsed recording time [min]')
    timeline.set_yticks([])

    for number, result in enumerate(results):
        linestyle = REGION_STYLES[number % len(REGION_STYLES)]
        alpha = 0.68 if number == 0 else 0.42
        for channel, colour in enumerate(COLOURS):
            label = f'channel {channel+1}' if number == 0 else None
            focus.loglog(result['full_frequency_hz'],
                         result['full_psd_mV2_per_Hz'][:, channel],
                         color=colour, linewidth=0.62, alpha=alpha,
                         linestyle=linestyle, label=label)
            full.loglog(result['full_frequency_hz'],
                        result['full_psd_mV2_per_Hz'][:, channel],
                        color=colour, linewidth=0.72, alpha=alpha,
                        linestyle=linestyle, label=label)

    focus.set_xlim(0.01, min(400.0, fs/2))
    focus.set_xlabel('frequency [Hz]')
    focus.set_ylabel('PSD [mV²/Hz]')
    focus.set_title('(b) Direct periodogram: frequencies up to 400 Hz',
                    loc='left', fontweight='bold')
    channel_legend = focus.legend(loc='lower left', ncol=3)
    if len(regions) > 1:
        region_handles = []
        for number, (start, end) in enumerate(regions):
            region_handles.append(Line2D(
                [0], [0], color=MUTED,
                linestyle=REGION_STYLES[number % len(REGION_STYLES)],
                linewidth=1.1,
                label=f'{start/60:.2f}–{end/60:.2f} min'))
        focus.add_artist(channel_legend)
        focus.legend(handles=region_handles, loc='upper right',
                     title='continuous record')

    full.set_xlim(0.01, fs/2)
    full.set_xlabel('frequency [Hz]')
    full.set_ylabel('PSD [mV²/Hz]')
    full.set_title('(c) Full measured band', loc='left', fontweight='bold')
    full.legend(loc='best')

    for axis in (timeline, focus, full):
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
    fig.text(
        0.065, 0.035,
        'Each continuous interval is mean-detrended once. Both panels show '
        'the direct periodogram averaged only into narrow logarithmic frequency '
        'display bins. No time-domain segmentation or averaging is applied.',
        ha='left', fontsize=7.5, color=MUTED,
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--summary', type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument('--exclusions', type=Path, default=DEFAULT_EXCLUSIONS)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--output-stem', default='seis_05_direct_periodogram')
    parser.add_argument('--title', default='SEIS table measurement')
    parser.add_argument('--dpi', type=int, default=150)
    args = parser.parse_args()

    style()
    with np.load(args.summary, allow_pickle=False) as summary:
        fs = 1.0/float(summary['dt'])
        duration_s = float(summary['duration_s'])
    exclusions = active_exclusions(args.exclusions, args.source.name)
    regions = kept_regions(duration_s, exclusions)
    arrays = load_record_regions(args.source, fs, regions)

    results = []
    for number, values in enumerate(arrays):
        print(f'periodogram region {number+1}/{len(arrays)}: '
              f'{len(values):,} samples', flush=True)
        results.append(periodogram_region(values, fs))
        arrays[number] = None

    for number, result in enumerate(results):
        for channel, (frequency, density) in enumerate(result['dominant']):
            print(f'region {number+1}, channel {channel+1}: dominant '
                  f'{frequency:.6f} Hz, {density:.6g} mV^2/Hz')

    args.out.mkdir(parents=True, exist_ok=True)
    fig = render(results, exclusions, regions, duration_s, fs, args.title)
    png = args.out / f'{args.output_stem}.png'
    fig.savefig(png, dpi=args.dpi)
    plt.close(fig)
    print(f'wrote {png}')


if __name__ == '__main__':
    main()
