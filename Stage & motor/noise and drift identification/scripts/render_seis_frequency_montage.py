#!/usr/bin/env python3
'''Render a clean table-versus-floor frequency-domain montage.

The left column is a direct single-record FFT periodogram expressed as ASD.
The right column is the gap-aware Welch PSD. Faint coloured lines retain the
three channels; the heavy black line is a smoothed geometric-mean trend used
only as a visual guide.
'''

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from scipy.fft import rfft, rfftfreq
from scipy.ndimage import gaussian_filter1d

from compute_seis_welch import (
    WelchAccumulator,
    active_exclusions,
    kept_regions,
)
from render_seis_sequence import csv_blocks


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
TABLE_SOURCE = ROOT / 'data' / '2026-10-15_seis' / 'TvalSeis.csv'
TABLE_SUMMARY = ROOT / 'analysis' / 'cache' / 'TvalSeis_sequence_50ms.npz'
FLOOR_SOURCE = ROOT / 'data' / '2026-10-15_seis_floor' / 'TvalSeisfloor.csv'
FLOOR_SUMMARY = ROOT / 'analysis' / 'cache' / 'TvalSeisfloor_sequence_50ms.npz'
EXCLUSIONS = ROOT / 'analysis' / 'exclusion_windows.csv'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'

COLOURS = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#171717'
MUTED = '#696762'
GRID = '#dfded9'
ALERT = '#bd3a32'


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


def binned_mean(frequency, density, edges):
    keep = (frequency >= edges[0]) & (frequency <= edges[-1])
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


def smooth_trace(values, sigma=22.0):
    '''Smooth one positive trace on an equally log-spaced frequency grid.'''
    trace = np.asarray(values, dtype=np.float64)
    valid = np.isfinite(trace) & (trace > 0)
    index = np.arange(len(trace))
    filled = np.interp(index, index[valid], np.log(trace[valid]))
    return np.exp(gaussian_filter1d(filled, sigma=sigma, mode='nearest'))


def smooth_channels(values, sigma=22.0):
    return np.column_stack([
        smooth_trace(values[:, channel], sigma=sigma)
        for channel in range(values.shape[1])
    ])


def analyse_recording(source, summary_path, exclusion_path, edges, nperseg):
    with np.load(summary_path, allow_pickle=False) as summary:
        fs = 1.0/float(summary['dt'])
        duration_s = float(summary['duration_s'])
        samples = int(summary['samples'])
    exclusions = active_exclusions(exclusion_path, source.name)
    regions = kept_regions(duration_s, exclusions)
    index_regions = [(int(np.ceil(start*fs)), int(np.floor(end*fs)))
                     for start, end in regions]
    longest_number = int(np.argmax([end-start for start, end in index_regions]))
    direct_start, direct_end = index_regions[longest_number]
    direct = np.empty((direct_end-direct_start, 3), dtype=np.float32)
    direct_filled = 0

    welch = WelchAccumulator(fs, nperseg)
    active_region = None
    cursor = 0
    for block_number, table in enumerate(csv_blocks(source), 1):
        table = np.atleast_2d(table)
        values = table[:, 1:4]
        block_start, block_end = cursor, cursor+len(values)

        left = max(block_start, direct_start)
        right = min(block_end, direct_end)
        if left < right:
            count = right-left
            direct[direct_filled:direct_filled+count] = (
                values[left-block_start:right-block_start])
            direct_filled += count

        for region_number, (start, end) in enumerate(index_regions):
            left = max(block_start, start)
            right = min(block_end, end)
            if left >= right:
                continue
            if active_region != region_number:
                welch.reset_gap()
                active_region = region_number
            welch.feed(values[left-block_start:right-block_start])
        cursor = block_end
        print(f'{source.name}: chunk {block_number}, {cursor:,}/{samples:,} '
              f'samples, {welch.count} Welch windows', flush=True)

    if direct_filled != len(direct):
        raise RuntimeError(
            f'direct record incomplete: {direct_filled:,}/{len(direct):,}')

    welch_frequency, welch_psd = welch.finish()
    welch_binned = np.column_stack([
        binned_mean(welch_frequency, welch_psd[:, channel], edges)
        for channel in range(3)
    ])
    del welch_psd

    direct_frequency = rfftfreq(len(direct), d=1.0/fs)
    direct_asd = np.empty((len(edges)-1, 3), dtype=np.float64)
    for channel in range(3):
        signal = direct[:, channel]
        signal -= np.mean(signal, dtype=np.float64)
        spectrum = rfft(signal, workers=-1)
        density = ((spectrum.real*spectrum.real + spectrum.imag*spectrum.imag)
                   / (fs*len(signal)))
        if len(signal) % 2 == 0:
            density[1:-1] *= 2.0
        else:
            density[1:] *= 2.0
        direct_asd[:, channel] = np.sqrt(
            binned_mean(direct_frequency, density, edges))
        del spectrum, density
        print(f'{source.name}: direct FFT channel {channel+1}/3', flush=True)
    del direct

    return dict(
        fs=fs,
        duration_s=duration_s,
        admitted_s=sum(end-start for start, end in regions),
        exclusions=exclusions,
        regions=regions,
        direct_region=regions[longest_number],
        direct_asd=direct_asd,
        direct_trend=smooth_channels(direct_asd),
        welch_psd=welch_binned,
        welch_trend=smooth_channels(welch_binned),
        welch_windows=welch.count,
    )


def plot_panel(axis, frequency, estimate, trend, title, ylabel, panel, colour,
               ylim):
    axis.loglog(frequency, estimate, color=colour,
                linewidth=0.55, alpha=0.25)
    axis.loglog(frequency, trend, color=INK, linewidth=1.65,
                label='smoothed per-axis trace')
    axis.axvspan(400, frequency[-1], color='#777777', alpha=0.035,
                 linewidth=0)
    axis.axvline(400, color=MUTED, linewidth=0.65, linestyle='--')
    axis.axvline(50, color='#8b4d4d', linewidth=0.6, linestyle=':')
    axis.set_xlim(0.02, frequency[-1])
    axis.set_ylim(*ylim)
    axis.set_xlabel('frequency [Hz]')
    axis.set_ylabel(ylabel)
    axis.set_title(f'({panel}) {title}', loc='left', fontweight='bold')
    axis.spines['top'].set_visible(False)
    axis.spines['right'].set_visible(False)


def render(table, floor, frequency):
    fig, axes = plt.subplots(3, 4, figsize=(19.0, 12.0), sharex=True)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.82, bottom=0.09,
                        wspace=0.25, hspace=0.30)
    fig.suptitle('SEIS frequency-domain montage — measurements and axes separated',
                 x=0.075, y=0.972, ha='left', fontsize=16,
                 fontweight='bold', color=INK)
    fig.text(
        0.075, 0.938,
        'One panel per recorded axis · faint colour: direct estimate · heavy '
        'black: smoothed trace of that same axis · full band shown · '
        'dashed marker: 400 Hz · '
        'dotted marker: 50 Hz',
        ha='left', fontsize=8.5, color=MUTED,
    )
    fig.text(0.292, 0.875, 'SEIS / TABLE RECORDING — TvalSeis.csv',
             ha='center', fontsize=11, fontweight='bold', color='#315f86')
    fig.text(0.752, 0.875, 'FLOOR RECORDING — TvalSeisfloor.csv',
             ha='center', fontsize=11, fontweight='bold', color='#7a513a')
    fig.text(0.18, 0.845, 'DIRECT FFT ASD', ha='center', fontsize=8.5,
             fontweight='bold', color=MUTED)
    fig.text(0.407, 0.845, 'WELCH PSD', ha='center', fontsize=8.5,
             fontweight='bold', color=MUTED)
    fig.text(0.64, 0.845, 'DIRECT FFT ASD', ha='center', fontsize=8.5,
             fontweight='bold', color=MUTED)
    fig.text(0.87, 0.845, 'WELCH PSD', ha='center', fontsize=8.5,
             fontweight='bold', color=MUTED)

    letters = iter('abcdefghijkl')
    direct_ylim = (1e-5, 10.0)
    psd_ylim = (3e-10, 30.0)
    for channel, colour in enumerate(COLOURS):
        plot_panel(
            axes[channel, 0], frequency, table['direct_asd'][:, channel],
            table['direct_trend'][:, channel], 'Direct FFT ASD',
            'ASD [mV/√Hz]', next(letters), colour, direct_ylim)
        plot_panel(
            axes[channel, 1], frequency, table['welch_psd'][:, channel],
            table['welch_trend'][:, channel], 'Welch PSD',
            'PSD [mV²/Hz]', next(letters), colour, psd_ylim)
        plot_panel(
            axes[channel, 2], frequency, floor['direct_asd'][:, channel],
            floor['direct_trend'][:, channel], 'Direct FFT ASD',
            'ASD [mV/√Hz]', next(letters), colour, direct_ylim)
        plot_panel(
            axes[channel, 3], frequency, floor['welch_psd'][:, channel],
            floor['welch_trend'][:, channel], 'Welch PSD',
            'PSD [mV²/Hz]', next(letters), colour, psd_ylim)
        row_y = (0.69, 0.455, 0.22)[channel]
        fig.text(0.018, row_y, f'AXIS / CHANNEL {channel+1}', rotation=90,
                 va='center', ha='center', fontsize=9.5, fontweight='bold',
                 color=colour)

    legend_handles = [
        Line2D([0], [0], color='#777777', linewidth=1.1, alpha=0.4,
               label='axis estimate'),
        Line2D([0], [0], color=INK, linewidth=1.8,
               label='smoothed per-axis trace'),
    ]
    fig.legend(handles=legend_handles, loc='upper right',
               bbox_to_anchor=(0.985, 0.938), ncol=2)

    table_start, table_end = table['direct_region']
    floor_start, floor_end = floor['direct_region']
    fig.text(
        0.292, 0.835,
        f'FFT: {table_start/60:.2f}–{table_end/60:.2f} min uninterrupted · '
        f'Welch: {table["welch_windows"]} windows, 26–30 min excluded',
        ha='center', fontsize=7.2, color=MUTED)
    fig.text(
        0.752, 0.835,
        f'FFT: {floor_start/60:.2f}–{floor_end/60:.2f} min · '
        f'Welch: {floor["welch_windows"]} windows · no exclusion',
        ha='center', fontsize=7.2, color=MUTED)

    fig.text(
        0.075, 0.025,
        'The left column is an FFT-derived amplitude spectral density, not a '
        'true FRF: no input/reference channel was recorded. Log-frequency '
        'binning is used for display only. Each black curve is derived only '
        'from the axis shown in that panel; axes and recordings are not '
        'combined. Voltage is not acceleration-calibrated.',
        ha='left', fontsize=7.5, color=MUTED,
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--dpi', type=int, default=150)
    parser.add_argument('--bins', type=int, default=12000)
    parser.add_argument('--nperseg', type=int, default=2**20)
    args = parser.parse_args()

    style()
    edges = np.geomspace(0.01, 10_000.0, args.bins+1)
    frequency = np.sqrt(edges[:-1]*edges[1:])
    table = analyse_recording(TABLE_SOURCE, TABLE_SUMMARY, EXCLUSIONS,
                              edges, args.nperseg)
    floor = analyse_recording(FLOOR_SOURCE, FLOOR_SUMMARY, EXCLUSIONS,
                              edges, args.nperseg)
    fig = render(table, floor, frequency)
    args.out.mkdir(parents=True, exist_ok=True)
    output = args.out / 'seis_frequency_montage.png'
    fig.savefig(output, dpi=args.dpi)
    plt.close(fig)
    print(f'wrote {output}')


if __name__ == '__main__':
    main()
