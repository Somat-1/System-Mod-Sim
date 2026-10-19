#!/usr/bin/env python3
'''Render the complete SEIS sequence and nominate the dominant disturbance.

The 20 kHz recording is summarized for display by exact per-bin minima and
maxima, together with RMS and peak vector magnitude relative to the first
10 seconds. No filtering or calibration is applied. The source CSV is never
modified; a proposed exclusion is written as metadata with status=candidate.
'''

import argparse
import csv
import io
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_SOURCE = ROOT / 'data' / '2026-10-15_seis' / 'TvalSeis.csv'
DEFAULT_CACHE = ROOT / 'analysis' / 'cache' / 'TvalSeis_sequence_50ms.npz'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'
DEFAULT_EXCLUSIONS = ROOT / 'analysis' / 'exclusion_windows.csv'

COLOURS = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#161616'
MUTED = '#696762'
GRID = '#dfded9'
ALERT = '#bd3a32'
READ_BLOCK = 64 * 1024 * 1024
CACHE_VERSION = 2


def csv_blocks(path):
    '''Yield complete numeric row blocks without splitting a CSV line.'''
    with path.open('rb') as stream:
        header = stream.readline().decode('ascii').strip()
        if header != 's,mV,mV,mV':
            raise ValueError(f'unexpected SEIS header: {header!r}')
        tail = b''
        while True:
            chunk = stream.read(READ_BLOCK)
            if not chunk:
                break
            chunk = tail + chunk
            cut = chunk.rfind(b'\n') + 1
            tail, complete = chunk[cut:], chunk[:cut]
            if complete:
                yield np.loadtxt(io.BytesIO(complete), delimiter=',',
                                 dtype=np.float32)
        if tail.strip():
            yield np.loadtxt(io.BytesIO(tail), delimiter=',',
                             dtype=np.float32)


def summarize(path, bin_ms=50.0, baseline_s=10.0):
    '''Stream the multi-gigabyte CSV into a compact, extrema-preserving cache.'''
    minima, maxima, means, magnitude_peak, magnitude_rms = [], [], [], [], []
    carry = np.empty((0, 3), dtype=np.float32)
    baseline = None
    dt = None
    total = 0
    first_time = None
    last_time = None

    for number, table in enumerate(csv_blocks(path), 1):
        table = np.atleast_2d(table)
        values = table[:, 1:4]
        time = table[:, 0].astype(np.float64, copy=False)
        if first_time is None:
            first_time = float(time[0])
            dt = float(np.median(np.diff(time[:min(20_000, len(time))])))
            take = min(len(values), max(1, int(round(baseline_s / dt))))
            baseline = np.median(values[:take], axis=0).astype(np.float32)
        last_time = float(time[-1])
        total += len(values)
        bin_samples = max(1, int(round(bin_ms / 1000.0 / dt)))
        combined = np.vstack((carry, values)) if carry.size else values
        complete = (len(combined) // bin_samples) * bin_samples
        if complete:
            blocks = combined[:complete].reshape(-1, bin_samples, 3)
            minima.append(np.min(blocks, axis=1))
            maxima.append(np.max(blocks, axis=1))
            means.append(np.mean(blocks, axis=1))
            magnitude = np.sqrt(np.sum((blocks - baseline) ** 2, axis=2))
            magnitude_peak.append(np.max(magnitude, axis=1))
            magnitude_rms.append(np.sqrt(np.mean(magnitude ** 2, axis=1)))
        carry = combined[complete:].copy()
        print(f'chunk {number}: {total:,} samples', flush=True)

    if carry.size:
        minima.append(np.min(carry, axis=0, keepdims=True))
        maxima.append(np.max(carry, axis=0, keepdims=True))
        means.append(np.mean(carry, axis=0, keepdims=True))
        magnitude = np.sqrt(np.sum((carry - baseline) ** 2, axis=1))
        magnitude_peak.append(np.array([np.max(magnitude)], dtype=np.float32))
        magnitude_rms.append(np.array([
            np.sqrt(np.mean(magnitude ** 2))], dtype=np.float32))

    minimum = np.vstack(minima).astype(np.float32)
    maximum = np.vstack(maxima).astype(np.float32)
    mean = np.vstack(means).astype(np.float32)
    peak = np.concatenate(magnitude_peak).astype(np.float32)
    rms = np.concatenate(magnitude_rms).astype(np.float32)
    # Recover the timebase from the complete capture. This avoids accumulating
    # the small timestamp round-off introduced by the float32 block reader.
    duration_s = last_time - first_time
    dt = duration_s / (total - 1)
    bin_s = bin_samples * dt
    centre_s = (np.arange(len(peak), dtype=np.float64) + 0.5) * bin_s
    centre_s[-1] = min(centre_s[-1], last_time)
    return dict(centre_s=centre_s, minimum=minimum, maximum=maximum,
                mean=mean, magnitude_peak=peak, magnitude_rms=rms,
                baseline=baseline, dt=np.array(dt), bin_s=np.array(bin_s),
                samples=np.array(total), duration_s=np.array(duration_s),
                cache_version=np.array(CACHE_VERSION),
                source_size=np.array(path.stat().st_size),
                source_mtime_ns=np.array(path.stat().st_mtime_ns))


def load_or_build(source, cache, force=False):
    valid = False
    if cache.exists() and not force:
        with np.load(cache) as saved:
            valid = ('cache_version' in saved.files and
                     int(saved['cache_version']) == CACHE_VERSION and
                     int(saved['source_size']) == source.stat().st_size and
                     int(saved['source_mtime_ns']) == source.stat().st_mtime_ns)
        if valid:
            print(f'loading cache {cache}')
            with np.load(cache) as saved:
                return {key: saved[key] for key in saved.files}
    print(f'summarizing {source} (first pass may take several minutes)')
    result = summarize(source)
    cache.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(cache, **result)
    print(f'wrote {cache}')
    return result


def proposed_window(summary, half_width_s=10.0):
    '''Nominate the largest raw vector-magnitude excursion for confirmation.'''
    peak = summary['magnitude_peak']
    index = int(np.nanargmax(peak))
    centre = float(summary['centre_s'][index])
    duration = float(summary['duration_s'])
    return max(0.0, centre-half_width_s), min(duration, centre+half_width_s), centre


def write_candidate(path, source, start, end, peak_time, summary):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w', newline='', encoding='utf-8') as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            'source', 'start_s', 'end_s', 'status', 'reason',
            'detected_peak_s', 'detected_peak_mV'])
        writer.writeheader()
        writer.writerow(dict(source=source.name, start_s=f'{start:.3f}',
                             end_s=f'{end:.3f}', status='candidate',
                             reason='operator entered room and touched table',
                             detected_peak_s=f'{peak_time:.3f}',
                             detected_peak_mV=f'{np.max(summary["magnitude_peak"]):.6g}'))


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


def shade_candidate(axes, start, end):
    for axis in axes:
        axis.axvspan(start, end, color=ALERT, alpha=0.12, linewidth=0)
        axis.axvline(start, color=ALERT, linewidth=0.65, linestyle='--')
        axis.axvline(end, color=ALERT, linewidth=0.65, linestyle='--')


def overview(summary, start, end, peak_time):
    time_min = summary['centre_s'] / 60.0
    start_min, end_min = start / 60.0, end / 60.0
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(11.69, 8.27))
    fig.subplots_adjust(left=0.09, right=0.985, top=0.875, bottom=0.09,
                        hspace=0.18)
    fig.suptitle('SEIS noise-and-drift recording — raw sequence overview',
                 x=0.09, ha='left', fontsize=15, fontweight='bold')
    fig.text(0.09, 0.915,
             f'{int(summary["samples"]):,} samples · '
             f'{1 / float(summary["dt"]):,.1f} samples/s · '
             f'{float(summary["duration_s"]) / 60:.2f} min · units as recorded',
             ha='left', fontsize=8.5, color=MUTED)
    for channel, axis in enumerate(axes[:3]):
        axis.fill_between(time_min, summary['minimum'][:, channel],
                          summary['maximum'][:, channel],
                          color=COLOURS[channel], alpha=0.78, linewidth=0)
        axis.plot(time_min, summary['mean'][:, channel],
                  color=COLOURS[channel], linewidth=0.35, alpha=0.8)
        axis.set_ylabel(f'channel {channel+1}\n[mV]')
        axis.set_title(f'({chr(97+channel)}) Raw channel {channel+1}', loc='left',
                       fontweight='bold')
    axes[3].semilogy(time_min, summary['magnitude_peak'], color=ALERT,
                     linewidth=0.65, label='peak per 50 ms display bin')
    axes[3].semilogy(time_min, summary['magnitude_rms'], color='#555555',
                     linewidth=0.55, label='RMS per display bin')
    axes[3].set_ylabel('magnitude relative\nto initial baseline [mV]')
    axes[3].set_xlabel('elapsed recording time [min]')
    axes[3].set_title('(d) Three-channel disturbance magnitude', loc='left',
                      fontweight='bold')
    axes[3].legend(loc='upper right')
    shade_candidate(axes, start_min, end_min)
    axes[0].text(0.995, 0.94,
                 f'candidate table-touch exclusion\n{start:.3f}–{end:.3f} s '
                 f'(peak {peak_time:.3f} s)',
                 transform=axes[0].transAxes, ha='right', va='top',
                 color=ALERT, fontsize=7.5)
    for axis in axes:
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
    fig.text(0.09, 0.025,
             'Exact raw min–max envelope in 50 ms display bins. Candidate is '
             'highlighted only; no samples have been removed or corrected.',
             ha='left', color=MUTED, fontsize=7.5)
    return fig


def detail(summary, start, end, peak_time, context_s=20.0):
    time = summary['centre_s']
    keep = (time >= peak_time-context_s) & (time <= peak_time+context_s)
    local_t = time[keep] - peak_time
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(11.69, 8.27))
    fig.subplots_adjust(left=0.09, right=0.985, top=0.875, bottom=0.09,
                        hspace=0.18)
    fig.suptitle('Dominant SEIS transient — confirmation view', x=0.09,
                 ha='left', fontsize=15, fontweight='bold')
    fig.text(0.09, 0.915,
             f'Candidate peak at {peak_time:.3f} s · proposed mask '
             f'{start:.3f}–{end:.3f} s · status: CANDIDATE',
             ha='left', fontsize=8.5, color=ALERT)
    for channel, axis in enumerate(axes[:3]):
        axis.fill_between(local_t, summary['minimum'][keep, channel],
                          summary['maximum'][keep, channel],
                          color=COLOURS[channel], alpha=0.78, linewidth=0)
        axis.plot(local_t, summary['mean'][keep, channel],
                  color=COLOURS[channel], linewidth=0.45)
        axis.set_ylabel(f'channel {channel+1}\n[mV]')
        axis.set_title(f'({chr(97+channel)}) Raw channel {channel+1}', loc='left',
                       fontweight='bold')
    axes[3].semilogy(local_t, summary['magnitude_peak'][keep], color=ALERT,
                     linewidth=0.9, label='peak magnitude')
    axes[3].semilogy(local_t, summary['magnitude_rms'][keep], color='#555555',
                     linewidth=0.7, label='RMS magnitude')
    axes[3].set_ylabel('relative magnitude [mV]')
    axes[3].set_xlabel('time relative to detected peak [s]')
    axes[3].set_title('(d) Three-channel disturbance magnitude', loc='left',
                      fontweight='bold')
    axes[3].legend(loc='upper right')
    rel_start, rel_end = start-peak_time, end-peak_time
    shade_candidate(axes, rel_start, rel_end)
    for axis in axes:
        axis.axvline(0, color=ALERT, linewidth=0.9)
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
    fig.text(0.09, 0.025,
             'Use this page to confirm or adjust the proposed exclusion before '
             'any noise, drift, PSD, or stability calculation is run.',
             ha='left', color=MUTED, fontsize=7.5)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--cache', type=Path, default=DEFAULT_CACHE)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--exclusions', type=Path, default=DEFAULT_EXCLUSIONS)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--pdf', action='store_true',
                        help='also write the multi-page PDF report')
    args = parser.parse_args()

    style()
    summary = load_or_build(args.source, args.cache, force=args.force)
    start, end, peak_time = proposed_window(summary)
    write_candidate(args.exclusions, args.source, start, end, peak_time, summary)
    args.out.mkdir(parents=True, exist_ok=True)

    figures = [overview(summary, start, end, peak_time),
               detail(summary, start, end, peak_time)]
    names = ('seis_01_recorded_sequence.png',
             'seis_02_disturbance_candidate.png')
    for fig, name in zip(figures, names):
        fig.savefig(args.out / name, dpi=220)
        print(f'wrote {args.out / name}')
    if args.pdf:
        pdf_path = args.out / 'seis_recorded_sequence.pdf'
        with PdfPages(pdf_path) as pdf:
            for fig in figures:
                pdf.savefig(fig)
        print(f'wrote {pdf_path}')
    for fig in figures:
        plt.close(fig)
    print(f'wrote {args.exclusions}')
    print(json.dumps(dict(candidate_start_s=start, candidate_end_s=end,
                          detected_peak_s=peak_time,
                          detected_peak_mV=float(np.max(summary['magnitude_peak']))),
                     indent=2))


if __name__ == '__main__':
    main()
