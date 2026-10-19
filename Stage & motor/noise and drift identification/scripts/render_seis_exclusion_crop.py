#!/usr/bin/env python3
'''Render the raw-data context around an active SEIS exclusion.

The raw CSV is not edited. This plot reads the extrema-preserving 50 ms
display cache and leaves active exclusion intervals blank, making the samples
admitted to later analyses unambiguous. It also produces a two-page reader
containing the existing Welch PSD followed by the raw-data crop.
'''

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from compute_seis_welch import (
    ALERT,
    COLOURS,
    INK,
    MUTED,
    active_exclusions,
    kept_regions,
    render as render_welch,
    style,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_SOURCE = ROOT / 'data' / '2026-10-15_seis' / 'TvalSeis.csv'
DEFAULT_SUMMARY = ROOT / 'analysis' / 'cache' / 'TvalSeis_sequence_50ms.npz'
DEFAULT_RESULT = ROOT / 'analysis' / 'cache' / 'TvalSeis_welch_psd.npz'
DEFAULT_EXCLUSIONS = ROOT / 'analysis' / 'exclusion_windows.csv'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'


def load_npz(path):
    with np.load(path, allow_pickle=False) as saved:
        return {key: saved[key] for key in saved.files}


def admitted_mask(time_s, exclusions):
    admitted = np.ones(time_s.shape, dtype=bool)
    for start, end, _ in exclusions:
        admitted &= ~((time_s >= start) & (time_s < end))
    return admitted


def render_crop(summary, exclusions, context_s=120.0):
    if not exclusions:
        raise ValueError('no active exclusion interval found')

    first_start = min(item[0] for item in exclusions)
    last_end = max(item[1] for item in exclusions)
    duration = float(summary['duration_s'])
    crop_start = max(0.0, first_start - context_s)
    crop_end = min(duration, last_end + context_s)

    time_s = summary['centre_s']
    visible = (time_s >= crop_start) & (time_s <= crop_end)
    admitted = admitted_mask(time_s, exclusions)
    show = visible & admitted
    time_min = time_s / 60.0

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(11.69, 8.27))
    fig.subplots_adjust(left=0.09, right=0.985, top=0.86, bottom=0.105,
                        hspace=0.19)
    fig.suptitle('SEIS raw-data crop - active disturbance exclusion',
                 x=0.09, y=0.965, ha='left', fontsize=15,
                 fontweight='bold', color=INK)
    fig.text(
        0.09,
        0.915,
        f'Displayed interval {crop_start/60:.1f}-{crop_end/60:.1f} min | '
        'raw min-max envelope in 50 ms bins | excluded samples are not plotted',
        ha='left', fontsize=8.5, color=MUTED,
    )

    for channel, axis in enumerate(axes[:3]):
        minimum = np.where(show, summary['minimum'][:, channel], np.nan)
        maximum = np.where(show, summary['maximum'][:, channel], np.nan)
        mean = np.where(show, summary['mean'][:, channel], np.nan)
        axis.fill_between(time_min, minimum, maximum,
                          where=np.isfinite(minimum), color=COLOURS[channel],
                          alpha=0.78, linewidth=0)
        axis.plot(time_min, mean, color=COLOURS[channel], linewidth=0.45)
        axis.set_ylabel(f'channel {channel + 1}\n[mV]')
        axis.set_title(f'({chr(97 + channel)}) Raw channel {channel + 1}',
                       loc='left', fontweight='bold')

    peak = np.where(show, summary['magnitude_peak'], np.nan)
    rms = np.where(show, summary['magnitude_rms'], np.nan)
    axes[3].semilogy(time_min, peak, color=ALERT, linewidth=0.75,
                     label='peak per 50 ms bin')
    axes[3].semilogy(time_min, rms, color='#555555', linewidth=0.65,
                     label='RMS per 50 ms bin')
    axes[3].set_ylabel('relative magnitude\n[mV]')
    axes[3].set_xlabel('elapsed recording time [min]')
    axes[3].set_title('(d) Three-channel disturbance magnitude', loc='left',
                      fontweight='bold')
    axes[3].legend(loc='upper right')

    for start, end, _ in exclusions:
        if end < crop_start or start > crop_end:
            continue
        for axis in axes:
            axis.axvspan(start / 60.0, end / 60.0, color=ALERT,
                         alpha=0.12, linewidth=0)
            axis.axvline(start / 60.0, color=ALERT, linewidth=0.7,
                         linestyle='--')
            axis.axvline(end / 60.0, color=ALERT, linewidth=0.7,
                         linestyle='--')
        axes[0].text((start + end) / 120.0, 0.5,
                     'EXCLUDED FROM ANALYSIS\n26-30 min',
                     ha='center', va='center', color=ALERT,
                     fontsize=8, fontweight='bold',
                     transform=axes[0].get_xaxis_transform())

    for axis in axes:
        axis.set_xlim(crop_start / 60.0, crop_end / 60.0)
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)

    fig.text(
        0.09,
        0.035,
        'The blank band is the active 1560-1800 s mask. Raw samples remain '
        'unchanged on disk; only their admission to analysis and this view is '
        'suppressed.',
        ha='left', fontsize=7.5, color=MUTED,
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--summary', type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument('--result', type=Path, default=DEFAULT_RESULT)
    parser.add_argument('--exclusions', type=Path, default=DEFAULT_EXCLUSIONS)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--context-s', type=float, default=120.0)
    parser.add_argument('--pdf', action='store_true',
                        help='also write standalone and combined PDFs')
    args = parser.parse_args()

    style()
    summary = load_npz(args.summary)
    exclusions = active_exclusions(args.exclusions, args.source.name)
    crop = render_crop(summary, exclusions, context_s=args.context_s)

    args.out.mkdir(parents=True, exist_ok=True)
    crop_png = args.out / 'seis_04_raw_exclusion_crop.png'
    crop_pdf = args.out / 'seis_04_raw_exclusion_crop.pdf'
    reader_pdf = args.out / 'seis_welch_with_raw_crop.pdf'
    crop.savefig(crop_png, dpi=240)
    if args.pdf:
        result = load_npz(args.result)
        regions = kept_regions(float(result['duration_s']), exclusions)
        crop.savefig(crop_pdf)
        welch = render_welch(result, exclusions, regions)
        with PdfPages(reader_pdf) as pdf:
            pdf.savefig(welch)
            pdf.savefig(crop)
        plt.close(welch)
    plt.close(crop)

    print(f'wrote {crop_png}')
    if args.pdf:
        print(f'wrote {crop_pdf}')
        print(f'wrote {reader_pdf}')


if __name__ == '__main__':
    main()
