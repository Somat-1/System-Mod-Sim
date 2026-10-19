#!/usr/bin/env python3
'''Build an extrema-preserving cache and render the full floor recording.'''

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

from render_seis_sequence import load_or_build


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_SOURCE = ROOT / 'data' / '2026-10-15_seis_floor' / 'TvalSeisfloor.csv'
DEFAULT_CACHE = ROOT / 'analysis' / 'cache' / 'TvalSeisfloor_sequence_50ms.npz'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'

COLOURS = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#161616'
MUTED = '#696762'
GRID = '#dfded9'


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


def render(summary):
    time_min = summary['centre_s'] / 60.0
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(15.5, 8.5))
    fig.subplots_adjust(left=0.07, right=0.985, top=0.86, bottom=0.105,
                        hspace=0.19)
    fig.suptitle('SEIS floor measurement — raw sequence overview',
                 x=0.07, y=0.965, ha='left', fontsize=15,
                 fontweight='bold', color=INK)
    fig.text(
        0.07, 0.915,
        f'{int(summary["samples"]):,} samples · '
        f'{1 / float(summary["dt"]):,.1f} samples/s · '
        f'{float(summary["duration_s"]) / 60:.2f} min · units as recorded',
        ha='left', fontsize=8.5, color=MUTED,
    )

    for channel, axis in enumerate(axes[:3]):
        axis.fill_between(time_min, summary['minimum'][:, channel],
                          summary['maximum'][:, channel],
                          color=COLOURS[channel], alpha=0.78, linewidth=0,
                          rasterized=True)
        axis.plot(time_min, summary['mean'][:, channel],
                  color=COLOURS[channel], linewidth=0.35, alpha=0.85)
        axis.set_ylabel(f'channel {channel + 1}\n[mV]')
        axis.set_title(f'({chr(97 + channel)}) Raw channel {channel + 1}',
                       loc='left', fontweight='bold')

    axes[3].semilogy(time_min, summary['magnitude_peak'], color='#bd3a32',
                     linewidth=0.65, label='peak per 50 ms bin')
    axes[3].semilogy(time_min, summary['magnitude_rms'], color='#555555',
                     linewidth=0.55, label='RMS per 50 ms bin')
    axes[3].set_ylabel('relative magnitude\n[mV]')
    axes[3].set_xlabel('elapsed recording time [min]')
    axes[3].set_title('(d) Three-channel disturbance magnitude', loc='left',
                      fontweight='bold')
    axes[3].legend(loc='upper right')

    for axis in axes:
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)

    fig.text(
        0.07, 0.035,
        'Exact raw min–max envelope in 50 ms display bins; centre line is the '
        'bin mean. No samples or time intervals are excluded.',
        ha='left', fontsize=7.5, color=MUTED,
    )
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source', type=Path, default=DEFAULT_SOURCE)
    parser.add_argument('--cache', type=Path, default=DEFAULT_CACHE)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--dpi', type=int, default=180)
    parser.add_argument('--pdf', action='store_true',
                        help='also write a PDF copy')
    args = parser.parse_args()

    style()
    summary = load_or_build(args.source, args.cache, force=args.force)
    args.out.mkdir(parents=True, exist_ok=True)
    fig = render(summary)
    png = args.out / 'seis_floor_01_recorded_sequence.png'
    pdf = args.out / 'seis_floor_01_recorded_sequence.pdf'
    fig.savefig(png, dpi=args.dpi)
    if args.pdf:
        fig.savefig(pdf)
    plt.close(fig)
    print(f'wrote {png}')
    if args.pdf:
        print(f'wrote {pdf}')


if __name__ == '__main__':
    main()
