#!/usr/bin/env python3
'''One noise figure per recording: the table and the floor, side by side in
style but each standing on its own.

Each figure has two panels:

  (a) the spectrum itself, with the instrument noise floor marked, so it is
      visible where the recording stops measuring the ground and starts
      measuring the sensor;
  (b) cumulative RMS integrated from low frequency upward, which answers the
      question a PSD does not: how much of the total motion comes from
      below any given frequency.

Panel (b) is the useful one for deciding what to fix. A PSD peak that carries
no area is not worth chasing; the steps in the cumulative curve are.
'''

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from compute_seis_transmissibility import load, octave_smooth

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_OUT = ROOT / 'analysis' / 'plots'
CACHE = ROOT / 'analysis' / 'cache'

AXIS_COLOUR = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#161616'
MUTED = '#6b6963'
GRID = '#e2e1dc'
WARM = '#bd3a32'

BANDS = ((0.2, 1), (1, 10), (10, 100), (100, 400))
RECORDINGS = (
    ('table', 'TvalSeis_welch_10s.npz',
     'Table (SEIS) — accelerometer on the optical table beside the stage'),
    ('floor', 'TvalSeisfloor_welch_10s.npz',
     'Floor — same accelerometer on the laboratory floor'),
)


def style():
    mpl.rcParams.update({
        'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
        'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
        'axes.linewidth': 0.8, 'axes.spines.top': False,
        'axes.spines.right': False, 'figure.facecolor': 'white'})


def noise_floor(frequency, psd, top=(2000.0, 9000.0), margin_db=3.0):
    '''The level the spectrum settles to at high frequency, and where it
    first comes within 3 dB of it.'''
    window = (frequency >= top[0]) & (frequency <= top[1])
    level = float(np.median(psd[window]))
    limit = level * 10.0 ** (margin_db / 10.0)
    reached = psd <= limit
    onset = float(frequency[np.argmax(reached)]) if reached.any() else np.nan
    return level, onset


def band_rms(frequency, psd):
    '''RMS in each band, as sqrt of the integrated one-sided PSD.'''
    out = []
    for lo, hi in BANDS:
        mask = (frequency >= lo) & (frequency <= hi)
        out.append(float(np.sqrt(np.trapezoid(psd[mask], frequency[mask]))))
    return out


def figure(record, data, fraction, out):
    frequency = data['frequency']
    fig, (ax, ax2) = plt.subplots(2, 1, figsize=(8.0, 6.6), sharex=True,
                                  constrained_layout=True)
    show = (frequency >= 0.2) & (frequency <= 9000)

    onsets, levels, table_rows = [], [], []
    for index in range(3):
        smooth = octave_smooth(frequency, data['psd'][:, index], fraction)
        level, onset = noise_floor(frequency, smooth)
        onsets.append(onset)
        levels.append(level)
        ax.loglog(frequency[show], smooth[show], color=AXIS_COLOUR[index],
                  lw=1.3, label=f'axis {index + 1}')
        ax.axhline(level, color=AXIS_COLOUR[index], lw=0.7, ls=':')

        # Cumulative RMS from the low-frequency end upward.
        inside = frequency >= 0.2
        cumulative = np.sqrt(np.concatenate((
            [0.0], np.cumsum(np.diff(frequency[inside])
                             * 0.5 * (data['psd'][inside, index][1:]
                                      + data['psd'][inside, index][:-1])))))
        ax2.semilogx(frequency[inside], cumulative * 1000.0,
                     color=AXIS_COLOUR[index], lw=1.3,
                     label=f'axis {index + 1}')
        table_rows.append(band_rms(frequency, data['psd'][:, index]))

    edge = float(np.nanmin(onsets))
    for panel in (ax, ax2):
        panel.axvspan(edge, 9000, color=WARM, alpha=0.07, lw=0)
        panel.grid(True, which='both', color=GRID, lw=0.5)
    ax.annotate(f'sensor noise floor\n(above {edge:.0f} Hz this is the\n'
                'instrument, not the ground)', (edge * 1.15, levels[0] * 40),
                fontsize=7.5, color=WARM, va='bottom')

    ax.set_ylabel('PSD [mV$^2$/Hz]')
    ax.set_title(f'(a)  Noise spectrum — {record}', loc='left', fontsize=10,
                 fontweight='bold')
    ax.legend(frameon=False, ncol=3, fontsize=9, loc='upper right')

    ax2.set_ylabel('cumulative RMS from 0.2 Hz [µV]')
    ax2.set_xlabel('frequency [Hz]')
    ax2.set_title('(b)  Cumulative RMS — where the motion actually comes from',
                  loc='left', fontsize=10, fontweight='bold')
    ax2.set_xlim(0.2, 9000)

    note = '  |  '.join(
        f'{lo:g}–{hi:g} Hz: ' + '/'.join(f'{r[i] * 1000:.0f}' for r in table_rows)
        for i, (lo, hi) in enumerate(BANDS))
    fig.suptitle(f'{record}  —  band RMS in µV, axes 1/2/3\n{note}',
                 fontsize=8.5, x=0.01, ha='left', color=MUTED)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out, table_rows, edge


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--fraction', type=float, default=12.0)
    args = parser.parse_args()

    style()
    args.out.mkdir(parents=True, exist_ok=True)
    for number, (key, cache, title) in enumerate(RECORDINGS, start=10):
        data = load(CACHE / cache)
        path, rows, edge = figure(
            title, data, args.fraction,
            args.out / f'seis_{number}_noise_{key}.png')
        print(f'wrote {path.name}   '
              f'({data["segments"]} averages, noise-limited above '
              f'{edge:.0f} Hz)')
        print('   %-12s %9s %9s %9s' % ('band [Hz]', 'axis 1', 'axis 2',
                                        'axis 3'))
        for i, (lo, hi) in enumerate(BANDS):
            print('   %-12s %9.1f %9.1f %9.1f'
                  % (f'{lo:g}-{hi:g}', *[r[i] * 1000 for r in rows]))
        print('   (values are RMS in microvolts)')
        print()


if __name__ == '__main__':
    main()
