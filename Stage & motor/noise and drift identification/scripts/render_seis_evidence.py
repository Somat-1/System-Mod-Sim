#!/usr/bin/env python3
'''The three-figure case: raw signals, spectra, transmissibility.

  seis_12_raw_time.png   floor and table waveforms, same axes and same scale
  seis_13_asd.png        ASD of both overlaid, with the sensor noise floor
                         as a third line and the floor's block-to-block
                         spread as a shaded band
  seis_14_transmissibility.png
                         the ratio in dB against log frequency, 0 dB marked,
                         and the band where the measurement is not valid
                         greyed out

Units. The recordings are voltages. Pass --sensitivity in V per m/s^2 to
relabel the ASD figure in (m/s^2)/sqrt(Hz); without it the axis stays in
mV/sqrt(Hz) and says so. A constant sensitivity cancels in the ratio, so the
transmissibility figure is unaffected either way.
'''

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from compute_seis_transmissibility import load, octave_smooth

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
CACHE = ROOT / 'analysis' / 'cache'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'

AXIS_COLOUR = ('#2a78d6', '#e96b36', '#1b9e77')
TABLE_C, FLOOR_C = '#2a78d6', '#8a8781'
INK, MUTED, GRID, WARM = '#161616', '#6b6963', '#e2e1dc', '#bd3a32'
LINE_LO, LINE_HI = 47.5, 51.5


def style():
    mpl.rcParams.update({
        'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
        'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
        'axes.linewidth': 0.8, 'axes.spines.top': False,
        'axes.spines.right': False, 'figure.facecolor': 'white'})


def noise_floor(frequency, psd, top=(3000.0, 9000.0), margin_db=3.0):
    window = (frequency >= top[0]) & (frequency <= top[1])
    level = float(np.median(psd[window]))
    reached = psd <= level * 10.0 ** (margin_db / 10.0)
    onset = float(frequency[np.argmax(reached)]) if reached.any() else np.nan
    return level, onset


def validity(frequency, table_psd, floor_psd):
    '''Where the ratio can still be believed.

    Each recording has its own noise floor, and the ratio is only as good as
    the weaker of the two signal-to-noise ratios. A recording sitting 10 dB
    above its floor is biased high by 0.4 dB, which is tolerable; at 3 dB the
    bias is 1.8 dB, which is not. Above the second boundary the ratio simply
    reports the ratio of the two instrument noise floors.
    '''
    t_level, _ = noise_floor(frequency, table_psd)
    f_level, _ = noise_floor(frequency, floor_psd)
    snr = np.minimum(table_psd / t_level, floor_psd / f_level)
    snr_db = 10.0 * np.log10(snr)
    high = frequency > 100.0

    def first_below(limit):
        bad = high & (snr_db < limit)
        return float(frequency[np.argmax(bad)]) if bad.any() else np.nan

    return dict(marginal=first_below(10.0), invalid=first_below(3.0),
                table_level=t_level, floor_level=f_level,
                ratio_db=10.0 * np.log10(t_level / f_level))


def scale_and_label(sensitivity):
    '''Volts to acceleration, if the sensitivity is known.

    PSD is in mV^2/Hz. Dividing the amplitude by 1000*S converts mV to
    m/s^2, so the ASD scale factor is 1/(1000*S).
    '''
    if sensitivity is None:
        return 1.0, 'ASD [mV/$\\sqrt{Hz}$]', ' (uncalibrated — volts)'
    return (1.0 / (1000.0 * sensitivity),
            'ASD [(m/s$^2$)/$\\sqrt{Hz}$]',
            f' (sensitivity {sensitivity:g} V per m/s²)')


def figure_raw(table, floor, out):
    fig, panels = plt.subplots(3, 1, figsize=(8.2, 5.4), sharex=True,
                               sharey=True, constrained_layout=True)
    t_snip, f_snip = table['snippet'], floor['snippet']
    t0, f0 = t_snip[0, 0], f_snip[0, 0]
    for index, ax in enumerate(panels):
        ax.plot(f_snip[:, 0] - f0, f_snip[:, index + 1] - f_snip[:, index + 1].mean(),
                color=FLOOR_C, lw=0.6, label='floor' if index == 0 else None)
        ax.plot(t_snip[:, 0] - t0, t_snip[:, index + 1] - t_snip[:, index + 1].mean(),
                color=AXIS_COLOUR[index], lw=0.6,
                label='table' if index == 0 else None)
        ax.set_ylabel(f'axis {index + 1}\n[mV]')
        ax.grid(True, color=GRID, lw=0.5)
        if index == 0:
            ax.legend(frameon=False, ncol=2, fontsize=9, loc='upper right')
    panels[-1].set_xlabel('time [s]')
    panels[0].set_title('Raw signals, same scale — the floor moves visibly '
                        'more than the table', loc='left', fontsize=10.5,
                        fontweight='bold')
    fig.suptitle(f'5 s from each recording (table from {t0:.0f} s, floor '
                 f'from {f0:.0f} s), means removed, no filtering',
                 fontsize=8.5, x=0.01, ha='left', color=MUTED)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def figure_asd(frequency, table, floor, floor_blocks, sensitivity, out):
    scale, ylabel, note = scale_and_label(sensitivity)
    fig, panels = plt.subplots(1, 3, figsize=(10.4, 4.0), sharex=True,
                               sharey=True, constrained_layout=True)
    show = (frequency >= 0.2) & (frequency <= 9000)

    for index, ax in enumerate(panels):
        t_psd = octave_smooth(frequency, table['psd'][:, index])
        f_psd = octave_smooth(frequency, floor['psd'][:, index])
        level, onset = noise_floor(frequency, t_psd)

        if floor_blocks is not None:
            band = np.array([octave_smooth(frequency, b[:, index])
                             for b in floor_blocks])
            ax.fill_between(frequency[show],
                            np.sqrt(band.min(axis=0)[show]) * scale,
                            np.sqrt(band.max(axis=0)[show]) * scale,
                            color=FLOOR_C, alpha=0.30, lw=0,
                            label='floor, block spread' if index == 0 else None)

        ax.loglog(frequency[show], np.sqrt(f_psd[show]) * scale,
                  color=FLOOR_C, lw=1.3,
                  label='floor' if index == 0 else None)
        ax.loglog(frequency[show], np.sqrt(t_psd[show]) * scale,
                  color=AXIS_COLOUR[index], lw=1.3,
                  label='table' if index == 0 else None)
        # The two recordings do NOT share a noise floor even though they
        # share an instrument, so both are drawn. The table sits higher,
        # which is why the ratio turns positive at high frequency.
        f_level, _ = noise_floor(frequency, f_psd)
        ax.axhline(np.sqrt(level) * scale, color=WARM, lw=1.0, ls='--',
                   label='table noise floor' if index == 0 else None)
        ax.axhline(np.sqrt(f_level) * scale, color=WARM, lw=1.0, ls=':',
                   label='floor noise floor' if index == 0 else None)
        ax.axvspan(onset, 9000, color=WARM, alpha=0.07, lw=0)

        ax.set_title(f'axis {index + 1}', loc='left', fontsize=9.5,
                     fontweight='bold')
        ax.set_xlabel('frequency [Hz]')
        ax.grid(True, which='both', color=GRID, lw=0.5)
        if index == 0:
            ax.set_ylabel(ylabel)
            ax.legend(frameon=False, fontsize=8, loc='lower left')

    fig.suptitle(
        'Floor against table, amplitude spectral density' + note +
        '\nGrey band: spread across 6 equal blocks of the floor recording — '
        'within-record stationarity.  Red shading: the table is at its own '
        'noise floor, comparison not valid.'
        '\nThe two noise floors differ (table 59–74, floor 37–40 nV/√Hz) '
        'although the instrument is the same — the table recording picked up '
        'more electrical noise.',
        fontsize=8, x=0.01, ha='left', color=MUTED)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def figure_transmissibility(frequency, table, floor, out):
    fig, ax = plt.subplots(figsize=(8.6, 4.6), constrained_layout=True)
    show = (frequency >= 0.2) & (frequency <= 9000)

    marginal, invalid, floor_ratio = [], [], []
    for index in range(3):
        t_psd = octave_smooth(frequency, table['psd'][:, index])
        f_psd = octave_smooth(frequency, floor['psd'][:, index])
        check = validity(frequency, t_psd, f_psd)
        marginal.append(check['marginal'])
        invalid.append(check['invalid'])
        floor_ratio.append(check['ratio_db'])
        ax.semilogx(frequency[show], 10.0 * np.log10(t_psd / f_psd)[show],
                    color=AXIS_COLOUR[index], lw=1.3, label=f'axis {index + 1}')

    soft, hard = float(np.nanmin(marginal)), float(np.nanmin(invalid))
    ax.axhline(0, color=INK, lw=1.1)
    ax.axvspan(soft, 9000, color=WARM, alpha=0.06, lw=0)
    ax.axvspan(hard, 9000, color=WARM, alpha=0.13, lw=0)
    ax.axvspan(LINE_LO, LINE_HI, color=MUTED, alpha=0.20, lw=0)
    ax.annotate(f'marginal\nfrom {soft:.0f} Hz', (soft * 0.97, -27),
                fontsize=7.5, color=WARM, va='top', ha='right')
    ax.annotate(f'NOT VALID above {hard:.0f} Hz\nboth recordings sit at their '
                f'own\nnoise floors, so the curve here is\nonly the ratio of '
                f'those floors ({np.mean(floor_ratio):+.1f} dB)',
                (hard * 1.1, -27), fontsize=7.5, color=WARM, va='top',
                fontweight='bold')
    ax.annotate('49 Hz mains line\n(present on the floor only)', (53, -36),
                fontsize=7.5, color=MUTED, ha='left', va='bottom')
    ax.annotate('0 dB — table moves with the floor', (0.22, 1.2), fontsize=8,
                color=INK, va='bottom')

    ax.set_xlim(0.2, 9000)
    ax.set_ylim(-45, 30)
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('transmissibility  20·log₁₀√(PSD$_{table}$/PSD$_{floor}$)  [dB]')
    ax.set_title('Transmissibility of the optical table', loc='left',
                 fontsize=11, fontweight='bold')
    ax.legend(frameon=False, ncol=3, fontsize=9, loc='upper right')
    ax.grid(True, which='both', color=GRID, lw=0.5)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sensitivity', type=float, default=None,
                        help='V per m/s^2; omit to keep the ASD in mV')
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    style()
    table = load(CACHE / 'TvalSeis_welch_10s.npz')
    floor = load(CACHE / 'TvalSeisfloor_welch_10s.npz')
    frequency = table['frequency']

    detail = {}
    for key, name in (('table', 'TvalSeis_detail.npz'),
                      ('floor', 'TvalSeisfloor_detail.npz')):
        path = CACHE / name
        if path.exists():
            with np.load(path) as data:
                detail[key] = dict(snippet=data['snippet'],
                                   blocks=data['block_psd'])

    args.out.mkdir(parents=True, exist_ok=True)
    written = []
    if 'table' in detail and 'floor' in detail:
        written.append(figure_raw(detail['table'], detail['floor'],
                                  args.out / 'seis_12_raw_time.png'))
    written.append(figure_asd(
        frequency, table, floor,
        detail.get('floor', {}).get('blocks'), args.sensitivity,
        args.out / 'seis_13_asd.png'))
    written.append(figure_transmissibility(
        frequency, table, floor, args.out / 'seis_14_transmissibility.png'))

    for path in written:
        print(f'wrote {path.name}')
    if args.sensitivity is None:
        print('\nASD is in mV/sqrt(Hz). Rerun with --sensitivity <V per m/s^2> '
              'to label it in (m/s^2)/sqrt(Hz).')


if __name__ == '__main__':
    main()
