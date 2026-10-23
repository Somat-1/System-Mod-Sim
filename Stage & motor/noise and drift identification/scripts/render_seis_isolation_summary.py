#!/usr/bin/env python3
'''Three plain figures answering "is the table isolating the stage?".

The diagnostic figure produced by compute_seis_transmissibility.py carries
every caveat on its face and is hard to read quickly. These three are the
same data arranged so that one finding is visible per figure:

  seis_06_isolation       the answer: transmissibility against frequency
  seis_07_bands           the same thing as numbers, band by band
  seis_08_psd_axes        the two recordings themselves, one panel per axis

Channels 1-3 are the three axes of the accelerometer. Table is the SEIS
recording, floor is the TvalSeisfloor recording.
'''

import argparse
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from compute_seis_transmissibility import (load, octave_smooth, BAND_LO,
                                           noise_limited_mask)

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_TABLE = ROOT / 'analysis' / 'cache' / 'TvalSeis_welch_10s.npz'
DEFAULT_FLOOR = ROOT / 'analysis' / 'cache' / 'TvalSeisfloor_welch_10s.npz'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'
# The two per-recording caches are produced by compute_seis_welch.py. When
# they are absent the combined 2-14 min Welch holds the same thing - both
# recordings, same grid, 10 s segments - and its settings file carries the
# metadata, so nothing has to be assumed.
PWELCH_CACHE = ROOT / 'analysis' / 'cache' / 'seis_pwelch_02_14min.npz'
PWELCH_SETTINGS = ROOT / 'analysis' / 'seis_pwelch_02_14min_settings.json'

AXIS_COLOUR = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#161616'
MUTED = '#6b6963'
GRID = '#e2e1dc'
WARM = '#bd3a32'
COOL = '#2f7f5b'

# The mains line is in the floor recording and not the table recording, so
# it is a difference between the two occasions, not isolation.
LINE_LO, LINE_HI = 47.5, 51.5
NOISE_HZ = 1850.0
BANDS = ((0.2, 1), (1, 5), (5, 10), (10, 20), (20, 50), (50, 100), (100, 300))


def style():
    mpl.rcParams.update({
        'font.size': 9, 'axes.edgecolor': MUTED, 'axes.labelcolor': INK,
        'text.color': INK, 'xtick.color': MUTED, 'ytick.color': MUTED,
        'axes.linewidth': 0.8, 'axes.spines.top': False,
        'axes.spines.right': False, 'figure.facecolor': 'white'})


def load_pair(table_path, floor_path):
    """The two per-recording caches if they exist, else the combined one."""
    if table_path.is_file() and floor_path.is_file():
        return load(table_path), load(floor_path)
    if not (PWELCH_CACHE.is_file() and PWELCH_SETTINGS.is_file()):
        raise SystemExit(
            f'no PSD cache: expected {table_path.name} and {floor_path.name}, '
            f'or {PWELCH_CACHE.name} with its settings file')
    cfg = json.loads(PWELCH_SETTINGS.read_text(encoding='utf-8'))
    with np.load(PWELCH_CACHE) as d:
        base = dict(frequency=d['frequency_hz'],
                    fs=float(cfg['sampling_frequency_hz']),
                    nperseg=int(cfg['segment_length_samples']),
                    segments=int(cfg['segment_count_per_recording']),
                    duration=float(cfg['analysis_duration_s']))
        return (dict(base, psd=d['table_psd_mV2_per_Hz']),
                dict(base, psd=d['floor_psd_mV2_per_Hz']))


def prepare(table, floor, fraction):
    frequency = table['frequency']
    axes = []
    for channel in range(3):
        t = octave_smooth(frequency, table['psd'][:, channel], fraction)
        f = octave_smooth(frequency, floor['psd'][:, channel], fraction)
        axes.append(dict(table=t, floor=f,
                         ratio_db=10.0 * np.log10(t / f)))
    return frequency, axes


def band_values(frequency, axes):
    keep = ~((frequency > LINE_LO) & (frequency < LINE_HI))
    rows = []
    for lo, hi in BANDS:
        mask = (frequency >= lo) & (frequency <= hi) & keep
        rows.append([10.0 * np.log10(
            np.trapezoid(a['table'][mask], frequency[mask])
            / np.trapezoid(a['floor'][mask], frequency[mask])) for a in axes])
    return np.array(rows)


def figure_isolation(frequency, axes, out):
    fig, ax = plt.subplots(figsize=(8.0, 4.4), constrained_layout=True)
    show = (frequency >= 0.2) & (frequency <= 500)

    ax.axhspan(0, 30, color=WARM, alpha=0.055, lw=0)
    ax.axhspan(-45, 0, color=COOL, alpha=0.055, lw=0)
    ax.axvspan(LINE_LO, LINE_HI, color=MUTED, alpha=0.16, lw=0)

    for index, a in enumerate(axes):
        ax.semilogx(frequency[show], a['ratio_db'][show],
                    color=AXIS_COLOUR[index], lw=1.4, label=f'axis {index + 1}')
    ax.axhline(0, color=INK, lw=1.0)

    ax.annotate('table AMPLIFIES floor motion', (0.23, 26), fontsize=8.5,
                color=WARM, va='top', fontweight='bold')
    ax.annotate('table ISOLATES the stage', (0.23, -41), fontsize=8.5,
                color=COOL, va='bottom', fontweight='bold')
    ax.annotate('49 Hz mains line\n(not isolation — see notes)', (56, -2.5),
                fontsize=7, color=MUTED, ha='left', va='top')
    ax.annotate('suspension modes\n2–16 Hz', (5.0, 29), fontsize=8,
                color=INK, ha='center', va='top')
    ax.annotate('', (1.9, 24.0), (17.0, 24.0),
                arrowprops=dict(arrowstyle='<->', lw=0.9, color=INK))

    ax.set_xlim(0.2, 500)
    ax.set_ylim(-45, 30)
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('transmissibility  table / floor  [dB]')
    ax.set_title('Vibration isolation of the optical table, measured on three '
                 'axes', loc='left', fontsize=10.5, fontweight='bold')
    ax.legend(frameon=False, ncol=3, loc='upper right', fontsize=9)
    ax.grid(True, which='both', color=GRID, lw=0.5)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def figure_bands(values, out):
    fig, ax = plt.subplots(figsize=(8.0, 3.8), constrained_layout=True)
    x = np.arange(len(BANDS))
    width = 0.26
    for index in range(3):
        ax.bar(x + (index - 1) * width, values[:, index], width,
               color=AXIS_COLOUR[index], label=f'axis {index + 1}')
    ax.axhline(0, color=INK, lw=1.0)
    ax.set_xticks(x)
    ax.set_xticklabels([f'{lo:g}–{hi:g}' for lo, hi in BANDS])
    ax.set_xlabel('frequency band [Hz]')
    ax.set_ylabel('transmissibility [dB]')
    ax.set_title('Band by band: positive is amplification, negative is '
                 'isolation', loc='left', fontsize=10.5, fontweight='bold')
    ax.legend(frameon=False, ncol=3, fontsize=9)
    ax.grid(True, axis='y', color=GRID, lw=0.5)
    ax.set_axisbelow(True)
    ax.set_ylim(values.min() - 4.0, values.max() + 4.0)
    for index in range(3):
        for i, v in enumerate(values[:, index]):
            ax.annotate(f'{v:+.0f}', (i + (index - 1) * width, v),
                        textcoords='offset points',
                        xytext=(0, 3 if v >= 0 else -11),
                        ha='center', fontsize=6.5, color=MUTED)
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def figure_psd(frequency, axes, out):
    """One row per axis, floor and table overlaid, on a shared frequency axis.

    Stacked rather than side by side so the three rows line up in frequency:
    a feature at a given frequency can be followed down the column, and the
    gap between the two curves - which is the whole point - is read at the
    same scale in every row.
    """
    fig, panels = plt.subplots(3, 1, figsize=(7.2, 7.8), sharex=True,
                               sharey=True, constrained_layout=True)
    show = (frequency >= 0.2) & (frequency <= 1000)
    for index, (ax, a) in enumerate(zip(panels, axes)):
        ax.loglog(frequency[show], a['floor'][show], color=MUTED, lw=1.2,
                  label='floor', zorder=2)
        ax.loglog(frequency[show], a['table'][show],
                  color=AXIS_COLOUR[index], lw=1.3, label='table', zorder=3)
        ax.set_title(f'axis {index + 1}', loc='left', fontsize=9.5,
                     fontweight='bold')
        ax.grid(True, which='major', linestyle='-', color='#cfcfca',
                linewidth=0.62, alpha=0.78)
        ax.grid(True, which='minor', linestyle='-', color='#e1e0dc',
                linewidth=0.42, alpha=0.72)
        ax.set_axisbelow(True)
        ax.set_ylabel('PSD [mV$^2$/Hz]')
        if index == 0:
            ax.legend(frameon=False, fontsize=8.5, loc='lower left')
    panels[-1].set_xlabel('frequency [Hz]')
    fig.suptitle('What each sensor actually recorded\n'
                 'where the table curve sits above the floor curve, '
                 'the table is amplifying',
                 fontsize=10, x=0.012, ha='left')
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def noise_floor_onset(frequency, axes, margin_db=10.0):
    """Lowest frequency at which any axis comes within margin_db of the floor.

    The instrument's own noise shows as the level the PSD settles to once the
    signal has fallen away; 4-8 kHz is where all six channels sit on it, and
    the two recordings agree there to 1.2 %, which is what identifies it as
    the sensor rather than the ground. Below a decade of margin the ratio's
    uncertainty grows quickly, so that is where the band starts.
    """
    flat = (frequency >= 4000.0) & (frequency <= 8000.0)
    limit = 10.0 ** (margin_db / 10.0)
    onsets = []
    for a in axes:
        level = float(np.median(a['table'][flat]))
        below = np.where((a['table'] <= limit * level) & (frequency > 200.0))[0]
        if below.size:
            onsets.append(float(frequency[below[0]]))
    return min(onsets) if onsets else float('nan')


def figure_transmissibility(frequency, axes, out):
    """Amplitude transmissibility T(f) = sqrt(PSD_table / PSD_floor).

    The curves and nothing else: no shaded regimes, no callouts. Unity stays
    as a plain reference line because it is what the ratio is read against -
    above it the table amplifies the floor, below it the table isolates.

    Note for the caption: the dip near 50 Hz is not isolation. The mains line
    is present in the floor recording and absent from the table one, so the
    ratio collapses there for a reason that has nothing to do with the
    isolator. The band is no longer marked on the figure.

    The sensor noise floor is checked rather than assumed. Over 0.2-500 Hz the
    table recording sits 321x to 1818x above the instrument's own floor at its
    quietest, so every point here is a real measurement; the faded styling
    only appears if a future band does run into that limit.
    """
    fig, ax = plt.subplots(figsize=(8.0, 4.9), constrained_layout=True)
    fmax = 2000.0
    show = (frequency >= 0.2) & (frequency <= fmax)

    onset = noise_floor_onset(frequency, axes)
    if np.isfinite(onset) and onset < fmax:
        ax.axvspan(onset, fmax, color=WARM, alpha=0.10, lw=0, zorder=0)
        ax.annotate('within 10 dB of the sensor noise floor',
                    ((onset * fmax) ** 0.5, 2.6e-2), fontsize=7.5,
                    color=WARM, ha='center', va='center', rotation=90)

    for index, a in enumerate(axes):
        t = np.sqrt(a['table'] / a['floor'])
        limited, _ = noise_limited_mask(a['table'], frequency)
        ax.loglog(frequency, np.where(show & limited, t, np.nan),
                  color=AXIS_COLOUR[index], lw=0.9, alpha=0.35, zorder=3)
        ax.loglog(frequency, np.where(show & ~limited, t, np.nan),
                  color=AXIS_COLOUR[index], lw=1.4, zorder=4,
                  label=f'axis {index + 1}')

    ax.axhline(1.0, color=INK, lw=1.1, zorder=5)
    ax.grid(True, which='major', linestyle='-', color='#cfcfca',
            linewidth=0.62, alpha=0.78)
    ax.grid(True, which='minor', linestyle='-', color='#e1e0dc',
            linewidth=0.42, alpha=0.72)
    ax.set_axisbelow(True)

    ax.set_xlim(0.2, fmax)
    ax.set_ylim(1.2e-3, 30)
    ax.set_xlabel('frequency [Hz]')
    ax.set_ylabel('transmissibility  $T=\\sqrt{S_\\mathrm{table}/'
                  'S_\\mathrm{floor}}$')
    ax.legend(frameon=False, ncol=3, loc='lower left', fontsize=9)
    ax.set_title('Vibration transmissibility of the optical table, '
                 'three axes', loc='left', fontsize=10.5, fontweight='bold')
    fig.savefig(out, dpi=200)
    plt.close(fig)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, default=DEFAULT_TABLE)
    parser.add_argument('--floor', type=Path, default=DEFAULT_FLOOR)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--fraction', type=float, default=12.0)
    args = parser.parse_args()

    style()
    table, floor = load_pair(args.table, args.floor)
    frequency, axes = prepare(table, floor, args.fraction)
    values = band_values(frequency, axes)

    args.out.mkdir(parents=True, exist_ok=True)
    for path in (figure_isolation(frequency, axes,
                                  args.out / 'seis_06_isolation.png'),
                 figure_bands(values, args.out / 'seis_07_bands.png'),
                 figure_psd(frequency, axes,
                            args.out / 'seis_08_psd_axes.png'),
                 figure_transmissibility(
                     frequency, axes,
                     args.out / 'seis_09_transmissibility.png')):
        print(f'wrote {path.name}')

    print()
    print('%-12s %8s %8s %8s' % ('band [Hz]', 'axis 1', 'axis 2', 'axis 3'))
    for (lo, hi), row in zip(BANDS, values):
        print('%-12s %+7.1f %+7.1f %+7.1f'
              % (f'{lo:g}-{hi:g}', row[0], row[1], row[2]))


if __name__ == '__main__':
    main()
