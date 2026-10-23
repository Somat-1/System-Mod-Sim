#!/usr/bin/env python3
'''Table-over-floor transmissibility from two separate SEIS Welch PSDs.

T(f) = sqrt(PSD_table / PSD_floor), reported in dB as 20*log10(T).

IMPORTANT, and stated again on the figure: the table and the floor were
recorded one after the other with the same instrument, not simultaneously.
A true transmissibility is a dual-channel measurement, where the ratio is
formed from the cross-spectrum of two sensors seeing the same excitation at
the same moment, and is qualified by coherence. What this script computes is
the ratio of two ambient spectra taken at different times. It is only a
transmissibility to the extent that the floor excitation was statistically
the same during both recordings, which these data cannot confirm. Narrow
features that come and go -- a pump cycling, a lift, a 49 Hz line -- break
that assumption badly; broadband behaviour survives it reasonably well.

The instrument is the same in both recordings, so wherever both are limited
by the sensor's own noise the ratio must approach 0 dB regardless of what
the table is doing. That band is detected and shaded rather than read as
isolation performance.
'''

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from scipy.signal import find_peaks

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DEFAULT_TABLE = ROOT / 'analysis' / 'cache' / 'TvalSeis_welch_10s.npz'
DEFAULT_FLOOR = ROOT / 'analysis' / 'cache' / 'TvalSeisfloor_welch_10s.npz'
DEFAULT_OUT = ROOT / 'analysis' / 'plots'
DEFAULT_REPORT = ROOT / 'analysis' / 'SEIS_TRANSMISSIBILITY.md'

COLOURS = ('#2a78d6', '#e96b36', '#1b9e77')
INK = '#161616'
MUTED = '#696762'
GRID = '#dfded9'
ALERT = '#bd3a32'

# The band a pneumatic isolator is judged on. Below the first number the
# sensor rolls off and the estimate is not trustworthy; above the last the
# instrument noise dominates on the table side.
BAND_LO, BAND_HI = 0.2, 400.0
RESONANCE_BAND = (0.5, 8.0)


def octave_smooth(frequency, values, fraction=12.0, min_bins=3):
    '''Smooth on a log-frequency axis, which is how the result is read.

    A fixed number of bins over-smooths the decade that matters (1-3 Hz) and
    under-smooths the top one, so the kernel is a constant fraction of an
    octave with a floor of a few bins so the lowest frequencies still average.
    '''
    out = np.empty_like(values)
    half = 2.0 ** (0.5 / fraction)
    df = frequency[1] - frequency[0]
    for i, f in enumerate(frequency):
        if f <= 0:
            out[i] = values[i]
            continue
        lo, hi = f / half, f * half
        a = np.searchsorted(frequency, lo, 'left')
        b = np.searchsorted(frequency, hi, 'right')
        if b - a < min_bins:
            pad = (min_bins - (b - a) + 1) // 2
            a, b = max(0, i - pad), min(len(frequency), i + pad + 1)
        out[i] = np.mean(values[a:b])
    return out


def load(path):
    with np.load(path) as data:
        return dict(frequency=data['frequency_hz'],
                    psd=data['psd_mV2_per_Hz'],
                    fs=float(data['fs_hz']),
                    nperseg=int(data['nperseg']),
                    segments=int(data['segment_count']),
                    duration=float(data['duration_s']))


def noise_limited_mask(table_psd, frequency, top_decade=(2000.0, 9000.0),
                       margin_db=3.0):
    '''Where the table recording has stopped measuring the table.

    The instrument's own noise shows up as a floor the PSD cannot go below.
    Taking the level it settles to in the top decade and marking everything
    within a few dB of it flags the band where the ratio is an instrument
    artefact rather than a property of the isolator.
    '''
    window = (frequency >= top_decade[0]) & (frequency <= top_decade[1])
    if not np.any(window):
        return np.zeros_like(frequency, dtype=bool), np.nan
    floor_level = np.median(table_psd[window])
    limit = floor_level * 10.0 ** (margin_db / 10.0)
    return table_psd <= limit, floor_level


def analyse(table, floor, fraction):
    frequency = table['frequency']
    if len(floor['frequency']) != len(frequency):
        raise SystemExit('the two PSDs are on different grids; recompute both '
                         'with the same --nperseg')

    result = dict(frequency=frequency, channels=[])
    for channel in range(table['psd'].shape[1]):
        t_raw = table['psd'][:, channel]
        f_raw = floor['psd'][:, channel]
        t_s = octave_smooth(frequency, t_raw, fraction)
        f_s = octave_smooth(frequency, f_raw, fraction)
        with np.errstate(divide='ignore', invalid='ignore'):
            ratio_db = 10.0 * np.log10(t_s / f_s)
        # The table is the quieter recording, so it reaches the instrument
        # floor first and is normally what limits the usable band. Where it
        # does, the ratio is biased HIGH: the table PSD is reading the
        # sensor, not the table, so the true isolation is better than shown.
        t_noise, t_level = noise_limited_mask(t_s, frequency)
        f_noise, f_level = noise_limited_mask(f_s, frequency)
        noise_mask = t_noise | f_noise
        floor_level = t_level
        onset = dict(
            table=float(frequency[np.argmax(t_noise)]) if t_noise.any()
            else float('nan'),
            floor=float(frequency[np.argmax(f_noise)]) if f_noise.any()
            else float('nan'))

        band = (frequency >= RESONANCE_BAND[0]) & (frequency <= RESONANCE_BAND[1])
        peak_index = np.argmax(np.where(band, ratio_db, -np.inf))
        usable = (frequency >= BAND_LO) & (frequency <= BAND_HI) & ~noise_mask

        # The first crossing into attenuation above the resonance, which is
        # where the isolator starts doing its job.
        crossing = np.nan
        above = frequency > frequency[peak_index]
        candidates = np.where(above & (ratio_db < 0) & usable)[0]
        if candidates.size:
            crossing = float(frequency[candidates[0]])

        result['channels'].append(dict(
            index=channel, table=t_s, floor=f_s, ratio_db=ratio_db,
            table_raw=t_raw, floor_raw=f_raw,
            noise_mask=noise_mask, noise_floor_level=floor_level,
            floor_noise_level=f_level, noise_onset=onset,
            peak_hz=float(frequency[peak_index]),
            peak_db=float(ratio_db[peak_index]),
            crossing_hz=crossing,
            at=dict((f, float(np.interp(f, frequency, ratio_db)))
                    for f in (1.0, 2.0, 5.0, 10.0, 20.0, 50.0, 100.0, 200.0))))
    return result


def render(result, table, floor, fraction):
    mpl.rcParams.update({'font.size': 8.5, 'axes.edgecolor': MUTED,
                         'axes.labelcolor': INK, 'text.color': INK,
                         'xtick.color': MUTED, 'ytick.color': MUTED,
                         'axes.linewidth': 0.7})
    fig, axes = plt.subplots(3, 1, figsize=(7.4, 8.6), sharex=True,
                             constrained_layout=True)
    frequency = result['frequency']
    ax_psd, ax_t, ax_n = axes

    for channel in result['channels']:
        colour = COLOURS[channel['index']]
        label = f"ch {channel['index'] + 1}"
        ax_psd.loglog(frequency, channel['table'], color=colour, lw=1.0,
                      label=f'{label} table')
        ax_psd.loglog(frequency, channel['floor'], color=colour, lw=0.9,
                      ls='--', alpha=0.75, label=f'{label} floor')
        ax_t.semilogx(frequency, channel['ratio_db'], color=colour, lw=1.1,
                      label=label)
        ax_n.loglog(frequency, channel['table'], color=colour, lw=1.0)

    ax_psd.set_ylabel('PSD [mV$^2$/Hz]')
    ax_psd.set_title('(a)  Welch PSD, table (solid) and floor (dashed)',
                     loc='left', fontsize=9)
    ax_psd.legend(ncol=3, fontsize=6.5, frameon=False)
    ax_psd.grid(True, which='both', color=GRID, lw=0.4)

    ax_t.axhline(0, color=INK, lw=0.8)
    ax_t.set_ylabel('transmissibility [dB]')
    ax_t.set_title('(b)  T = 20·log₁₀ √(PSD$_{table}$/PSD$_{floor}$) — '
                   'above 0 dB the table amplifies', loc='left', fontsize=9)
    ax_t.grid(True, which='both', color=GRID, lw=0.4)
    for channel in result['channels']:
        ax_t.plot(channel['peak_hz'], channel['peak_db'], 'o', ms=4,
                  color=COLOURS[channel['index']], mec='white', mew=0.6,
                  zorder=5)

    # The instrument floor, shaded on every panel it compromises.
    mask = np.any([c['noise_mask'] for c in result['channels']], axis=0)
    if np.any(mask):
        edge = frequency[np.argmax(mask)] if mask.any() else np.nan
        for ax in (ax_t, ax_n):
            ax.axvspan(edge, frequency[-1], color=ALERT, alpha=0.08, lw=0)
        ax_t.annotate('table at instrument noise floor:\nratio here is not '
                      'isolation', (edge, ax_t.get_ylim()[0]),
                      xytext=(6, 8), textcoords='offset points',
                      fontsize=7, color=ALERT, va='bottom')

    for channel in result['channels']:
        if np.isfinite(channel['noise_floor_level']):
            ax_n.axhline(channel['noise_floor_level'],
                         color=COLOURS[channel['index']], lw=0.6, ls=':')
    ax_n.set_ylabel('PSD [mV$^2$/Hz]')
    ax_n.set_xlabel('frequency [Hz]')
    ax_n.set_title('(c)  Table PSD against its own high-frequency asymptote '
                   '(dotted) — the sensor noise floor', loc='left', fontsize=9)
    ax_n.grid(True, which='both', color=GRID, lw=0.4)
    ax_n.set_xlim(frequency[1], frequency[-1])

    fig.suptitle('SEIS table-over-floor spectral ratio — '
                 f"{table['segments']} table and {floor['segments']} floor "
                 f"Welch averages, {table['nperseg'] / table['fs']:.0f} s "
                 f"Hann segments, 50 % overlap, 1/{fraction:.0f}-octave "
                 'smoothing\nSEQUENTIAL recordings: this is a spectral ratio, '
                 'not a coherence-qualified transmissibility',
                 fontsize=9, x=0.01, ha='left')
    return fig


def write_report(path, result, table, floor, fraction):
    lines = ['# SEIS table-over-floor transmissibility', '',
             '## Method', '',
             f"- Welch, periodic Hann, 50 % overlap, "
             f"{table['nperseg']} samples = "
             f"{table['nperseg'] / table['fs']:.1f} s per segment",
             f"- frequency resolution {table['fs'] / table['nperseg']:.4f} Hz",
             f"- table: {table['segments']} averages over "
             f"{table['duration']:.1f} s; floor: {floor['segments']} averages "
             f"over {floor['duration']:.1f} s",
             f"- displayed with 1/{fraction:.0f}-octave smoothing",
             '- T(f) = sqrt(PSD_table / PSD_floor), reported as 20*log10(T)',
             '',
             '> The two recordings were made sequentially with one instrument,',
             '> not simultaneously. This is therefore a ratio of ambient',
             '> spectra, valid as a transmissibility only insofar as the floor',
             '> excitation was statistically unchanged between them. No',
             '> coherence can be computed from single-sensor data.', '',
             '## Per channel', '',
             '| Channel | Resonance [Hz] | Peak T [dB] | 0 dB crossing [Hz] |'
             ' T(10 Hz) | T(20 Hz) | T(50 Hz) | T(100 Hz) |',
             '|---|---:|---:|---:|---:|---:|---:|---:|']
    for channel in result['channels']:
        at = channel['at']
        crossing = ('%.2f' % channel['crossing_hz']
                    if np.isfinite(channel['crossing_hz']) else 'n/a')
        lines.append(
            f"| {channel['index'] + 1} | {channel['peak_hz']:.2f} | "
            f"{channel['peak_db']:+.1f} | {crossing} | "
            f"{at[10.0]:+.1f} | {at[20.0]:+.1f} | {at[50.0]:+.1f} | "
            f"{at[100.0]:+.1f} |")

    # The mains line sits in one recording and not the other, so every
    # band statistic steps around it rather than averaging it in.
    frequency = result['frequency']
    keep = ~((frequency > 47.5) & (frequency < 51.5))
    lines += ['', '## Rigid-body modes, 0.5-25 Hz', '',
              'Peaks in the smoothed ratio, by prominence. Two axes seeing',
              'the same mode is expected; a mode on one axis only is a',
              'rotation the other axis does not project onto.', '',
              '| Channel | Peaks [Hz (dB)] |', '|---|---|']
    for channel in result['channels']:
        band = (frequency >= 0.5) & (frequency <= 25.0)
        db = channel['ratio_db']
        idx, _ = find_peaks(db[band], prominence=2.5)
        sub_f, sub_db = frequency[band], db[band]
        peaks = '  '.join(f'{sub_f[i]:.2f} ({sub_db[i]:+.1f})' for i in idx)
        lines.append(f"| {channel['index'] + 1} | {peaks} |")

    lines += ['', '## Band-integrated ratio', '',
              'Ratio of band-integrated PSD, 47.5-51.5 Hz excluded because',
              'the 49.02 Hz line is present in the floor recording and',
              'essentially absent from the table recording, which is a',
              'difference between the two occasions rather than isolation.',
              '', '| Band [Hz] | Ch 1 [dB] | Ch 2 [dB] | Ch 3 [dB] |',
              '|---|---:|---:|---:|']
    for lo, hi in ((0.2, 1), (1, 5), (5, 10), (10, 20), (20, 50),
                   (50, 100), (100, 300), (300, 1000)):
        mask = (frequency >= lo) & (frequency <= hi) & keep
        row = []
        for channel in result['channels']:
            row.append(10.0 * np.log10(
                np.trapezoid(channel['table'][mask], frequency[mask])
                / np.trapezoid(channel['floor'][mask], frequency[mask])))
        lines.append(f"| {lo:g}-{hi:g} | " +
                     ' | '.join(f'{v:+.1f}' for v in row) + ' |')

    lines += ['', '## Horizontal symmetry', '',
              'Channels 1 and 2 share every mode frequency, so they are two',
              'axes of the same rigid body. A defective leg would be expected',
              'to break that symmetry.', '',
              '| Band [Hz] | Ch 1 - Ch 2 [dB] |', '|---|---:|']
    for lo, hi in ((1, 5), (5, 10), (10, 20), (20, 100)):
        mask = (frequency >= lo) & (frequency <= hi) & keep
        delta = float(np.mean(result['channels'][0]['ratio_db'][mask]
                              - result['channels'][1]['ratio_db'][mask]))
        lines.append(f'| {lo:g}-{hi:g} | {delta:+.1f} |')

    lines += ['', '## Instrument noise floor', '',
              'The frequency above which each recording is within 3 dB of its',
              'own high-frequency asymptote, i.e. is measuring the sensor',
              'rather than the ground. The table is the quieter recording and',
              'reaches its floor first; above that frequency the ratio is',
              'biased high and the true isolation is better than plotted.', '',
              '| Channel | Table noise-limited above [Hz] | Floor '
              'noise-limited above [Hz] | Table floor level [mV^2/Hz] |',
              '|---|---:|---:|---:|']
    for channel in result['channels']:
        onset = channel['noise_onset']
        lines.append(
            f"| {channel['index'] + 1} | {onset['table']:.0f} | "
            f"{onset['floor']:.0f} | {channel['noise_floor_level']:.3g} |")
    path.write_text('\n'.join(lines) + '\n', encoding='utf-8')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--table', type=Path, default=DEFAULT_TABLE)
    parser.add_argument('--floor', type=Path, default=DEFAULT_FLOOR)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--report', type=Path, default=DEFAULT_REPORT)
    parser.add_argument('--fraction', type=float, default=12.0,
                        help='smoothing kernel, in 1/N octave')
    parser.add_argument('--output-stem', default='seis_09_transmissibility_diagnostic')
    parser.add_argument('--dpi', type=int, default=180)
    args = parser.parse_args()

    table, floor = load(args.table), load(args.floor)
    result = analyse(table, floor, args.fraction)

    args.out.mkdir(parents=True, exist_ok=True)
    fig = render(result, table, floor, args.fraction)
    png = args.out / f'{args.output_stem}.png'
    fig.savefig(png, dpi=args.dpi)
    plt.close(fig)
    write_report(args.report, result, table, floor, args.fraction)

    print(f'wrote {png}')
    print(f'wrote {args.report}')
    for channel in result['channels']:
        at = channel['at']
        print(f"  ch {channel['index'] + 1}: resonance "
              f"{channel['peak_hz']:.2f} Hz at {channel['peak_db']:+.1f} dB | "
              f"0 dB crossing {channel['crossing_hz']:.2f} Hz | "
              f"10 Hz {at[10.0]:+.1f} | 50 Hz {at[50.0]:+.1f} | "
              f"100 Hz {at[100.0]:+.1f} dB | noise-limited above "
              f"{channel['noise_onset']['table']:.0f} Hz")


if __name__ == '__main__':
    main()
