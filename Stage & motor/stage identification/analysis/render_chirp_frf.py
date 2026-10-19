#!/usr/bin/env python3
'''Render absolute and level-normalised chirp response curves per recording.

There is no measured excitation channel in the Dewesoft files.  The absolute
curves are therefore an FRF proxy: narrow-band RMS acceleration at the known
drive frequency, not acceleration divided by measured force or displacement.
The normalised view divides every sweep and accelerometer axis by its own
median response over the measured 250--1000 Hz band and reports the result in
dB.  It compares response shape while making no claim about input level.
'''

import argparse
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

import vizstyle as vs


HERE = Path(__file__).resolve().parent
AXES = ('axis 1', 'axis 2', 'axis 3')


def curves_for(recording):
    curves = []
    for run_number, run in enumerate(recording.get('runs', []), 1):
        for leg in run.get('legs', []):
            frequency = np.asarray(leg['frequency'], dtype=float)
            for axis in AXES:
                response = np.asarray([
                    np.nan if value is None else value
                    for value in leg['response'][axis]
                ], dtype=float)
                curves.append(dict(run=run_number,
                                   direction=leg['direction'], axis=axis,
                                   frequency=frequency, response=response))
    return curves


def finite_positive(curve):
    keep = (np.isfinite(curve['frequency']) &
            np.isfinite(curve['response']) &
            (curve['response'] > 0))
    return curve['frequency'][keep], curve['response'][keep]


def page(recording, normalised=False):
    name = recording['recording']
    curves = curves_for(recording)
    runs = sorted({curve['run'] for curve in curves})
    fig, axes = plt.subplots(1, 3, figsize=(11.69, 8.27), sharex=True,
                             sharey=normalised)
    fig.subplots_adjust(left=0.075, right=0.985, top=0.83, bottom=0.17,
                        wspace=0.17)

    kind = 'Level-normalised FRF proxy' if normalised else 'Absolute FRF proxy'
    fig.suptitle(f'{kind} — {name}', x=0.075, ha='left', fontsize=16,
                 fontweight='bold')
    fig.text(0.075, 0.865,
             f'{len(runs)} detected HIGH-profile run'
             f'{"s" if len(runs) != 1 else ""} · three accelerometer axes · '
             '250–1000 Hz',
             ha='left', fontsize=9, color=vs.INK_SOFT)

    direction_style = {
        'up': (vs.FORWARD, 'up sweep'),
        'down': (vs.REVERSE, 'down sweep'),
    }
    run_styles = ('-', '--', ':', '-.')
    handles = []
    labels = []
    all_values = []

    for axis_number, (axis_name, ax) in enumerate(zip(AXES, axes), 1):
        for curve in curves:
            if curve['axis'] != axis_name:
                continue
            frequency, response = finite_positive(curve)
            if not response.size:
                continue
            colour, direction_label = direction_style[curve['direction']]
            linestyle = run_styles[(curve['run'] - 1) % len(run_styles)]
            if normalised:
                reference = float(np.nanmedian(response))
                values = 20.0 * np.log10(response / reference)
            else:
                values = response
            all_values.extend(values[np.isfinite(values)].tolist())
            line, = ax.plot(frequency, values, color=colour,
                            linestyle=linestyle, linewidth=1.15, alpha=0.9)
            label = f'run {curve["run"]}, {direction_label}'
            if axis_number == 1 and label not in labels:
                handles.append(line)
                labels.append(label)

        ax.set_xlim(250, 1000)
        ax.set_xlabel('drive frequency [Hz]')
        ax.set_title(axis_name)
        if normalised:
            ax.axhline(0, color=vs.INK_MUTED, linewidth=0.75, zorder=0)
        else:
            ax.set_yscale('log')
        vs.despine(ax)

    if normalised:
        axes[0].set_ylabel('response / own sweep median [dB]')
        finite = np.asarray(all_values, dtype=float)
        if finite.size:
            lo, hi = np.nanpercentile(finite, [0.5, 99.5])
            pad = max(2.0, 0.06 * (hi - lo))
            for ax in axes:
                ax.set_ylim(np.floor(lo - pad), np.ceil(hi + pad))
        footer = ('Shape-only comparison. Every sweep and axis is divided by '
                  'its own median response over 250–1000 Hz; 0 dB is that '
                  'median. Normalisation does not reconstruct the missing '
                  'input measurement.')
    else:
        axes[0].set_ylabel('acceleration at drive frequency [m/s² RMS]')
        footer = ('Unnormalised drive-frequency response extracted from the '
                  'measured acceleration. No excitation channel was recorded, '
                  'so this is an FRF proxy—not acceleration divided by a '
                  'measured input.')

    fig.legend(handles, labels, loc='upper right', bbox_to_anchor=(0.985, 0.89),
               ncol=max(1, min(2, len(runs))), frameon=False, fontsize=8)
    fig.text(0.075, 0.07, footer, ha='left', va='bottom', fontsize=8,
             color=vs.INK_MUTED, wrap=True)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--results', type=Path,
                        default=HERE / 'chirp_results.json')
    parser.add_argument('--out', type=Path,
                        default=HERE / 'plots' / 'chirp_frf')
    args = parser.parse_args()

    vs.use()
    payload = json.loads(args.results.read_text(encoding='utf-8'))
    args.out.mkdir(parents=True, exist_ok=True)
    pdf_path = args.out / 'chirp_frf_raw_and_normalised.pdf'

    written = []
    with PdfPages(pdf_path) as pdf:
        for recording in payload.get('recordings', []):
            name = recording['recording']
            for normalised, suffix in ((False, 'raw'),
                                       (True, 'normalised')):
                fig = page(recording, normalised=normalised)
                path = args.out / f'chirp_frf_{name}_{suffix}.png'
                fig.savefig(path, dpi=200)
                pdf.savefig(fig, dpi=200)
                plt.close(fig)
                written.append(path)
                print(f'wrote {path.name}')

    print(f'wrote {pdf_path.name} ({len(written)} pages)')


if __name__ == '__main__':
    main()
