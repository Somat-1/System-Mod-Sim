#!/usr/bin/env python3
'''Process the chirp accelerometer recordings into results and figures.

Neither recording carries a timestamp tying it to the planned profile, so the
profile each one contains is identified from the excitation itself: the
instantaneous frequency of the drive is tracked through a spectrogram, the
sweep legs are found as the intervals where that track moves monotonically at
the planned rate, and each leg is then matched against the profiles in
`chirp/sequence/chirp_plan.csv` by band and duration.

What comes out per sweep leg is the response at the frequency being driven -
the acceleration in a narrow band that follows the sweep - which is the curve
the campaign is after. Up and down legs are kept apart throughout: a sweep
that agrees with itself in both directions is measuring the structure, and one
that does not is measuring the sweep.
'''

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy import ndimage, signal

import matplotlib.pyplot as plt

import vizstyle as vs

HERE = Path(__file__).resolve().parent
STAGE_DIR = HERE.parent
NPERSEG = 4096
NOVERLAP = 3072
TRACK_FLOOR_HZ = 40.0        # below this the drive is lost in the noise floor
SWEEP_RATE_TOL = 0.45        # fraction of the planned rate a leg may differ by
AXES = ('axis 1', 'axis 2', 'axis 3')


def load_profiles(path):
    rows = list(csv.DictReader(path.open(encoding='utf-8')))
    profiles = {}
    for r in rows:
        p = profiles.setdefault(r['profile'], [])
        row = {k: (float(v) if k in ('t_start', 't_end', 'duration_s',
                                     'f_start', 'f_end',
                                     'amplitude_microsteps') else v)
               for k, v in r.items()}
        p.append(row)
    out = {}
    for name, items in profiles.items():
        sweeps = [i for i in items if i['kind'] == 'sweep']
        out[name] = dict(
            rows=items, sweeps=sweeps,
            duration=max(i['t_end'] for i in items),
            f_lo=min(min(i['f_start'], i['f_end']) for i in sweeps),
            f_hi=max(max(i['f_start'], i['f_end']) for i in sweeps),
            rate=float(np.median([abs(i['f_end'] - i['f_start']) / i['duration_s']
                                  for i in sweeps if i['duration_s'] > 0])))
    return out


# --------------------------------------------------------------------------
# finding the excitation
# --------------------------------------------------------------------------

def spectro(x, fs):
    f, t, S = signal.spectrogram(x, fs=fs, nperseg=NPERSEG, noverlap=NOVERLAP,
                                 window='hann')
    level = 10 * np.log10(S + 1e-16)
    # How far each bin stands above its own spectral neighbourhood, which is
    # what a driven tone does. Normalising against the whole slice instead
    # would just measure the 1/f tilt, and every low-frequency bin would look
    # like a drive whether anything was driving it or not.
    relief = level - ndimage.median_filter(level, size=(41, 1), mode='nearest')
    return dict(f=f, t=t, S=S, level=level, relief=relief)


def find_sweeps(sp, sweeps, min_relief=8.0):
    '''Locate each planned sweep leg by matched line search.

    A peak-picking tracker does not survive this rig: the structure has a loud
    band near 680 Hz that outweighs the drive wherever the sweep is not passing
    through it, so the peak bin follows the resonance instead of the
    excitation. The drive is found instead by sliding each planned leg - whose
    start frequency, end frequency and duration are all known - along the
    recording and scoring the spectrogram under that line.
    '''
    f, t, relief = sp['f'], sp['t'], sp['relief']
    dt = float(np.median(np.diff(t)))
    df = float(f[1] - f[0])
    found = []
    for leg in sweeps:
        n = int(round(leg['duration_s'] / dt))
        if n < 10 or n >= t.size:
            continue
        ramp = np.linspace(leg['f_start'], leg['f_end'], n)
        rows = np.clip(np.rint(ramp / df).astype(int), 0, f.size - 1)
        cols = np.arange(n)
        score = np.array([relief[rows, cols + i0].mean()
                          for i0 in range(t.size - n)])
        # Keep every start that is a clear local best, not just the global one:
        # a recording may hold the same profile more than once.
        order = np.argsort(score)[::-1]
        taken = []
        for i0 in order:
            if score[i0] < min_relief:
                break
            if any(abs(i0 - j) < n * 0.8 for j in taken):
                continue
            taken.append(int(i0))
        for i0 in taken:
            found.append(dict(
                i0=i0, i1=i0 + n, t0=float(t[i0]), t1=float(t[i0 + n - 1]),
                f0=leg['f_start'], f1=leg['f_end'],
                direction='up' if leg['f_end'] > leg['f_start'] else 'down',
                label=leg['label'], profile_hint=leg['profile'],
                seconds=float(leg['duration_s']),
                rate_hz_per_s=float((leg['f_end'] - leg['f_start'])
                                    / leg['duration_s']),
                relief_db=float(score[i0])))
    found.sort(key=lambda d: d['t0'])

    # Two planned legs can describe the same stretch of recording (the high leg
    # of FULL and of HIGH are identical); keep the better-scoring one.
    kept = []
    for leg in found:
        clash = next((k for k in kept
                      if leg['t0'] < k['t1'] - 2 and k['t0'] < leg['t1'] - 2), None)
        if clash is None:
            kept.append(leg)
        elif leg['relief_db'] > clash['relief_db']:
            kept[kept.index(clash)] = leg
    return kept


def group_runs(legs, gap_s=40.0):
    '''An up leg and the down leg that follows it are one profile run.'''
    runs, current = [], []
    for leg in legs:
        if current and leg['t0'] - current[-1]['t1'] > gap_s:
            runs.append(current)
            current = []
        current.append(leg)
    if current:
        runs.append(current)
    return runs


def identify(run, profiles):
    '''Which planned profile this run is, from the legs that were found.

    A profile is only called matched when every sweep leg it specifies is
    present. That matters here: the high leg of FULL and the whole of HIGH
    cover the same band, and the thing that tells them apart is whether the
    1-150 Hz legs are in the recording as well.
    '''
    lo = min(min(l['f0'], l['f1']) for l in run)
    hi = max(max(l['f0'], l['f1']) for l in run)
    # Identify on the bands that were swept, not on the labels: the 250-1000 Hz
    # leg appears in both profiles, so a label would not separate them. What
    # separates them is whether the 1-150 Hz legs are there as well.
    found = {(l['f0'], l['f1']) for l in run}
    best, score = None, -1.0
    for name, p in profiles.items():
        want = {(s['f_start'], s['f_end']) for s in p['sweeps']}
        frac = len(found & want) / len(want) - 0.25 * len(found - want)
        if frac > score:
            best, score = name, frac
    want = {(s['f_start'], s['f_end']) for s in profiles[best]['sweeps']}
    missing = sorted(f'{a:g}-{b:g} Hz' for a, b in want - found)
    labels = sorted(f'{a:g}-{b:g} Hz' for a, b in found)
    return best, float(score), dict(f_lo=float(lo), f_hi=float(hi),
                                    swept_seconds=float(sum(l['seconds']
                                                            for l in run)),
                                    legs_found=labels, legs_missing=missing,
                                    complete=not missing)


# --------------------------------------------------------------------------
# response along the sweep
# --------------------------------------------------------------------------

def response(accel, fs, leg, sp, n_bins=180):
    '''Response at the driven frequency, binned on that frequency.

    The sweep passes each frequency once, so the response at a frequency is
    read out of the slice of the spectrogram that was being driven there, in a
    band around the drive rather than at a single bin: the drive is not exactly
    sinusoidal and the bin alone would throw away the skirt.
    '''
    f, t = sp['f'], sp['t']
    # Each leg is bracketed by a 3 s dwell at 250 Hz. Drop half a second at
    # each end so that dwell cannot leak through the analysis window and be
    # read as a response at the frequency the sweep starts from.
    guard = int(round(0.5 / float(np.median(np.diff(t)))))
    i0, i1 = leg['i0'] + guard, leg['i1'] - guard
    sl = slice(i0, i1)
    span = (leg['f1'] - leg['f0']) / (leg['i1'] - leg['i0'])
    drive = leg['f0'] + span * np.arange(i0 - leg['i0'], i1 - leg['i0'])
    order = np.argsort(drive)
    edges = np.linspace(drive.min(), drive.max(), n_bins + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])

    out = {}
    for ax_i in range(accel.shape[1]):
        _, _, S = signal.spectrogram(accel[:, ax_i].astype(np.float64), fs=fs,
                                     nperseg=NPERSEG, noverlap=NOVERLAP,
                                     window='hann')
        S = S[:, sl]
        amp = np.empty(drive.size)
        for k in range(drive.size):
            lo, hi = drive[k] * 0.93, drive[k] * 1.07
            band = (f >= lo) & (f <= hi)
            amp[k] = np.sqrt(S[band, k].sum())
        binned = np.full(n_bins, np.nan)
        idx = np.clip(np.digitize(drive[order], edges) - 1, 0, n_bins - 1)
        for b in range(n_bins):
            vals = amp[order][idx == b]
            if vals.size:
                binned[b] = float(np.median(vals))
        out[AXES[ax_i]] = binned
    del t
    return centres, out


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

def fig_spectrogram(name, sp, runs, idents, profiles, out):
    fig, ax = plt.subplots(figsize=(12, 4.4), constrained_layout=True)
    f, t, S = sp['f'], sp['t'], sp['S']
    keep = f <= 1150
    ax.pcolormesh(t, f[keep], 10 * np.log10(S[keep] + 1e-14), shading='auto',
                  cmap='magma', vmin=-110, vmax=-45, rasterized=True)
    for run, (prof, _, meas) in zip(runs, idents):
        for leg in run:
            ax.plot([leg['t0'], leg['t1']], [leg['f0'], leg['f1']],
                    color='#ffffff', lw=0.9, ls='--', alpha=0.85)
        mid = 0.5 * (run[0]['t0'] + run[-1]['t1'])
        ax.annotate(f'{prof.upper()}', (mid, 1090), ha='center', va='top',
                    color='#ffffff', fontsize=9, fontweight='bold')
    for edge in (150, 250):
        ax.axhline(edge, color='#ffffff', lw=0.7, alpha=0.5)
    ax.annotate('150-250 Hz: never driven, by design', (t[-1] * 0.02, 200),
                ha='left', va='center', color='#ffffff', fontsize=7.5)
    ax.set_xlabel('recording time (s)')
    ax.set_ylabel('frequency (Hz)')
    ax.set_title(f'{name}: what is actually in the recording')
    ax.set_ylim(0, 1150)
    # No SVG for this one: the mesh is an image, and the vector file would be
    # the same pixels in a wrapper eight times the size.
    return vs.save(fig, out, f'chirp_01_spectrogram_{name}')


def fig_inventory(inventory, curves, profiles, out):
    '''What was recorded against what was planned, and whether runs agree.'''
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(12, 4.2),
                                  constrained_layout=True)
    for row, rec in enumerate(inventory):
        y = len(inventory) - 1 - row
        ax.barh(y, rec['seconds'], height=0.46, color=vs.GRID,
                edgecolor=vs.SURFACE, lw=1.2)
        for n, run in enumerate(rec['runs']):
            for leg in run['legs']:
                ax.barh(y, leg['t1'] - leg['t0'], left=leg['t0'], height=0.46,
                        color=vs.FORWARD if leg['direction'] == 'up'
                        else vs.REVERSE, edgecolor=vs.SURFACE, lw=1.0)
            mid = 0.5 * (run['legs'][0]['t0'] + run['legs'][-1]['t1'])
            ax.annotate(run['profile'].upper(), (mid, y + 0.27), ha='center',
                        va='bottom', fontsize=8, color=vs.INK_SOFT,
                        fontweight='bold')
    ax.set_yticks(range(len(inventory)),
                  [r['recording'] for r in inventory][::-1])
    ax.set_ylim(-0.75, len(inventory) - 0.4)
    ax.set_xlabel('recording time (s)')
    ax.set_title('What each recording contains')
    vs.note(ax, 'blue: up sweep   orange: down sweep   grey: recorded, '
                'nothing driven\n'
                f"the planned FULL profile ({profiles['full']['duration']:.0f} s, "
                'down to 1 Hz) is in neither file', 'lower left')
    vs.despine(ax)

    runs = []
    for c in curves:
        if c['run'] not in runs:
            runs.append(c['run'])
    grid = None
    shapes, levels = [], []
    for key in runs:
        grid, one, _, _, _ = _stack([c for c in curves if c['run'] == key],
                                    'up', 'axis 1')
        level = float(np.nanmedian(one))
        levels.append(level)
        shapes.append(one / level)
    # Geometric mean, not the median: with three runs the median is one of
    # them, and that run would be pinned to zero over much of the axis.
    reference = np.exp(np.nanmean(np.log(np.vstack(shapes)), axis=0))
    for n, (key, shape) in enumerate(zip(runs[:3], shapes[:3])):
        ax2.plot(grid, 20 * np.log10(shape / reference), color=vs.SERIES[n],
                 lw=1.3, label=key)
    ax2.axhline(0, color=vs.INK_MUTED, lw=0.8)
    ax2.set_xlim(250, 1000)
    ax2.set_ylim(-12, 12)
    ax2.set_xlabel('drive frequency (Hz)')
    ax2.set_ylabel('departure from the common shape (dB)')
    ax2.set_title('Do the three runs agree, once level is removed?')
    ax2.legend(loc='upper right')
    spread = 20 * np.log10(max(levels) / min(levels))
    vs.note(ax2, 'each run divided by its own median level first:\n'
                 f'the overall levels themselves differ by {spread:.0f} dB',
            'lower left')
    vs.despine(ax2)
    return vs.save(fig, out, 'chirp_02_inventory')


def level_reference(x, fs, run, sp):
    '''Measure the two fixed-amplitude excitations every HIGH run contains.

    The sync markers (3 bursts at 750 Hz, 32 microsteps) and the 2 s dwell at
    1 kHz at the turnaround (48 microsteps) are specified identically for every
    run. Whatever they read, they should read the same each time. Where they do
    not, the difference is in the drive or the sensing, not in the stage, and
    nothing about response amplitude can be compared across those runs.
    '''
    out = {}
    sos = signal.butter(4, [730, 770], 'bandpass', fs=fs, output='sos')
    for tag, centre in (('marker_start', run[0]['t0'] - 3.5),
                        ('marker_end', run[-1]['t1'] + 3.0)):
        a, b = int((centre - 2.0) * fs), int((centre + 2.0) * fs)
        if a < 0 or b > x.size:
            out[tag] = None
            continue
        out[tag] = float(np.sqrt(np.mean(signal.sosfilt(sos, x[a:b]) ** 2)))
    # The gap between the up and down legs is the planned 2 s mid dwell. It is
    # recorded because it is a structural fingerprint of the run: two files
    # whose turnaround differs did not come off the same configuration,
    # whatever their response curves look like.
    out['turnaround_s'] = float(run[-1]['t0'] - run[0]['t1'])
    del sp
    return out


def _stack(curves, direction, axis):
    '''Median and spread across every leg going one way, on a common grid.'''
    picked = [c for c in curves if c['direction'] == direction]
    if not picked:
        return None, None, None, None, 0
    grid = picked[0]['freq']
    rows = np.vstack([np.interp(grid, c['freq'],
                                np.array([np.nan if v is None else v
                                          for v in c['response'][axis]], float))
                      for c in picked])
    return (grid, np.nanmedian(rows, axis=0), np.nanmin(rows, axis=0),
            np.nanmax(rows, axis=0), rows.shape[0])


def fig_response(curves, out):
    """Response against drive frequency: one panel per accelerometer axis.

    Two series, up and down, each the median over the repeats with the spread
    across them as a band. Drawing all six legs separately would put six
    near-identical curves on one axis and say less than the band does.
    """
    fig, axes = plt.subplots(1, 3, figsize=(13.5, 4.3), constrained_layout=True,
                             sharey=True)
    for ax_i, axis in enumerate(AXES):
        ax = axes[ax_i]
        for direction, colour in (('up', vs.FORWARD), ('down', vs.REVERSE)):
            grid, med, lo, hi, n = _stack(curves, direction, axis)
            if med is None:
                continue
            ax.fill_between(grid, lo, hi, color=colour, alpha=0.18, lw=0)
            ax.semilogy(grid, med, color=colour, lw=1.3,
                        label=f'{direction} sweep (n={n})' if ax_i == 0 else None)
        ax.set_xlabel('drive frequency (Hz)')
        ax.set_title(f'Response, {axis}')
        ax.set_xlim(250, 1000)
        if ax_i == 0:
            ax.set_ylabel('acceleration at the drive (m/s2, rms)')
            ax.legend(loc='lower right')
            vs.note(ax, 'band spans the repeats')
        vs.despine(ax)
    return vs.save(fig, out, 'chirp_03_response')


def fig_updown(curves, out):
    """The 250 Hz question, and whether the repeats agree with each other."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(11.5, 4.3),
                                  constrained_layout=True)
    runs = []
    for c in curves:
        if c['run'] not in runs:
            runs.append(c['run'])
    for n, key in enumerate(runs[:3]):
        grid, med, _, _, _ = _stack([c for c in curves if c['run'] == key],
                                    'up', 'axis 1')
        if med is None:
            continue
        ax.semilogy(grid, med, color=vs.SERIES[n], lw=1.3, label=key)
    ax.axvline(250, color=vs.INK_MUTED, lw=1.0, ls=':')
    ax.set_xlim(250, 420)
    ax.set_xlabel('drive frequency (Hz)')
    ax.set_ylabel('acceleration at the drive (m/s2, rms)')
    ax.set_title('Around 250 Hz, where the sweep starts')
    ax.legend(loc='upper right')
    vs.note(ax, 'up sweeps only, one per run')
    vs.despine(ax)

    # Each run compared against itself. Pooling every up leg and dividing by
    # the pooled down legs would mix runs whose overall levels differ, and the
    # result would describe that mismatch rather than the sweep direction.
    for n, key in enumerate(runs[:3]):
        mine = [c for c in curves if c['run'] == key]
        grid, up, _, _, _ = _stack(mine, 'up', 'axis 1')
        _, down, _, _, _ = _stack(mine, 'down', 'axis 1')
        if up is None or down is None:
            continue
        ax2.plot(grid, 20 * np.log10(up / down), color=vs.SERIES[n], lw=1.3,
                 label=key)
    ax2.axhline(0, color=vs.INK_MUTED, lw=0.8)
    ax2.set_xlim(250, 1000)
    ax2.set_ylim(-12, 12)
    ax2.set_xlabel('drive frequency (Hz)')
    ax2.set_ylabel('up minus down, within one run (dB)')
    ax2.set_title('Direction difference: each run against itself')
    ax2.legend(loc='upper right')
    vs.note(ax2, 'a structural resonance sits at the same frequency\n'
                 'both ways; a sweep-rate artefact leans with the direction',
            'lower right')
    vs.despine(ax2)
    return vs.save(fig, out, 'chirp_04_direction')


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, default=HERE / 'cache')
    ap.add_argument('--plan', type=Path,
                    default=STAGE_DIR / 'chirp' / 'sequence' / 'chirp_plan.csv')
    ap.add_argument('--out', type=Path, default=HERE / 'plots')
    ap.add_argument('--recordings', nargs='*',
                    default=['tempTvalm', 'tempTvalm2'])
    args = ap.parse_args()

    vs.use()
    profiles = load_profiles(args.plan)
    all_sweeps = [s for p in profiles.values() for s in p['sweeps']]
    print('planned profiles: ' + ', '.join(
        f"{k} ({v['f_lo']:g}-{v['f_hi']:g} Hz, {v['duration']:.0f} s)"
        for k, v in profiles.items()))

    results, curves, all_runs = [], [], []
    for name in args.recordings:
        path = args.cache / f'{name}.npz'
        if not path.exists():
            print(f'{name}: no cache, skipped')
            continue
        z = np.load(path)
        accel, dt = z['accel'], float(z['dt'])
        fs = 1.0 / dt
        print(f'\n{name}: {accel.shape[0]} samples x {accel.shape[1]} axes @ '
              f'{fs / 1000:.0f} kHz = {accel.shape[0] * dt:.1f} s')

        sp = spectro(accel[:, 0].astype(np.float64), fs)
        legs = find_sweeps(sp, all_sweeps)
        runs = group_runs(legs)
        idents = [identify(run, profiles) for run in runs]
        print(f'  {len(legs)} sweep legs in {len(runs)} run(s)')

        rec = dict(recording=name, seconds=float(accel.shape[0] * dt),
                   sample_hz=fs, axes=accel.shape[1], runs=[])
        for k, (run, (prof, score, meas)) in enumerate(zip(runs, idents), 1):
            label = f'{name} run {k}' if len(runs) > 1 else name
            print(f'  run {k}: {prof.upper()} (match {score:.3f}) '
                  f'{meas["f_lo"]:.0f}-{meas["f_hi"]:.0f} Hz, '
                  f'{run[0]["t0"]:.1f}-{run[-1]["t1"]:.1f} s, '
                  f'{len(run)} legs')
            all_runs.append((label, run, prof))
            ref = level_reference(accel[:, 0].astype(np.float64), fs, run, sp)
            run_rec = dict(label=label, profile=prof, match=score, legs=[],
                           level_reference=ref, **meas)
            print('     fixed-amplitude references: ' + ', '.join(
                f'{k} {v:.5f}' if v is not None else f'{k} missing'
                for k, v in ref.items()))
            for leg in run:
                freq, per_axis = response(accel, fs, leg, sp)
                curves.append(dict(run=label, direction=leg['direction'],
                                   freq=freq, response=per_axis))
                peak_axis = {a: dict(
                    peak_hz=float(freq[np.nanargmax(v)]),
                    peak_value=float(np.nanmax(v)),
                    median=float(np.nanmedian(v))) for a, v in per_axis.items()}
                run_rec['legs'].append(dict(
                    direction=leg['direction'], t0=leg['t0'], t1=leg['t1'],
                    f0=leg['f0'], f1=leg['f1'], seconds=leg['seconds'],
                    rate_hz_per_s=leg['rate_hz_per_s'], peaks=peak_axis,
                    frequency=freq.tolist(),
                    response={a: np.where(np.isfinite(v), v, None).tolist()
                              for a, v in per_axis.items()}))
                print(f'     {leg["direction"]:>4}: {leg["f0"]:6.1f} -> '
                      f'{leg["f1"]:6.1f} Hz in {leg["seconds"]:5.1f} s '
                      f'({leg["rate_hz_per_s"]:+.2f} Hz/s), '
                      f'peak axis 1 at {peak_axis["axis 1"]["peak_hz"]:.0f} Hz')
            rec['runs'].append(run_rec)
        results.append(rec)

        print('  figures:')
        fig_spectrogram(name, sp, runs, idents, profiles, args.out)

    if all_runs:
        print('figures:')
        fig_inventory(results, curves, profiles, args.out)
        fig_response(curves, args.out)
        fig_updown(curves, args.out)

    out = HERE / 'chirp_results.json'
    out.write_text(json.dumps(dict(
        planned={k: {kk: vv for kk, vv in v.items() if kk not in ('rows', 'sweeps')}
                 for k, v in profiles.items()},
        recordings=results), indent=1), encoding='utf-8')
    print(f'wrote {out.name}')


if __name__ == '__main__':
    main()
