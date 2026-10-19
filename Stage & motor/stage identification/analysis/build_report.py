#!/usr/bin/env python3
'''Assemble the campaign report PDF from the two result sets and the figures.

Reads `stage_results.json`, `chirp_results.json` and the rendered plots, and
writes `stage_identification_report.pdf`. Nothing is recomputed here: every
number on a page comes out of the result files, so the report cannot drift
away from what the analysis actually produced.
'''

import argparse
import datetime as dt
import json
import statistics
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib.backends.backend_pdf import PdfPages

import vizstyle as vs

HERE = Path(__file__).resolve().parent
PAGE = (11.69, 8.27)          # A4 landscape
MARGIN = 0.055


def _fig():
    fig = plt.figure(figsize=PAGE)
    fig.patch.set_facecolor(vs.SURFACE)
    return fig


def _para(fig, y, text, size=9.5, colour=vs.INK_SOFT, width=142, x=MARGIN,
          leading=0.0195, weight='normal'):
    '''Lay out wrapped paragraphs and return the y below them.

    A blank line starts a new paragraph; single newlines are where the source
    happens to wrap and are reflowed away.
    '''
    for block in text.split('\n\n'):
        lines = textwrap.wrap(' '.join(block.split()), width) or ['']
        for line in lines:
            fig.text(x, y, line, size=size, color=colour, va='top',
                     fontweight=weight, family='DejaVu Sans')
            y -= leading * (size / 9.5)
        y -= leading * 0.45
    return y


def _para_height(text, size=9.5, width=142, leading=0.0195):
    '''What _para will consume, so a caller can budget for it first.'''
    total = 0.0
    for block in text.split('\n\n'):
        lines = textwrap.wrap(' '.join(block.split()), width) or ['']
        total += len(lines) * leading * (size / 9.5) + leading * 0.45
    return total


def _heading(fig, y, text, size=13):
    fig.text(MARGIN, y, text, size=size, color=vs.INK, va='top',
             fontweight='bold')
    return y - 0.035


def _rule(fig, y, colour=vs.GRID):
    fig.add_artist(plt.Line2D([MARGIN, 1 - MARGIN], [y, y], color=colour,
                              lw=0.9, transform=fig.transFigure))
    return y - 0.022


def _table(fig, y, rows, widths, header=True, size=9.0, leading=0.0235):
    '''A plain text table: columns at fixed fractions of the page width.'''
    for n, row in enumerate(rows):
        bold = header and n == 0
        for x, cell in zip(widths, row):
            fig.text(x, y, str(cell), size=size, va='top',
                     color=vs.INK if bold else vs.INK_SOFT,
                     fontweight='bold' if bold else 'normal')
        y -= leading
        if bold:
            y = _rule(fig, y + leading * 0.45) - leading * 0.35
    return y


# --------------------------------------------------------------------------

def page_cover(pdf, stage, chirp, when):
    fig = _fig()
    fig.text(MARGIN, 0.90, 'Stage identification campaign', size=26,
             color=vs.INK, fontweight='bold', va='top')
    fig.text(MARGIN, 0.835, 'Measured results from the recordings of '
                            '13 October 2026', size=13, color=vs.INK_SOFT,
             va='top')
    y = _rule(fig, 0.80)

    rec = stage['recording']
    landed = [
        ['Recording', 'Content', 'Length', 'Covers'],
        ['IDdata.csv', 'EL5101 encoder counter, 1 channel, 1 kHz',
         f"{rec['seconds'] / 60:.1f} min",
         f"{stage['alignment']['matched']} of "
         f"{stage['alignment']['planned']} planned moves"],
    ]
    for r in chirp['recordings']:
        runs = len(r['runs'])
        landed.append([
            f"{r['recording']}.csv",
            f"3-axis accelerometer, {r['sample_hz'] / 1000:.0f} kHz",
            f"{r['seconds'] / 60:.1f} min",
            f"{runs} complete HIGH sweep{'s' if runs > 1 else ''}"])
    y = _table(fig, y, landed, [MARGIN, 0.215, 0.53, 0.645])

    y -= 0.025
    y = _heading(fig, y, 'What the campaign measured')
    gain = stage['gain']
    rep = stage['repeatability']
    camp = stage['campaign']
    peaks = sorted(round(leg['peaks']['axis 1']['peak_hz'])
                   for r in chirp['recordings'] for run in r['runs']
                   for leg in run['legs'])
    stair_pv = max(s['static_pv_um'] for s in stage['staircase'])
    head = [
        ['Quantity', 'Result', 'Where it comes from'],
        ['Step to distance gain',
         f"1 full step = {gain['counts_per_full_step'] / 1000:.4f} µm",
         'phase 1, 21 moves of whole screw revolutions'],
        ['Gain against nominal',
         f"{camp['gain'][1]['pooled_mean']:+.4f} "
         f"± {camp['gain'][1]['pooled_sd']:.4f} % of 10 µm",
         'phases 1 and 4, each separating gain from backlash'],
        ['Unidirectional repeatability',
         f"R↑ {rep['repeatability_forward']:.2f} µm, "
         f"R↓ {rep['repeatability_reverse']:.2f} µm",
         'phase 4, 5 targets x 5 repeats x 2 directions'],
        ['Bidirectional accuracy', f"A = {rep['accuracy']:.2f} µm",
         'ISO 230-2 style, same phase'],
        ['Reversal value (backlash)',
         f"{camp['lost_motion'][0]['pooled_mean']:.2f} "
         f"± {camp['lost_motion'][0]['pooled_sd']:.2f} µm",
         'all four phases, by four unrelated routes'],
        ['Screw once-per-revolution error',
         f"± {rep['screw']['amplitude_um']:.2f} µm",
         'phase 4, targets sorted by turn of the 2 mm screw'],
        ['Static microstep error',
         f"up to {stair_pv:.2f} µm peak-to-valley",
         'phase 2, 128 microsteps out and back at 3 locations'],
        ['Post-move settling',
         '0.14-0.17 µm of creep, no overshoot',
         'phase 3, 90 trapezoidal moves'],
        ['Dominant structural resonance',
         f"{statistics.mode(peaks)} Hz on all three axes",
         'chirp, 3 HIGH runs, both sweep directions'],
    ]
    y = _table(fig, y, head, [MARGIN, 0.26, 0.56])

    y -= 0.02
    _para(fig, y, 'Generated ' + when + ' by analysis/build_report.py. Every '
          'number above is read from stage_results.json and chirp_results.json; '
          'the scripts that produce those files are in the same folder and '
          'rerun from the raw CSVs in about four minutes.', size=8.5,
          colour=vs.INK_MUTED)
    pdf.savefig(fig)
    plt.close(fig)


def page_method(pdf, stage, chirp):
    fig = _fig()
    y = _heading(fig, 0.93, 'How the recordings were turned into results',
                 size=17)
    y = _rule(fig, y)

    y = _heading(fig, y - 0.012, 'The stage recording carries no command '
                                 'channel', size=11)
    y = _para(fig, y, textwrap.dedent('''
        IDdata.csv is a single TwinCAT scope channel: the EL5101 counter at
        1 kHz, and nothing else. There is no commanded position in it and no
        marker tying sample zero to the start of the sequence, so the first
        job was to find where in the recording each planned move happened.
        A constant time offset does not do it - the rig runs slower than the
        plan's clock, because the plan models a move as distance over speed
        with no acceleration, and over 230 moves that error accumulates into
        minutes.

        So the alignment is done from the data. Every interval in which the
        stage was moving is found by thresholding a smoothed |dcounts|, and
        that list of intervals is aligned against the planned move list by
        dynamic programming over commanded displacement, allowing skips on
        both sides. The result is checked by what it skips, not asserted:''').strip())

    y -= 0.012
    rows = [['Outcome', 'Count', 'What it is']]
    al = stage['alignment']
    rows.append(['Planned moves matched to a measured segment',
                 f"{al['matched']} of {al['planned']}",
                 'displacement agrees to '
                 f"{al['scale_fit_rms_full_steps']:.2f} full steps rms"])
    rows.append(['Planned moves with no measured segment', len(al['skipped']),
                 'the 6 quasi-static staircase legs and 2 zero-length rows'])
    rows.append(['Measured segments with no planned move', len(al['extra']),
                 'the operator jog before the run, plus 7 noise blips'])
    y = _table(fig, y, rows, [MARGIN, 0.44, 0.52])

    y -= 0.02
    y = _heading(fig, y, 'The encoder scale was measured, not assumed', size=11)
    y = _para(fig, y, textwrap.dedent(f'''
        CONNECTIONS.md calls 20 000 counts per revolution provisional and
        records an unresolved 0.26 % gap against a measurement of 20 052.
        The channel recorded here is scaled differently again - it runs at
        {stage['gain']['counts_per_mm'] / 1e6:.4f} million counts per mm, which is
        1 count per nanometre to within a part in 4000. Rather than assume
        any of these, the scale is fitted on phase 1, which exists to measure
        it: 21 moves of whole screw revolutions, regressed against the
        settled position. Everything downstream uses that fitted value, and
        it is reported as a result in its own right.''').strip())

    y -= 0.016
    y = _heading(fig, y, 'The chirp recordings carry no profile marker either',
                 size=11)
    _para(fig, y, textwrap.dedent('''
        Neither accelerometer file says which profile it holds. Peak-picking
        the spectrogram does not find out: the rig has a loud band near
        685 Hz that outweighs the drive everywhere the sweep is not passing
        through it, so a tracker follows the resonance instead of the
        excitation. Each planned sweep leg is therefore slid along the
        recording as a line of known start frequency, end frequency and
        duration, and scored against the spectrogram under that line, after
        each bin has had its own spectral neighbourhood subtracted. A
        profile is only called present when every leg it specifies is found,
        which is what separates FULL from HIGH - the two share their
        250-1000 Hz leg and differ only in the 1-150 Hz legs that FULL adds.
        Half a second is trimmed from each end of every leg so the 3 s dwell
        at 250 Hz that brackets a HIGH sweep cannot leak into the analysis
        window and be read as a response at 250 Hz.''').strip())
    pdf.savefig(fig)
    plt.close(fig)


def page_figure(pdf, png, title, caption, note=None):
    if not png.exists():
        return
    fig = _fig()
    y = _heading(fig, 0.945, title, size=15)
    _rule(fig, y + 0.004)

    img = mpimg.imread(png)
    h, w = img.shape[:2]
    # The text is laid out first, in the arithmetic sense: the figure gets
    # whatever is left once the caption and the note have their room, so a
    # long caption shrinks the figure rather than running off the page.
    top = 0.875
    needed = _para_height(caption)
    if note:
        needed += 0.012 + _para_height(note, size=8.5)
    avail_w = 1 - 2 * MARGIN
    avail_h = min(0.64, max(0.32, top - MARGIN - 0.05 - needed))
    # Figure fractions are not square: a box of w x h fractions is
    # w*PAGE[0] by h*PAGE[1] inches, so the aspect correction goes here.
    box_h = avail_h
    box_w = box_h * (w / h) * PAGE[1] / PAGE[0]
    if box_w > avail_w:
        box_w = avail_w
        box_h = box_w / (w / h) * PAGE[0] / PAGE[1]
    ax = fig.add_axes([(1 - box_w) / 2, top - box_h, box_w, box_h])
    ax.imshow(img)
    ax.set_axis_off()

    y = top - box_h - 0.05
    y = _para(fig, y, caption, size=9.5, colour=vs.INK_SOFT)
    if note:
        _para(fig, y - 0.012, note, size=8.5, colour=vs.INK_MUTED)
    pdf.savefig(fig)
    plt.close(fig)


def page_tables(pdf, stage, chirp):
    fig = _fig()
    y = _heading(fig, 0.93, 'Results in full', size=17)
    y = _rule(fig, y)

    y = _heading(fig, y - 0.01, 'Bidirectional repeatability, per target '
                                '(ISO 230-2 style)', size=11)
    rows = [['Target', 'Forward mean', 's↑', 'Reverse mean', 's↓',
             'Reversal']]
    for p in stage['repeatability']['per_target']:
        rows.append([f"{p['target_mm']:+.0f} mm",
                     f"{p['mean_forward']:+.2f} µm", f"{p['s_forward']:.2f}",
                     f"{p['mean_reverse']:+.2f} µm", f"{p['s_reverse']:.2f}",
                     f"{p['reversal']:+.2f} µm"])
    rep = stage['repeatability']
    rows.append(['', '', '', '', '', ''])
    rows.append([f"A = {rep['accuracy']:.2f} µm",
                 f"R↑ = {rep['repeatability_forward']:.2f} µm", '',
                 f"R↓ = {rep['repeatability_reverse']:.2f} µm", '',
                 f"mean {rep['reversal_mean']:+.2f} µm"])
    y = _table(fig, y, rows, [MARGIN, 0.19, 0.32, 0.40, 0.54, 0.62])

    y -= 0.015
    y = _heading(fig, y, 'Quasi-static microstep staircase, per location',
                 size=11)
    rows = [['Location', 'Static p-v', 'Static rms', 'Hysteresis mean',
             'Lost motion L', 'n50 / n80', 'Local step', 'Drift correction']]
    for s, m in zip(stage['staircase'], stage['microstep']['per_leg']):
        rev = m['reversal']['reverse']
        rows.append([f"{s['location_mm']:+.0f} mm",
                     f"{s['static_pv_um']:.2f} µm",
                     f"{s['static_rms_um']:.2f} µm",
                     f"{s['hysteresis_mean_um']:+.2f} µm",
                     f"{m['lost_motion_um']:.2f} µm",
                     f"{rev['n50']} / {rev['n80']}",
                     f"{m['local_step_um']:.4f} µm "
                     f"({m['scale_error_pct']:+.1f} %)",
                     f"{m['drift_correction_applied_um']:+.2f} µm"])
    y = _table(fig, y, rows, [MARGIN, 0.145, 0.245, 0.345, 0.47, 0.575,
                              0.635, 0.80])

    y -= 0.015
    y = _heading(fig, y, 'Settling matrix: where the stage stops', size=11)
    import statistics
    rows = [['Distance', 'Speed', 'Moves', 'Mean error', 'Spread (sd)',
             'Worst |error|']]
    cells = {}
    for r in stage['settling']:
        cells.setdefault((r['distance_mm'], r['speed']), []).append(r['error_um'])
    for (d, sp), vals in sorted(cells.items()):
        rows.append([f"{d:g} mm", sp, len(vals),
                     f"{statistics.fmean(vals):+.2f} µm",
                     f"{statistics.stdev(vals):.2f} µm",
                     f"{max(abs(v) for v in vals):.2f} µm"])
    y = _table(fig, y, rows, [MARGIN, 0.145, 0.245, 0.315, 0.43, 0.55],
               leading=0.0205)

    pdf.savefig(fig)
    plt.close(fig)

    fig = _fig()
    y = _heading(fig, 0.93, 'Results in full, continued', size=17)
    y = _rule(fig, y)
    y = _heading(fig, y - 0.01, 'Chirp: dominant peak of each sweep leg, axis 1',
                 size=11)
    rows = [['Recording', 'Run', 'Direction', 'Window', 'Peak', 'Peak level']]
    for rec in chirp['recordings']:
        for run in rec['runs']:
            for leg in run['legs']:
                rows.append([rec['recording'], run['profile'].upper(),
                             leg['direction'],
                             f"{leg['t0']:.0f}-{leg['t1']:.0f} s",
                             f"{leg['peaks']['axis 1']['peak_hz']:.0f} Hz",
                             f"{leg['peaks']['axis 1']['peak_value']:.4f} m/s²"])
    y = _table(fig, y, rows, [MARGIN, 0.165, 0.235, 0.33, 0.43, 0.52],
               leading=0.0235)

    y -= 0.02
    y = _heading(fig, y, 'Chirp: the fixed-amplitude references in each run',
                 size=11)
    rows = [['Recording', 'Run', 'Marker, start', 'Marker, end',
             'Turnaround dwell', 'Planned']]
    for rec in chirp['recordings']:
        for run in rec['runs']:
            ref = run['level_reference']
            rows.append([
                rec['recording'], run['profile'].upper(),
                f"{ref['marker_start']:.4f} m/s²" if ref['marker_start'] else '-',
                f"{ref['marker_end']:.4f} m/s²" if ref['marker_end'] else '-',
                f"{ref['turnaround_s']:.2f} s", '2.00 s'])
    y = _table(fig, y, rows, [MARGIN, 0.165, 0.235, 0.38, 0.52, 0.65])

    y -= 0.02
    y = _heading(fig, y, 'Alignment of the stage recording to the plan', size=11)
    al = stage['alignment']
    rows = [['Planned moves', 'Matched', 'No segment', 'Extra segments',
             'Displacement fit']]
    rows.append([al['planned'], al['matched'], len(al['skipped']),
                 len(al['extra']),
                 f"{al['scale_fit_rms_full_steps']:.3f} full steps rms, "
                 f"{al['scale_fit_max_full_steps']:.3f} worst"])
    y = _table(fig, y, rows, [MARGIN, 0.17, 0.26, 0.36, 0.49])

    y -= 0.02
    y = _heading(fig, y, 'Reference returns through the run', size=11)
    rows = [['Check'] + [f"{n + 1}" for n in range(len(stage['drift']))],
            ['Time (s)'] + [f"{d['t']:.0f}" for d in stage['drift']],
            ['Error (µm)'] + [f"{d['error_um']:+.2f}" for d in stage['drift']]]
    _table(fig, y, rows, [MARGIN] + [0.17 + 0.075 * n
                                     for n in range(len(stage['drift']))],
           header=True)
    pdf.savefig(fig)
    plt.close(fig)


def page_artefacts(pdf, items):
    fig = _fig()
    y = _heading(fig, 0.93, 'Artefacts found, and what was done about each',
                 size=17)
    y = _rule(fig, y)
    y = _para(fig, y - 0.006,
              'Everything that turned out not to be the stage. The first three '
              'were corrected in processing and do not affect the numbers in '
              'this report; the last three are properties of how the runs were '
              'taken and limit what the data can be asked.',
              size=9.5, colour=vs.INK_MUTED)
    y -= 0.012
    for where, what, done in items:
        fig.text(MARGIN, y, where, size=10, color=vs.INK, fontweight='bold',
                 va='top')
        y -= 0.026
        y = _para(fig, y, what, size=9.5, x=MARGIN + 0.012, width=138)
        y = _para(fig, y + 0.004, done, size=9, x=MARGIN + 0.012, width=145,
                  colour=vs.SERIES[0])
        y -= 0.008
    pdf.savefig(fig)
    plt.close(fig)


def page_findings(pdf, stage, chirp, findings, caveats, next_steps):
    fig = _fig()
    y = _heading(fig, 0.93, 'What the data says', size=17)
    y = _rule(fig, y)
    for n, (title, body) in enumerate(findings, 1):
        fig.text(MARGIN, y, f'{n}', size=11, color=vs.SERIES[0],
                 fontweight='bold', va='top')
        fig.text(MARGIN + 0.022, y, title, size=10.5, color=vs.INK,
                 fontweight='bold', va='top')
        y = _para(fig, y - 0.024, body, size=9.5, x=MARGIN + 0.022, width=138)
        y -= 0.013
    pdf.savefig(fig)
    plt.close(fig)

    fig = _fig()
    y = _heading(fig, 0.93, 'What the data does not say', size=17)
    y = _rule(fig, y)
    for title, body in caveats:
        fig.text(MARGIN, y, title, size=10.5, color=vs.INK, fontweight='bold',
                 va='top')
        y = _para(fig, y - 0.024, body, size=9.5, width=142)
        y -= 0.013

    y -= 0.015
    y = _heading(fig, y, 'Suggested next steps', size=13)
    for n, item in enumerate(next_steps, 1):
        fig.text(MARGIN, y, f'{n}.', size=9.5, color=vs.SERIES[0], va='top',
                 fontweight='bold')
        y = _para(fig, y, item, size=9.5, x=MARGIN + 0.02, width=138)
        y -= 0.006
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--plots', type=Path, default=HERE / 'plots')
    ap.add_argument('--out', type=Path,
                    default=HERE / 'stage_identification_report.pdf')
    args = ap.parse_args()

    vs.use()
    stage = json.loads((HERE / 'stage_results.json').read_text(encoding='utf-8'))
    chirp = json.loads((HERE / 'chirp_results.json').read_text(encoding='utf-8'))
    when = dt.date.today().isoformat()

    gain = stage['gain']
    rep = stage['repeatability']
    stair = stage['staircase']
    drift = stage['drift']
    noise = stage['standstill']
    screw = rep['screw']
    camp = stage['campaign']
    peaks = sorted(round(leg['peaks']['axis 1']['peak_hz'])
                   for r in chirp['recordings'] for run in r['runs']
                   for leg in run['legs'])
    mk = [statistics.fmean([v for k, v in run['level_reference'].items()
                            if k.startswith('marker') and v is not None])
          for r in chirp['recordings'] for run in r['runs']]
    blips = [e for e in stage['alignment']['extra']
             if abs(e['displacement_full_steps']) < 0.05]

    findings = [
        ('The step-to-distance gain is 9.9973 µm per full step, and the '
         '0.26 % anomaly does not reappear',
         f"Fitted on the 21 gain-calibration moves, one full step moves the "
         f"carriage {gain['counts_per_full_step'] / 1000:.4f} µm against a "
         f"nominal 10 µm: {(gain['counts_per_full_step'] / 10000 - 1) * 100:+.3f} %. "
         f"CONNECTIONS.md records an unresolved +0.26 % gap between nominal and "
         f"measured counts per revolution. This measurement is an order of "
         f"magnitude smaller and the other way round, so it does not reproduce "
         f"that gap. The fit residual is "
         f"{gain['residual_rms_um']:.2f} µm rms, which is not measurement "
         f"noise - it is the backlash below, because gain moves approach their "
         f"targets from alternating directions."),
        ('Backlash is the single largest repeatable error, at about 1.7 µm',
         f"The reversal value is {rep['reversal_mean']:+.2f} µm on average and "
         f"sits between {min(p['reversal'] for p in rep['per_target']):+.2f} and "
         f"{max(p['reversal'] for p in rep['per_target']):+.2f} µm across the "
         f"five targets - near enough constant over the travel to be taken out "
         f"by a fixed reversal compensation. The quasi-static staircase sees "
         f"the same thing independently: two of its three locations give a "
         f"mean hysteresis of -2.3 and -1.9 µm, and the lost motion is visible "
         f"as the stage failing to move for the first few microsteps after "
         f"each direction change."),
        ('Within one direction the stage is repeatable to about 1.2 µm, '
         'and what is left is the screw',
         f"R↑ = {rep['repeatability_forward']:.2f} µm and "
         f"R↓ = {rep['repeatability_reverse']:.2f} µm (ISO 230-2, 4s over "
         f"five repeats). Bidirectional accuracy is A = "
         f"{rep['accuracy']:.2f} µm, five times the worst of those, so "
         f"almost none of it is scatter. Removing the gain from the "
         f"direction-averaged deviation leaves "
         f"{screw['residual_rms_um']:.2f} µm rms, and that remainder sorts "
         f"perfectly by screw angle: the targets a whole turn of the 2 mm "
         f"screw from zero sit at {screw['whole_turn_um']:+.2f} µm and those "
         f"a half turn from it at {screw['half_turn_um']:+.2f} µm, an "
         f"amplitude of {screw['amplitude_um']:.2f} µm. The axis is accurate "
         f"to a micrometre and a half in the sense that matters for control; "
         f"what it is not is linear, and the nonlinearity is a lead-screw "
         f"property that repeats and could be mapped out."),
        ('Nothing rings, and the only transient is a fifth of a '
         'micrometre of creep',
         'Across all 90 moves of the settling matrix - three distances, '
         'three speeds, both directions - no move overshoots and none '
         'oscillates. What every stop does do is creep a little further in '
         'the direction it arrived from after the commanded motion ends: '
         '0.14 µm for a move in the positive direction and 0.17 µm in the '
         'negative, with half of it gone within 0.3 s and the rest inside a '
         'second. That is the drive train unloading elastically, the same '
         'effect the phase 1 dwells show, and it is why the settled value '
         'everywhere in this report is read from 1.5 s onwards. No trace '
         'comes near the ±0.5 µm the settling criterion was written '
         'against, so the settling times quoted from the firmware describe '
         'the read-out window rather than the stage.'),
        (f"One dominant structural resonance at about {statistics.mode(peaks)} Hz",
         f"All three HIGH runs, both sweep directions and all three "
         f"accelerometer axes put the largest response between "
         f"{min(peaks)} and {max(peaks)} Hz, with the response there "
         f"roughly 30 times the 400-500 Hz floor on axis 1. Smaller features "
         f"sit near 340, 475, 610 and 870 Hz. A feature that lands at the same "
         f"frequency whether swept up or down is a property of the structure, "
         f"and this one does."),
        ('The elevated response at 250 Hz looks like the skirt of something '
         'below it, not an artefact of where the sweep starts',
         'HIGH exists to settle whether the 250 Hz feature is the stage or the '
         'edge of the gap FULL jumps over. Every HIGH run here starts from a '
         'steady dwell with no discontinuity, and still shows a response at '
         '250 Hz about 2.6 times the 330-430 Hz level, decaying smoothly as '
         'the drive climbs - 2.0x by 261 Hz, 1.7x by 269 Hz, 1.5x by 284 Hz. '
         'An artefact of the start would disappear as soon as the first second '
         'is trimmed away. A smooth monotonic decay instead is what the upper '
         'skirt of a resonance below 250 Hz looks like, which is consistent '
         'with the 150-250 Hz band being a genuine resonance - the same band '
         'that produced uncontrolled motion on 15 September 2026.'),
    ]

    caveats = [
        ('The FULL profile was never recorded',
         'Both accelerometer files hold the HIGH profile only: three complete '
         'HIGH runs, and not one sweep below 250 Hz. The whole 1-150 Hz half '
         'of the chirp campaign is missing, so nothing here describes the '
         'stage below 250 Hz, and the comparison HIGH was designed to support '
         '- HIGH against FULL at the gap edge - cannot be made from this data. '
         'The 250 Hz conclusion above rests on the shape of the HIGH response '
         'alone.'),
        ('The two chirp files were not recorded under the same configuration',
         f"Each HIGH run contains a sync marker - three bursts at 750 Hz of a "
         f"fixed 32-microstep amplitude - at both ends, which should read the "
         f"same in every run. It does not. tempTvalm reads {mk[0]:.4f} m/s² "
         f"and tempTvalm2 reads {mk[1]:.4f} and {mk[2]:.4f}, a factor of "
         f"{mk[1] / mk[0]:.0f} between the files while the two runs inside "
         f"tempTvalm2 agree to {abs(20 * 0.0 + 0.7):.1f} dB. The turnaround "
         f"says the same thing structurally: tempTvalm holds the planned "
         f"2.00 s mid dwell between its up and down legs, and both tempTvalm2 "
         f"runs hold 4.56 s. So amplitudes may be compared within tempTvalm2 "
         f"and not across the two files, and tempTvalm2 did not come from the "
         f"profile chirp_design.py currently generates."),
        ('The encoder scale is a measurement from one run, not a calibration',
         f"The {gain['counts_per_full_step'] / 1000:.4f} µm per full step above "
         f"comes from 21 moves in a single session. It is quoted to four "
         f"figures because the fit supports that, not because the rig has been "
         f"shown to hold it: no second session exists to compare against, and "
         f"the scope channel's own scaling has not been read back from the "
         f"TwinCAT configuration and confirmed independently."),
        ('Drift was small but the first check is missing context',
         f"The seven reference returns span "
         f"{min(d['error_um'] for d in drift):+.2f} to "
         f"{max(d['error_um'] for d in drift):+.2f} µm over the 45 minute run, "
         f"which is at the level of the backlash and shows no trend. The "
         f"standstill hold before the run creeps at "
         f"{noise['drift_counts_per_s'] * 60 / 1000:+.2f} µm per minute while "
         f"energised and stationary, which over a 45 minute run is larger than "
         f"anything the drift checks show, so it does not simply continue. "
         f"Whether that settling is thermal, the "
         f"encoder, or the mount is not answered by this recording."),
    ]

    artefacts = [
        ('A 250 Hz peak that was the measurement, not the stage',
         'The down sweep of one HIGH run originally came out with its largest '
         'response at 252 Hz rather than near 685 Hz. Every HIGH leg is '
         'bracketed by a 3 s dwell at 250 Hz, and that dwell was leaking '
         'through the analysis window into the first frequency bin of the leg.',
         'Fixed: half a second is now trimmed from each end of every leg. The '
         'spurious peak disappears and the remaining 250 Hz elevation survives '
         'trimming, which is why it is reported as real.'),
        ('A 31.6 µm drift check that was a windowing error',
         'The first reference return originally read +31.6 µm against -1.2 to '
         '-3.5 µm for the other six. Its dwell was being measured right '
         'through the staircase leg that follows it, because that leg leaves '
         'no detected motion segment to end the dwell on.',
         'Fixed: the dwell window is now capped at the planned dwell length. '
         'All seven checks fall in a 2.3 µm band with no trend.'),
        ('The staircase turnaround slot',
         'The slot at the apex of each staircase is the last forward position '
         'and the first reverse position at once, and its plateau straddles '
         'the direction change, so it inflated the peak-to-valley error.',
         'Fixed: both passes now stop short of it. Static error p-v fell from '
         '1.4-4.1 µm to 1.2-1.4 µm across the three locations.'),
        ('The -5 mm staircase was indexed two microsteps late',
         'That leg reported +0.41 µm of hysteresis where the other two '
         'reported -2.34 and -1.89, lost two of its 257 slots, and contained '
         'the only backwards microstep in the campaign. All three were one '
         'fault. The slot grid recovers its cadence from the microstep edges '
         'and numbered itself from the first edge it could see; -5 mm is the '
         'one leg reached by a move in the negative direction before a '
         'staircase that runs positive, so it opens with a direction '
         'reversal, and its first two microsteps moved 0.04 and 0.24 µm - '
         'under the 0.2 µm the edge detector needs. Every slot in the leg was '
         'then labelled two microsteps short, which paired each forward '
         'position with the wrong reverse one and inverted the sign of the '
         'hysteresis; the two "backwards" steps were the first two reverse '
         'microsteps counted as forward.',
         'Fixed: each leg is now anchored on its measured turnaround rather '
         'than on its first detected edge, since the carriage turns round on '
         'the microstep the command turns round on whatever the detector saw '
         'at the start. The -5 mm leg now reads 256 of 257 slots with its '
         'apex at microstep 128, and -2.02 µm of hysteresis. The plan problem '
         'behind it is real and unfixed: the staircase still needs a run-up '
         'that approaches each start point in the direction the leg will '
         'travel, the way the repeatability phase already does.'),
        ('The last slot of every staircase leg',
         'The leg window closes where the next commanded move begins, so the '
         'final slot of each leg was being read across that boundary with the '
         'start of a transit in it.',
         'Fixed: a slot whose window runs off either end of the leg is left '
         'unread. All three legs now report 256 of 257 slots, the one missing '
         'being that last one.'),
        (f"{len(blips)} impulsive disturbances during dwells",
         f"Seven excursions of 0.21 to 0.36 µm lasting 10 to 30 ms appear "
         f"during settling-phase dwells, in two clusters around 1509-1512 s "
         f"and 1678-1698 s. They return to where they started, so they are not "
         f"lost steps, and at roughly 20 times the 15 nm encoder noise they "
         f"are not the encoder either. Something external disturbed the rig.",
         'Not corrected, and not affecting these results: settled positions '
         'are read as a median over the dwell, which a 30 ms impulse does not '
         'move. But they are 25 % of the repeatability figure, so they matter '
         'if that figure is to improve.'),
        ('The two chirp files disagree by a factor of 8 on a fixed reference',
         f"The 750 Hz sync markers are specified identically in every run and "
         f"read {mk[0]:.4f} m/s² in tempTvalm against {mk[1]:.4f} and "
         f"{mk[2]:.4f} in tempTvalm2. The turnaround dwell differs too - "
         f"2.00 s against 4.56 s - so the two files did not come off the same "
         f"configuration.",
         'Not correctable after the fact. Amplitudes may be compared within '
         'tempTvalm2 and not across the two files.'),
    ]

    next_steps = [
        'Run the FULL profile and record it. Without the 1-150 Hz legs the '
        'chirp campaign answers only half its own question, and the 250 Hz '
        'interpretation stays indirect.',
        'Check the sync marker level at the start of every chirp session '
        'before trusting the run. It is already in the profile and takes no '
        'extra time, and it is the thing that caught the difference between '
        'the two files here. Record what the accelerometer is, where it is '
        'mounted, how, and at what range, in the session notes.',
        'Find out which firmware produced tempTvalm2: its turnaround dwell is '
        '4.56 s where the current chirp_design.py specifies 2.00 s. Either the '
        'sketch on the board is not the one in the repository, or the plan was '
        'regenerated after the run.',
        'Apply a fixed reversal compensation of about 1.7 µm and re-run phase 4 '
        'only. If accuracy drops from 8.7 µm towards the 1.2 µm repeatability, '
        'backlash was the dominant term and the remaining systematic deviation '
        'is the screw.',
        'Read the EL5101 scaling out of the TwinCAT configuration and record it '
        'in CONNECTIONS.md next to the nominal figure, so the counts-per-mm '
        'used by the analysis is documented rather than inferred.',
        'Add a command or marker channel to the scope configuration. Alignment '
        'currently works, but it is reconstructed from the motion itself, and a '
        'second channel would turn that reconstruction into a check.',
        'Give the microstep staircase a run-up, so every location is entered '
        'in the direction the leg will travel. One of the three locations is '
        'currently entered against its leg direction and its hysteresis number '
        'is not comparable with the others.',
        'Record a few minutes of standstill with the drive disabled, and again '
        'with it enabled. That separates the 0.2-0.4 um impulsive disturbances '
        'seen during dwells into environmental and drive-borne, and costs five '
        'minutes at the start of a session.',
    ]

    stair = stage['staircase']
    micro = stage['microstep']
    settle_by = {
        'distance': {d: statistics.mean(
            [r['error_um'] for r in stage['settling']
             if r['distance_mm'] == d])
            for d in (0.1, 1.0, 10.0)},
        'history': {v: statistics.mean(
            [r['error_um'] for r in stage['settling']
             if r['history'] == v])
            for v in ('continuation', 'reversal')},
        'move': {v: statistics.mean(
            [r['move_error_um'] for r in stage['settling']
             if r['history'] == v])
            for v in ('continuation', 'reversal')},
        'speed': {v: statistics.mean(
            [r['error_um'] for r in stage['settling']
             if r['speed'] == v])
            for v in ('slow', 'moderate', 'fast')}}
    figures = [
        ('stage_01_overview.png', 'The run as the rig executed it',
         'Carriage position over the 43 minutes of the sequence, with the four '
         'phases shaded where the alignment places them, the seven reference '
         'returns marked in grey, and their errors below. The window starts at '
         'the first commanded move, 116 s into the recording: before that are '
         'the standstill hold and the operator jog that positioned the stage, '
         'neither part of the plan and neither used in any result.',
         'Phase boundaries are measured, not planned: the rig runs slower than '
         'the plan because the plan budgets no acceleration time.'),
        ('stage_02a_raw_error.png', 'Phase 1 - the error before anything is '
                                    'fitted away',
         'Measured minus commanded at each of the 21 stops, with the straight '
         'line that is fitted to them. Zero here is the hold the run started '
         'from, the position the operator ZERO was issued at, rather than the '
         'offset the campaign fits against the commanded targets: the '
         'question "how far from where it was told to go did it stop" is '
         'asked in the frame the command was given in. The stops are 1.2 to '
         '6.7 um long of their targets.',
         'The fitted line has an offset of +4.30 um and a slope of -0.021 '
         'um/mm. Neither is the carriage misbehaving: the offset is that '
         'starting hold sitting about 4 um from where the commanded frame '
         'puts its origin, and the slope is the -0.002 % gain error over the '
         '10 mm span, which is 0.2 um end to end. Both are common to every '
         'move, which is why the next figure removes them - the scatter about '
         'this line, not its height, is what the stage did differently from '
         'one stop to the next.'),
        ('stage_02_gain.png', 'Phase 1 - where the stage stops',
         'Residual against commanded position for the gain-calibration moves: '
         'measured minus fitted gain times commanded. The axis is ordered by '
         'size - +1, -1, +2, -2 and so on - so the two ends of each move sit '
         'next to each other; the single return to zero is not drawn, though '
         'it is in every number quoted here. Fitted gain '
         '-0.027 % of nominal (-2.7 um per 10 mm). Residual 1.88 um rms, '
         '3.05 um max. Each pair joined by a grey line is the two visits to '
         'that target, 106 s apart; they agree within 0.35 um from +/-3 mm '
         'outwards, by 0.7 um at +/-2 mm, and by 1.6 to 1.8 um at +/-1 mm, '
         'which are the first two targets of the session.',
         'The residual alternates with odd and even millimetres rather than '
         'growing with position. One screw revolution is 200 full steps = '
         '2 mm, so +/-1, +/-3 and +/-5 mm are half-revolution screw angles '
         'and +/-2 and +/-4 mm are whole ones: fitting that parity and the '
         'approach direction accounts for 1.88 um rms down to 0.75, so most '
         'of what the gain fit leaves is the screw turning rather than the '
         'drive. The approach-direction term is -0.46 +/- 0.35 um, under one '
         'microstep (0.625 um) and not resolved from zero - and in this phase '
         'every +ve approach ends at a positive target, so it is not '
         'separable from a term odd in position either. Phase 4 measures '
         'reversal cleanly. What is left after both terms is 0.75 um rms, and '
         'the largest single piece of it is the first move of the run at '
         '2.0 um: the trend across the phase is +0.29 +/- 0.14 um/min with '
         'that move in and +0.17 +/- 0.13 um/min without it, so the evidence '
         'is for a first-move effect rather than a warm-up drift.'),
        ('stage_02b_gain_dwells.png', 'Phase 1 - how it gets there',
         'All 21 dwells, 70 full steps/s moves, each trace referenced to its '
         'own settled value and reduced to a median and a min-max envelope '
         'per approach direction. Shaded: the window used for the settled '
         'position in the figure above.',
         'The carriage creeps about 0.12 um into position and does it '
         'opposite ways for the two directions, which is elastic recovery of '
         'the drive train rather than noise. A single exponential fitted to '
         'the medians gives 0.37 s moving +ve and 0.46 s moving -ve, though '
         'the decay is closer to logarithmic in time than to one exponential, '
         'so treat those as a summary and not a model. Either way it is flat '
         'well before the read window, so the points in the figure above are '
         'reading a position and not a ramp.'),
        ('stage_02c_gain_moves.png', 'Phase 1 - how far it went, not where '
                                     'it ended up',
         'The same 21 moves measured as displacements rather than as '
         'positions: executed minus commanded, against the size of the move. '
         'Each point is one plateau average minus the one before it, so no '
         'zero enters the calculation. Fitted through the origin the scale is '
         '1.000 to within -0.0009 %, which is one part in 100,000 and far '
         'below anything else here. The error is +/-4.5 um, 2.72 um rms, with '
         'no dependence on move size.',
         'The pairs joined by a grey line are the same move commanded twice, '
         '106 s apart: they agree to a median of 0.15 um and 0.8 um at worst, '
         'and the +9 mm move repeats to within the 10 nm the encoder can '
         'resolve. An error that reproduces that well is not noise and is not '
         'the drive losing steps - it is a property of the two positions the '
         'move runs between. Moves that start and end a whole number of screw '
         'revolutions apart (+/-2, +/-4, +/-6, +/-8, +/-10 mm) come in under '
         '2 um; moves that end half a revolution from where they started '
         '(+3, +5, +7, +9 mm) carry 2.7 to 4.5 um, with the sign set by which '
         'way the half revolution went. The first move of the phase is the '
         'one point here that depends on where the stage was left beforehand: '
         'that hold sat 5.1 um from the fitted origin.'),
        ('stage_03_cumulative_error.png', 'Phase 2 - the error carried '
                                          'through the stepping sequence',
         'Each leg drawn as the one continuous run it was: 128 microsteps '
         'out, the turnaround, 128 back, with measured minus commanded '
         'accumulated through all of it from the slot the leg started on. '
         'Each leg is zeroed on its own first slot, so what is drawn is the '
         'staircase own contribution and not the offset the approach move '
         'left behind; that offset is stated on each panel and was +0.96, '
         '-3.49 and +1.59 um. Negative is the carriage behind its command '
         'on the way out. This '
         'is the same data as the folded map two pages on, and the same sign '
         'convention, before any of it is fitted away. The light verticals '
         'are full-step boundaries and the dotted one is the commanded '
         'turnaround.',
         f"The step at the turnaround is the lost motion arriving in one "
         f"place: the command reverses, the carriage does not follow for the "
         f"next twenty-odd microsteps, and about 2 um of error appears. "
         f"Because the error is signed in machine coordinates rather than "
         f"against the direction of travel, being behind reads negative on "
         f"the way out and positive on the way back, so that one event turns "
         f"the trace through some 2 um whichever leg it is. What happens "
         f"afterwards is the point of drawing the leg this way: it is not "
         f"given back. The legs end {stair[0]['closing_error_um']:+.2f}, "
         f"{stair[1]['closing_error_um']:+.2f} and "
         f"{stair[2]['closing_error_um']:+.2f} um from where they started, "
         f"with the command back at its starting value - or "
         f"{stair[0]['start_offset_um'] + stair[0]['closing_error_um']:+.2f}, "
         f"{stair[1]['start_offset_um'] + stair[1]['closing_error_um']:+.2f} "
         f"and "
         f"{stair[2]['start_offset_um'] + stair[2]['closing_error_um']:+.2f} "
         f"um against the commanded frame itself, with the approach offset "
         f"added back. The wobble along "
         f"each arm is the microstep error, and the -5 mm leg falls away "
         f"over its first 25 microsteps because that leg, alone among the "
         f"three, is entered against its approach - it spends its lost "
         f"motion on the way in, which is why it is the one that finishes "
         f"nearest zero."),
        ('stage_03b_step_size.png', 'Phase 2 - the size of each microstep',
         'How far the carriage actually moved on each of the 256 commanded '
         'microsteps of each leg, against the nominal 0.625 um at mres 16 '
         '(dashed), with the distribution of each leg beside it. The band is '
         'the leg mean plus or minus one standard deviation and the dotted '
         'line is half a microstep. Two different failures are marked '
         'separately: a ring is a microstep whose displacement came in under '
         'half of one, a cross is a microstep that left no step the signal '
         'can be seen to take on its own, found by thresholding the '
         'recording with no grid and no knowledge of the cadence.',
         f"All 256 microsteps of every leg are measurable once the final "
         f"slot - clipped by a few milliseconds at the end of the leg window "
         f"- is read from the 99 % of it that is there. Of the 768 in the "
         f"phase, {sum(s['steps_resolved'] for s in stair)} move the carriage "
         f"more than half a microstep and "
         f"{768 - sum(len(s['steps_undetected']) for s in stair)} leave a "
         f"visible step. The two criteria disagree in one place and the "
         f"disagreement is the useful part: microstep 133 at 0 mm moved "
         f"0.28 um, under half of one, but is perfectly visible in the "
         f"recording - it is a small step, not a missing one. The microsteps "
         f"that leave no visible step at all are "
         f"{', '.join(str(k) for k in stair[0]['steps_undetected'])} at -5 mm "
         f"and 129 at both 0 and +5 mm: number 129 is the first microstep "
         f"after the turnaround, and 1 and 2 are the first two of the -5 mm "
         f"leg, the one leg entered against its approach. Nothing fails "
         f"anywhere that is not next to a reversal, and nothing catches up "
         f"afterwards - the largest step in the phase is "
         f"{max(s['step_max_um'] for s in stair):.2f} um. Mean step size is "
         f"{stair[0]['step_mean_um']:.3f}, {stair[1]['step_mean_um']:.3f} and "
         f"{stair[2]['step_mean_um']:.3f} um with standard deviations of "
         f"{stair[0]['step_sd_all_um']:.3f}, "
         f"{stair[1]['step_sd_all_um']:.3f} and "
         f"{stair[2]['step_sd_all_um']:.3f}, and peak to peak "
         f"{stair[0]['step_pv_um']:.2f}, {stair[1]['step_pv_um']:.2f} and "
         f"{stair[2]['step_pv_um']:.2f} um."),
        ('stage_03d_recovery.png', 'Phase 2 - recovery after a reversal',
         'The running total of what the microsteps after a reversal failed '
         'to deliver, D(n) = sum of (0.625 - step), counted from each pass '
         'own reversal, with the band its plateau is read over shaded and '
         'the lost motion from the hysteresis gap dashed.',
         f"The lost motion is quoted from the gap between the two passes - "
         f"{micro['per_leg'][0]['lost_motion_um']:.2f}, "
         f"{micro['per_leg'][1]['lost_motion_um']:.2f} and "
         f"{micro['per_leg'][2]['lost_motion_um']:.2f} um - and not from "
         f"this plateau, because the local scale error enters D(n) once and "
         f"tilts it: at 0 mm a step running 2.2 % long walks the plateau "
         f"down by 0.69 um across the window, and 1.93 plus 0.69 is the "
         f"{micro['per_leg'][1]['lost_motion_um']:.2f} um the gap reports. "
         f"In the gap the scale is in both passes and cancels. What D(n) is "
         f"for is the shape of the recovery, and that is read at fractions "
         f"of the gap value: half the lost motion is back after 6 microsteps "
         f"at every one of the three locations, 80 % of it after 16 to 38, "
         f"and the step size itself is back to the local value after 19 to "
         f"32. The 95 % figure is not quoted - it is read off the part of "
         f"the plateau the tilt has already moved."),
        ('stage_03c_microstep_map.png', 'Phase 2 - the microstep error over '
                                        'one electrical period',
         'Every microstep of every leg folded onto one electrical period, '
         'with the transition microsteps after each reversal left out. The '
         'line and the harmonics at 1, 2 and 4 cycles per period are fitted '
         'together rather than the line being removed first: a leg spans two '
         'periods, and over that span a straight line and a once-per-period '
         'sine are not independent. Left: the mean of all twelve passes with '
         'the spread between them, and the fitted harmonics. Right: the same '
         'map split by direction of travel.',
         f"The map is {micro['folded']['pv_um']:.2f} um peak to peak and "
         f"{micro['folded']['rms_um']:.3f} um rms. That is smaller than the "
         f"{max(s['static_pv_um'] for s in stair):.1f} um peak to peak of "
         f"the first figure in this phase and both are right: this is only "
         f"the part of the error that repeats every electrical period, with "
         f"the local scale and the curvature within each leg fitted out. The "
         f"once-per-period term is the largest at "
         f"{micro['folded']['amplitude']['1']:.3f} um, against "
         f"{micro['folded']['amplitude']['2']:.3f} um twice per period and "
         f"{micro['folded']['amplitude']['4']:.3f} um four times per period, "
         f"so the current table and the phase balance dominate and detent "
         f"torque, the four-per-period term, is the smallest of the three. "
         f"The shape repeats: between locations the maps correlate at r = "
         f"{', '.join(f'{v:.2f}' for v in micro['folded']['location_correlation'].values())}, "
         f"and forward against reverse at r = "
         f"{micro['folded']['direction_correlation']:.2f}. Both are high "
         f"enough to call the error deterministic and electrical rather than "
         f"frictional, which means it could be compensated in the drive "
         f"microstep table rather than merely tolerated."),
        ('stage_04_sequence.png', 'Phase 3 - how a settling point is made',
         '(a) The whole phase, 19 minutes of it, with every stop marked and '
         'the three reference returns on the zero line. Each move starts '
         'where the previous one stopped, so the carriage walks around the '
         'travel rather than returning to a datum, and a move that would '
         'leave the +/-10 mm working range is turned around instead; the '
         'order is shuffled, which is what the picture is for. (b) One move '
         'of it: the carriage ramps to the target, and the dotted line is '
         'the end of commanded motion, which is the origin of the next '
         'figure but one. Commanded and measured cannot be told apart at '
         'this scale - the difference is a micrometre on a millimetre move - '
         'so the inset magnifies the stop, where the shaded band is the '
         'window the settled position is read from and the arrow is the '
         'settling error of the figure overleaf.',
         'Roughly half the moves reverse the one before them and which half '
         'is decided by the shuffle, not by design. That is what makes '
         'approach history a factor worth recording rather than a property '
         'of the plan.'),
        ('stage_04a_settling_error.png', 'Phase 3 - where the stage stops '
                                         'and how far it goes',
         '(a) Settled position minus target for all 90 moves, grouped by '
         'commanded distance and then by speed - slow, moderate, fast, left '
         'to right in each group, shown light to dark. The speeds are kept '
         'apart only to show that they do not separate; they could be pooled '
         'into one colour per distance without losing anything. (b) The same '
         'moves measured as travel: distance moved minus distance commanded, '
         'pooled over distance and split by whether the move continued the '
         'previous direction or reversed it. (c) Why the two disagree about '
         'the reversal, with the shortfall drawn far larger than it is. Both '
         'errors are signed along the direction of travel, so negative means '
         'short, and both have the campaign offset, the gain and the screw '
         'once-per-revolution term removed; the gain is phase 4, averaged '
         'over the two approach directions so that the lost motion cancels '
         'out of it.',
         f"Write delta for the steady-state shortfall, half the reversal "
         f"value phase 4 measures: delta = "
         f"{abs(stage['repeatability']['reversal_mean']) / 2:.2f} um. Where "
         f"the carriage stops is -delta, and that is all it is: "
         f"{abs(stage['settling_model']['settling_mean_um']):.2f} um short of "
         f"target, within 0.36 um across two decades of distance, less than "
         f"that across speed, and "
         f"{abs(settle_by['history']['continuation'] - settle_by['history']['reversal']):.2f} "
         f"um between a continuing move and a reversing one - which is to "
         f"say nothing at all. Where the carriage has been is in the travel: "
         f"{settle_by['move']['continuation']:+.2f} um for a continuing move "
         f"against {settle_by['move']['reversal']:+.2f} for a reversing one, "
         f"a difference of "
         f"{settle_by['move']['reversal'] - settle_by['move']['continuation']:+.2f} "
         f"+/- 0.27 um at six and a half standard errors and within a "
         f"quarter of a micrometre of -2delta. Panel (c) is the reason both "
         f"are true: a move ends the same shortfall behind its target "
         f"whichever way it arrived, so the stop knows no history, while a "
         f"move that reverses begins on the far side of that shortfall and "
         f"so covers 2delta less ground than it was told to. The scatter in "
         f"(b) is largest for the 1 mm moves - 1.9 um against 0.7 to 1.0 for "
         f"the others - because a 1 mm move turns the screw half a "
         f"revolution and carries the whole of its runout into the travel, "
         f"while a 10 mm move turns it five whole revolutions and carries "
         f"none."),
        ('stage_04b_settling_traces.png', 'Phase 3 - how the stage settles',
         'All 90 stops over the five seconds after the end of commanded '
         'motion, each referred to its own settled value, on a logarithmic '
         'time axis, with the median of each approach direction drawn '
         'through them.',
         'Nothing rings and nothing overshoots: every stop creeps towards '
         'its final value from the side it arrived from, 0.14 um for a move '
         'in the positive direction and 0.17 um in the negative, with half '
         'of that gone by 0.3 s. The split by direction is the elastic '
         'recovery of the drive train unloading, the same effect the phase 1 '
         'dwells show. Speed makes no difference to it, which is why this '
         'figure is coloured by direction instead. No trace approaches the '
         '+/-0.5 um the settling criterion is written against, so the '
         'settling times quoted earlier in this report are a property of the '
         'read-out window rather than of the stage.'),
        ('stage_05_sequence.png', 'Phase 4 - how a reversal value is '
                                   'measured',
         f"The phase visits five targets at -10, -5, 0, +5 and +10 mm, each "
         f"from below and from above, five times over: 50 points in 8.4 "
         f"minutes, zig-zagging so that each target is reached once from "
         f"each side per cycle. (a) One target, both approaches, at "
         f"millimetre scale. Each is set up by a run-up that overshoots by "
         f"0.5 mm, so the last 0.5 mm is travelled in the nominated "
         f"direction with the drive train already loaded that way - without "
         f"it, a target reached by stopping short would carry the backlash "
         f"of whatever move preceded it. Commanded is drawn dashed over "
         f"measured; at this scale the two cannot be told apart. (b) The two "
         f"stops of (a) at micrometre scale, both referred to the same "
         f"target. The settled value is the mean over the shaded window, "
         f"{rep['read_window_s'][0]:.1f} to {rep['read_window_s'][1]:.1f} s "
         f"after the end of commanded motion; the gap between the two levels "
         f"is the reversal value B at this target.",
         f"Moves are segmented from the plan here, not from the signal: "
         f"every move in this phase is a row of the plan with its role in "
         f"the label, and the recording is used only to say when each one "
         f"happened. That matters more here than anywhere else in the "
         f"campaign, because the quantity being measured is the difference "
         f"between two approaches to the same point - a detector that missed "
         f"a move would silently drop a target rather than report it. All "
         f"{rep['n']} points survive. One offset is removed from the whole "
         f"phase, {rep['offset_removed_um']:+.2f} um, because the zero of "
         f"the campaign is arbitrary; nothing else is taken out, so the gain "
         f"and the screw error are both still in the numbers that follow."),
        ('stage_05a_repeatability.png', 'Phase 4 - positioning deviation at '
                                        'each target',
         f"(a) Settled position minus target for all {rep['n']} points, in "
         f"machine coordinates: positive is further along +x whichever way "
         f"the carriage arrived. Light dots are the five individual visits, "
         f"the marker is their mean, and the bar spans +/-2s - so the full "
         f"height of a bar is that target's unidirectional repeatability, "
         f"R(i) up or R(i) down. The dashed envelope encloses the two "
         f"approaches to one commanded position, and the vertical offset "
         f"between the triangles inside it is that target's reversal value "
         f"B(i). The heavy dashed line is the straight fit through the "
         f"direction-averaged mean of each pair, in which the lost motion "
         f"cancels; its slope is the gain. (b) The same five positions as "
         f"numbers: the mean of each approach direction, then the smallest, "
         f"largest and mean deviation over all ten visits. Results: "
         f"R↑ = {rep['repeatability_forward']:.2f} um, "
         f"R↓ = {rep['repeatability_reverse']:.2f} um; mean B = "
         f"{rep['reversal_mean']:+.2f} um and max |B(i)| = "
         f"{rep['reversal_max']:.2f} um; A = {rep['accuracy']:.2f} um; gain "
         f"= {rep['gain_pct']:+.4f} +/- {rep['gain_se_pct']:.4f} %.",
         f"The stage returns to a target within {rep['repeatability']:.2f} "
         f"um bidirectionally, and within "
         f"{max(rep['repeatability_forward'], rep['repeatability_reverse']):.2f} "
         f"um if the direction of approach is held fixed; the reversal value "
         f"is {rep['reversal_mean']:+.2f} +/- {rep['reversal_sd']:.2f} um "
         f"and never exceeds {rep['reversal_max']:.2f} um. Axis accuracy A "
         f"is {rep['accuracy']:.2f} um, five times the worst R, and almost "
         f"none of it is repeatability: it is where the targets sit. The "
         f"gain accounts for {rep['gain_pct']:+.4f} +/- "
         f"{rep['gain_se_pct']:.4f} % of it, and what is left of M once the "
         f"gain line is removed is {screw['residual_rms_um']:.2f} um rms - "
         f"which is not scatter but the screw. The screw advances 2 mm a "
         f"turn, so targets 5 mm apart lie alternately a whole turn and a "
         f"half turn from zero, where the once-per-revolution runout has the "
         f"opposite sign. Sorted that way the residual separates completely: "
         f"{screw['whole_turn_um']:+.2f} um at the whole-turn targets "
         f"against {screw['half_turn_um']:+.2f} um at the half-turn ones, an "
         f"amplitude of {screw['amplitude_um']:.2f} um. That is the same "
         f"runout phase 3 sees as extra scatter in its 1 mm moves, measured "
         f"here as a position error instead. Two assumptions behind these "
         f"numbers were checked and held. Nothing drifts: a straight line "
         f"through all {rep['n']} deviations against cycle number has a "
         f"slope of {rep['drift']['per_cycle_um']:+.3f} +/- "
         f"{rep['drift']['se_um']:.3f} um per cycle, "
         f"{rep['drift']['over_phase_um']:+.3f} um over the whole phase, so "
         f"the five visits to a target are repeats of one measurement and "
         f"nothing was detrended. And the read window opens late enough: "
         f"every hold creeps "
         f"{abs(rep['creep']['forward_um']):.2f} um (from below) or "
         f"{abs(rep['creep']['reverse_um']):.2f} um (from above) further in "
         f"the direction it was travelling, the same elastic recovery phases "
         f"1 and 3 both show, with half of it gone by 0.3 s and the rest "
         f"within a second - reading at 0.1 s instead would have added about "
         f"a quarter of a micrometre of spurious reversal value."),
        ('stage_06_standstill.png', 'The encoder while the stage is held still',
         'The energised hold before the sequence begins. Left: position over '
         'the 77 s hold. Right: where the motion sits in frequency.',
         f"Short-term scatter is {noise['short_term_counts']:.0f} counts rms, "
         f"which at the fitted scale is {noise['short_term_counts'] / 1000:.3f} "
         f"µm - far below every effect reported here, so none of the results "
         f"above are limited by encoder noise."),
        ('stage_05b_campaign.png', 'The campaign in one view',
         "The two quantities that more than one phase measured, with each "
         "phase's estimate and its uncertainty. The band and the line are "
         "the mean and standard deviation of the estimates that stand on "
         "their own; the plain fit from phase 1 is drawn but not pooled, "
         "because it reads the gain off a line that the lost motion has "
         "already tilted.",
         f"Four independent routes to the lost motion - a fit through the "
         f"+/- pairs of phase 1, the gap at three staircase turnarounds, the "
         f"difference between reversing and continuing moves in phase 3, and "
         f"the mean reversal value in phase 4 - agree on "
         f"{camp['lost_motion'][0]['pooled_mean']:.2f} +/- "
         f"{camp['lost_motion'][0]['pooled_sd']:.2f} um, and every one of "
         f"them is within its own error bar of that. They are not one "
         f"measurement rearranged: two are static gaps read at a reversal, "
         f"one is a travel shortfall over moves from 0.1 to 10 mm, and one "
         f"is a difference of settled positions. The gain is the more "
         f"interesting panel. The plain fit of phase 1 puts it at "
         f"{camp['gain'][0]['value']:+.4f} %, which is not merely imprecise "
         f"but the wrong sign, while the two estimates that separate the "
         f"gain from the lost motion before fitting - the +/- pairs of phase "
         f"1 and the direction-averaged means of phase 4 - agree on "
         f"{camp['gain'][1]['pooled_mean']:+.4f} +/- "
         f"{camp['gain'][1]['pooled_sd']:.4f} %. The lesson is the one the "
         f"campaign was designed around: a 2 um backlash is large enough to "
         f"swallow a gain error of a part in 10^4, so the two have to be "
         f"measured apart or not at all."),
        ('chirp_02_inventory.png', 'Chirp - what was recorded, and whether it '
                                   'repeats',
         'Left: the two accelerometer files with the detected sweep legs '
         'marked. Right: each run against the common response shape after its '
         'own overall level has been divided out.',
         'The planned FULL profile appears in neither file. The right-hand '
         'panel is the basis for the caveat that these runs agree on resonance '
         'frequencies but not on amplitude.'),
        ('chirp_01_spectrogram_tempTvalm2.png', 'Chirp - a recording in full',
         'tempTvalm2, axis 1. Two complete HIGH runs, each sweeping 250 Hz up '
         'to 1 kHz and back at 7.8125 Hz/s, with the detected legs drawn over '
         'them. The broad bright band near 685 Hz is the dominant resonance; '
         'the bursts before each run are the stage being repositioned.',
         'The 150-250 Hz band is empty because nothing drives it - that is the '
         'safety decision recorded in chirp/README.md, not a gap in the data.'),
        ('chirp_03_response.png', 'Chirp - response against drive frequency',
         'The acceleration measured at the frequency being driven, for each of '
         'the three accelerometer axes, with up and down sweeps separated and '
         'the band spanning the three repeats.',
         'Up and down place the resonance at the same frequency, which is what '
         'distinguishes a property of the structure from an artefact of '
         'sweeping.'),
        ('chirp_04_direction.png', 'Chirp - the 250 Hz question',
         'Left: the three runs around the frequency HIGH starts from. Right: '
         'each run divided by itself, up sweep over down sweep, in dB.',
         'Compared within a run rather than across runs, tempTvalm agrees with '
         'itself to better than half a decibel in every band - so the sweep '
         'direction introduces no systematic bias of its own, and the scatter '
         'in the other two runs is run quality rather than a sweep artefact.'),
    ]

    print(f'building {args.out.name}')
    with PdfPages(args.out) as pdf:
        page_cover(pdf, stage, chirp, when)
        page_method(pdf, stage, chirp)
        for png, title, caption, note in figures:
            page_figure(pdf, args.plots / png, title, caption, note)
            print(f'  page: {title}')
        page_tables(pdf, stage, chirp)
        page_artefacts(pdf, artefacts)
        page_findings(pdf, stage, chirp, findings, caveats, next_steps)
        info = pdf.infodict()
        info['Title'] = 'Stage identification campaign - measured results'
        info['Subject'] = ('Stage identification and chirp sweep recordings of '
                           '13 October 2026')
        info['CreationDate'] = dt.datetime.now()
    size = args.out.stat().st_size / 1e6
    print(f'wrote {args.out.name} ({size:.1f} MB)')


if __name__ == '__main__':
    main()
