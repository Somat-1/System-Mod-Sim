#!/usr/bin/env python3
'''Process the stage identification recording into results and figures.

The recording is a single EL5101 counter channel at 1 kHz with no command
channel in it, so the first job is to find where in the recording each planned
move actually happened. That is done from the data rather than from the clock:

1. threshold a smoothed |dcounts| to find every interval the stage was moving;
2. align those intervals to the planned move list by dynamic programming over
   commanded displacement, which tolerates planned moves that leave no
   detectable segment (the quasi-static staircase legs) and detected segments
   that belong to no planned move (the operator jog before the run);
3. fit the counts-per-full-step scale on the gain phase, which is the phase
   that exists to measure it;
4. read the settled position out of each dwell and derive the per-phase results.

Written for the rig: `python stage_analysis.py` after `cache_recordings.py`.
'''

import argparse
import csv
import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
import numpy as np

import vizstyle as vs

HERE = Path(__file__).resolve().parent
STAGE_DIR = HERE.parent
FS = 1000.0                       # EL5101 scope sample rate, Hz
FULL_STEP_UM = 10.0               # 100 full steps per mm
NOMINAL_COUNTS_PER_REV = 20_000   # the value CONNECTIONS.md calls provisional
NOMINAL_COUNTS_PER_FULL_STEP = 10_000  # 1 count = 1 nm, 10 um per full step
FULL_STEPS_PER_REV = 200

MOVE_THRESHOLD = 100              # counts per sample, on a 51-sample mean
MOVE_MERGE_MS = 60
MICROSTEP_EDGE = 200              # counts in one sample marks a commanded step
SETTLE_TOL_UM = (2.0, 1.0, 0.5, 0.2)   # settling time is reported at each


# --------------------------------------------------------------------------
# loading and alignment
# --------------------------------------------------------------------------

def load_plan(path):
    rows = list(csv.DictReader(path.open(encoding='utf-8')))
    for r in rows:
        for k in ('t_start', 't_end', 'dwell_s', 'micro_dwell_s'):
            r[k] = float(r[k])
        for k in ('start', 'target', 'speed', 'micro_steps', 'mres'):
            r[k] = int(r[k])
    return rows


def motion_segments(counts):
    '''Intervals during which the stage was moving, as (start, end) indices.'''
    speed = np.convolve(np.abs(np.diff(counts)), np.ones(51) / 51, 'same')
    moving = speed > MOVE_THRESHOLD
    edges = np.diff(moving.astype(np.int8))
    starts = np.flatnonzero(edges == 1) + 1
    ends = np.flatnonzero(edges == -1) + 1
    if moving[0]:
        starts = np.r_[0, starts]
    if moving[-1]:
        ends = np.r_[ends, moving.size]
    keep_s, keep_e = [starts[0]], []
    for i in range(1, starts.size):
        if starts[i] - ends[i - 1] < MOVE_MERGE_MS:
            continue
        keep_e.append(ends[i - 1])
        keep_s.append(starts[i])
    keep_e.append(ends[-1])
    return np.array(keep_s), np.array(keep_e)


def align(plan_delta, obs_delta, sigma=1.0, gap=-3.0, match_bonus=6.0):
    '''Needleman-Wunsch over commanded displacement.

    Both sides may skip: a staircase leg is 128 microsteps too small to trip
    the motion threshold and so has no segment, and the jog that positioned the
    stage before the run has no planned move.
    '''
    n, m = plan_delta.size, obs_delta.size
    score = np.full((n + 1, m + 1), -np.inf)
    back = np.zeros((n + 1, m + 1), np.int8)
    score[0, 0] = 0.0
    for i in range(n + 1):
        for j in range(m + 1):
            if i == 0 and j == 0:
                continue
            best, arg = -np.inf, 0
            if i and j:
                resid = (obs_delta[j - 1] - plan_delta[i - 1]) / sigma
                cand = score[i - 1, j - 1] + match_bonus - min(resid ** 2, 20.0)
                if cand > best:
                    best, arg = cand, 1
            if i and score[i - 1, j] + gap > best:
                best, arg = score[i - 1, j] + gap, 2
            if j and score[i, j - 1] + gap > best:
                best, arg = score[i, j - 1] + gap, 3
            score[i, j], back[i, j] = best, arg

    i, j = n, m
    pairs, skipped, extra = [], [], []
    while i or j:
        move = back[i, j]
        if move == 1:
            pairs.append((i - 1, j - 1))
            i, j = i - 1, j - 1
        elif move == 2:
            skipped.append(i - 1)
            i -= 1
        else:
            extra.append(j - 1)
            j -= 1
    return pairs[::-1], skipped[::-1], extra[::-1]


def fit_scale(counts, seg_s, seg_e, plan, shift, n=20):
    '''Counts per full step, from the gain phase, through the origin.'''
    obs = np.array([counts[seg_e[i + shift]] - counts[seg_s[i + shift]]
                    for i in range(n)], float)
    cmd = np.array([plan[i]['target'] - plan[i]['start'] for i in range(n)], float)
    gain = float(obs @ cmd / (cmd @ cmd))
    resid = obs / gain - cmd
    return gain, float(np.sqrt(np.mean(resid ** 2))), float(np.abs(resid).max())


def build_alignment(counts, plan):
    seg_s, seg_e = motion_segments(counts)
    # A coarse scale is enough to align on; it is refit properly afterwards.
    coarse = (counts.max() - counts.min()) / (
        max(r['target'] for r in plan) - min(r['start'] for r in plan))
    best = None
    for shift in (0, 1, 2):
        gain, rms, _ = fit_scale(counts, seg_s, seg_e, plan, shift)
        if best is None or rms < best[1]:
            best = (shift, rms, gain)
    shift, _, gain = best
    obs_delta = (counts[seg_e] - counts[seg_s]) / gain
    plan_delta = np.array([r['target'] - r['start'] for r in plan], float)
    pairs, skipped, extra = align(plan_delta, obs_delta)
    del coarse
    return dict(seg_s=seg_s, seg_e=seg_e, pairs=pairs, skipped=skipped,
                extra=extra, gain=gain, shift=shift)


# --------------------------------------------------------------------------
# per-move extraction
# --------------------------------------------------------------------------

def settled_positions(q, align_info, plan):
    '''For every matched move: the settle window and the position it held.

    The window runs from the end of the move to the start of the next detected
    motion, which is the honest definition of the dwell: it does not assume the
    planned dwell length, and the plan's clock has drifted from the rig's by
    this point in the run anyway.
    '''
    seg_s, seg_e = align_info['seg_s'], align_info['seg_e']
    out = {}
    for i, j in align_info['pairs']:
        stop = seg_e[j]
        nxt = seg_s[j + 1] if j + 1 < seg_s.size else q.size
        # Cap the window at the dwell the plan asked for. Without the cap, a
        # move followed by a staircase leg - which leaves no detected segment -
        # would have its "dwell" run on through the whole leg, and the position
        # read out of it would be taken from the middle of a moving ramp.
        nxt = min(nxt, stop + int((plan[i]['dwell_s'] + 2.0) * FS))
        if nxt - stop < 100:                      # shorter than 100 ms: no dwell
            continue
        window = q[stop:nxt]
        tail = window[int(0.6 * window.size):]
        out[i] = dict(
            move_start=int(seg_s[j]), move_end=int(stop), dwell_end=int(nxt),
            settled=float(np.median(tail)), noise=float(tail.std()),
            dwell_s=(nxt - stop) / FS,
            move_s=(stop - seg_s[j]) / FS,
            target=plan[i]['target'])
    return out


def settling_times(q, record):
    '''Time from the end of commanded motion until it stays inside a band.'''
    window = q[record['move_end']:record['dwell_end']]
    final = record['settled']
    times = {}
    for tol_um in SETTLE_TOL_UM:
        tol = tol_um / FULL_STEP_UM
        outside = np.flatnonzero(np.abs(window - final) > tol)
        times[tol_um] = (outside[-1] + 1) / FS if outside.size else 0.0
    return times


def detect_steps(y_um, jump_um=0.3, smooth=41, half=100, refractory=300):
    '''Step times found in the signal alone, with no grid and no cadence.

    A boxcar well under one slot, then the change across a fixed window
    either side of each sample. A microstep moves most of 0.625 um in a few
    milliseconds, so it stands clear of the 0.05 um the encoder carries at
    standstill. The run midpoint is taken rather than the largest sample:
    the difference is flat across the step, so its maximum is picked by
    noise, while the midpoint is the step itself.

    This answers a different question from the size of a step. A microstep
    can be detected here and still have moved the carriage less than half of
    one, and the two are reported separately because they fail separately:
    the first microstep after a reversal usually leaves no detectable step at
    all, while a few others are detectable but small.
    '''
    kernel = np.ones(smooth) / smooth
    sm = np.convolve(y_um, kernel, 'same')
    change = np.zeros_like(sm)
    change[half:-half] = sm[2 * half:] - sm[:-2 * half]

    times = []
    over = np.flatnonzero(np.abs(change) > jump_um)
    if over.size:
        for run in np.split(over, np.flatnonzero(np.diff(over) > 1) + 1):
            if run.size < 5:
                continue
            peak = int((run[0] + run[-1]) // 2)
            if times and peak - times[-1] < refractory:
                continue
            times.append(peak)
    return np.array(times, dtype=int)


def _read_slots(q, a, b, start, period, n_slots):
    '''Read every slot of a staircase grid: the median of the last half of
    the dwell that ends at that slot's edge.

    A slot is read from the last 40 % of its own dwell. The final slot of a
    leg is allowed to be clipped by the end of the leg window, since the
    window closes where the next commanded move begins and the last dwell
    runs right up to it - but only slightly: a slot missing more than a tenth
    of itself is left unread rather than read from whatever is there.
    '''
    pos = np.full(n_slots + 1, np.nan)
    noise = np.full(n_slots + 1, np.nan)
    for k in range(n_slots + 1):
        lo = int(round(start + (k - 1) * period)) if k else 0
        hi = int(round(start + k * period))
        if lo < 0 or hi - lo < 100:
            continue
        if hi > b - a:
            if hi - (b - a) > 0.1 * period:
                continue
            hi = b - a
        window = q[a + lo:a + hi][(hi - lo) // 2:]
        pos[k] = float(np.median(window))
        # The scatter the median was taken over, as a standard error. The
        # samples in a window are not independent - the carriage creeps
        # through it - so this is a floor on the read-out noise, not an
        # estimate of it.
        noise[k] = float(window.std() / np.sqrt(window.size))
    return pos, noise


def staircase_legs(q, counts, align_info, plan):
    '''Microstep-level position for each quasi-static staircase window.

    One window holds the forward leg and the reverse leg that retraces it, and
    the pair has to be read together: the hysteresis is the difference between
    the two passes over the same commanded positions.

    The microsteps are not found by thresholding each one. Some of them move
    the stage too little to clear the noise, and those are exactly the ones
    worth keeping, so dropping them would bias the result towards the steps
    that worked. Instead the commanded cadence is recovered - it is a fixed
    dwell per microstep - and every slot on that grid is read out, whether the
    stage moved in it or not.
    '''
    seg_s, seg_e = align_info['seg_s'], align_info['seg_e']
    row_of = {i: j for i, j in align_info['pairs']}
    legs, seen = [], set()
    for i in sorted(align_info['skipped']):
        if plan[i]['micro_steps'] == 0:
            continue
        before = max((k for k in row_of if k < i), default=None)
        after = min((k for k in row_of if k > i), default=None)
        if before is None or after is None or (before, after) in seen:
            continue
        seen.add((before, after))
        a, b = int(seg_e[row_of[before]]), int(seg_s[row_of[after]])

        edge = np.flatnonzero(np.abs(np.diff(counts[a:b])) > MICROSTEP_EDGE) + 1
        if edge.size < 100:
            continue
        edge = edge[np.r_[True, np.diff(edge) > 300]]      # 300 ms refractory
        period = float(np.median(np.diff(edge)))
        index = np.rint((edge - edge[0]) / period)
        period, start = np.polyfit(index, edge.astype(float), 1)

        # The cadence is recovered reliably; its phase is not. The grid above
        # is numbered from the first microstep large enough to leave an edge,
        # and a leg entered against its approach spends its first microsteps
        # taking up lost motion - they move too little to leave one. That is
        # not a lost sample or two: it renumbers every slot in the leg, pairs
        # each forward position with the wrong reverse one, and inverts the
        # hysteresis. The apex is the fixed point that recovers the phase.
        # The carriage turns round on the microstep the command turns round
        # on, so the measured turnaround is commanded microstep n, whatever
        # the detector saw at the start of the leg.
        n_up = plan[i]['micro_steps']
        expected = 2 * n_up
        provisional, _ = _read_slots(q, a, b, start, period,
                                     int(round((b - a - start) / period)))
        shift = 0
        if np.any(np.isfinite(provisional)):
            shift = int(np.nanargmax(provisional)) - n_up
            start += shift * period
        pos, noise = _read_slots(q, a, b, start, period, expected)

        # Which commanded microsteps left a step the signal can be seen to
        # take on its own. Grid edge k divides slot k from slot k + 1, so the
        # step across it is the one arriving at k + 1.
        found = {int(round((t - start) / period)) + 1
                 for t in detect_steps((q[a:b] - q[a]) * FULL_STEP_UM)}
        undetected = [k for k in range(1, expected + 1) if k not in found]

        legs.append(dict(window=(a, b), position=pos,
                         period_s=float(period / FS), slots=expected,
                         # Where the recovered grid starts, in samples from
                         # the window, so a slot can be put back on the clock.
                         grid_start=float(start), grid_shift=shift,
                         read=int(np.sum(np.isfinite(pos))),
                         slot_noise=noise, undetected=undetected,
                         expected_slots=expected,
                         micro_steps=n_up,
                         mres=plan[i]['mres'], location=plan[i]['start']))
    return legs


# --------------------------------------------------------------------------
# phase results
# --------------------------------------------------------------------------

def line_fit(x, y):
    '''Straight line with the standard error of its slope.

    A slope quoted without one cannot answer "is this a trend or is this
    scatter", which is the only question the time-ordered figure asks.
    '''
    slope, intercept = np.polyfit(x, y, 1)
    resid = y - (slope * x + intercept)
    s2 = (resid @ resid) / (x.size - 2)
    sxx = ((x - x.mean()) ** 2).sum()
    se = float(np.sqrt(s2 / sxx))
    se_intercept = float(np.sqrt(s2 * (1 / x.size + x.mean() ** 2 / sxx)))
    return dict(slope=float(slope), intercept=float(intercept), se=se,
                se_intercept=se_intercept)


def gain_result(records, plan, gain):
    idx = sorted(i for i in records if plan[i]['phase'] == 'gain')
    cmd = np.array([plan[i]['target'] for i in idx], float)
    obs = np.array([records[i]['settled'] for i in idx], float)
    slope, intercept = np.polyfit(cmd, obs, 1)
    resid = obs - (slope * cmd + intercept)
    counts_per_fs = gain * slope

    # Which way the carriage was moving into each target, and whether this is
    # the first or the second time that target was approached that way: the
    # phase visits every target twice, so the pair is a repeatability check
    # and the difference between the two directions carries the lost motion.
    direction, visit, seen = [], [], {}
    for i in idx:
        d = 1 if plan[i]['target'] > plan[i]['start'] else -1
        seen[(plan[i]['target'], d)] = seen.get((plan[i]['target'], d), 0) + 1
        direction.append(d)
        visit.append(seen[(plan[i]['target'], d)])
    direction = np.array(direction)
    resid_um = resid * FULL_STEP_UM

    # Two terms are expected in that residual before anything is said about
    # repeatability. The first is the lost motion, which flips sign with the
    # approach direction. The second is the screw's once-per-revolution
    # runout: one revolution is 200 full steps = 2 mm, so every target in this
    # phase sits either at a whole revolution (even millimetres) or at a half
    # one (odd millimetres), the quadrature part of the runout is unsampled,
    # and the term collapses to the parity of the millimetre. Fit both, and
    # what is left is the repeatability with no structure left to explain it.
    mm = cmd / 100.0
    parity = (-1.0) ** np.abs(np.rint(mm))
    design = np.column_stack([np.ones(mm.size), direction, parity])
    coef, *_ = np.linalg.lstsq(design, resid_um, rcond=None)
    clean = resid_um - design @ coef
    cov = ((clean @ clean) / (mm.size - design.shape[1])
           * np.linalg.inv(design.T @ design))
    se = np.sqrt(np.diag(cov))

    # Execution order and clock, for the warm-up question: the residual
    # against when the move happened rather than against where it went. The
    # fit is run twice, with and without the first move of the run, because a
    # single leveraged point at one end is exactly what a drift looks like.
    t = np.array([records[i]['move_end'] / FS for i in idx])
    trend = line_fit(t / 60.0, clean)
    trend_tail = line_fit(t[1:] / 60.0, clean[1:])

    mres = int(np.median([plan[i]['mres'] for i in idx]))

    return dict(
        n=len(idx), slope=float(slope), intercept=float(intercept),
        residual_rms_um=float(resid.std() * FULL_STEP_UM),
        residual_max_um=float(np.abs(resid).max() * FULL_STEP_UM),
        counts_per_full_step=float(counts_per_fs),
        counts_per_rev=float(counts_per_fs * FULL_STEPS_PER_REV),
        counts_per_mm=float(counts_per_fs * 100),
        row=[int(i) for i in idx], t_s=t.tolist(),
        direction=direction.tolist(), visit=visit,
        mres=mres, microstep_um=FULL_STEP_UM / mres,
        reversal_um=float(2 * coef[1]), reversal_se_um=float(2 * se[1]),
        screw_once_per_rev_um=float(coef[2]),
        screw_once_per_rev_se_um=float(se[2]),
        clean_resid_um=clean.tolist(),
        clean_rms_um=float(clean.std()), clean_max_um=float(np.abs(clean).max()),
        trend=trend, trend_no_first=trend_tail,
        cmd=cmd.tolist(), obs=obs.tolist(), resid_um=resid_um.tolist())


def gain_moves_result(res, records, plan, q):
    """Each gain move as a relative displacement, not as a position.

    The position figure asks where the carriage ended up against a zero; this
    asks how far it went, which does not involve a zero at all. The two differ
    by exactly the thing the gain fit cannot see: an error that belongs to the
    endpoints rather than to the distance survives here, and a scale error
    shows up as a slope.

    The first move has no previous settled position, so the hold in front of
    it stands in for one; it is the only point here that depends on where the
    stage was left before the run started.
    """
    first = records[res['row'][0]]
    hold = float(np.median(q[max(0, first['move_start'] - 2 * int(FS)):
                             first['move_start']]))

    cmd, obs, occurrence, seen = [], [], [], {}
    prev = hold
    for i in res['row']:
        settled = records[i]['settled']
        cmd.append((plan[i]['target'] - plan[i]['start']) * FULL_STEP_UM)
        obs.append((settled - prev) * FULL_STEP_UM)
        prev = settled
        # Counted by the move, not by its length: the phase contains three
        # moves of 5 mm that run between different pairs of positions, and
        # calling those a repeat would be the opposite of the point.
        leg = (plan[i]['start'], plan[i]['target'])
        seen[leg] = seen.get(leg, 0) + 1
        occurrence.append(seen[leg])
    cmd = np.array(cmd, float)
    obs = np.array(obs, float)

    # Scale, fitted through the origin: a relative move has no offset to fit.
    slope = float(obs @ cmd / (cmd @ cmd))
    err = obs - cmd
    return dict(hold_um=hold * FULL_STEP_UM, slope=slope,
                target_mm=[plan[i]['target'] / 100.0 for i in res['row']],
                scale_error_pct=float((slope - 1) * 100),
                cmd_um=cmd.tolist(), obs_um=obs.tolist(), err_um=err.tolist(),
                occurrence=occurrence,
                direction=[1 if c > 0 else -1 for c in cmd],
                err_rms_um=float(err.std()), err_max_um=float(np.abs(err).max()),
                err_mean_um=float(err.mean()))


SCREW_MM = 2.0            # one screw revolution
SETTLE_WINDOW_S = 5.0     # how much of each dwell the settling traces carry


def settling_result(records, plan, q, gain_um_per_mm):
    """The settling matrix, move by move, with the two errors it supports.

    The raw stopping error of a move is dominated by things the move had no
    part in: where the campaign's origin sits, the step-to-distance gain, and
    the screw's once-per-revolution runout at whatever angle the target
    happens to fall on. All three are properties of the position, not of the
    move.

    The gain matters more here than anywhere else in the campaign, and not
    because it is large. The error is signed along the direction of travel,
    and the 10 mm moves shuttle between the ends of the range, so their stop
    position times their direction averages +5.1 mm where the 0.1 mm moves
    average -0.05. Any gain left in therefore appears as an error that grows
    with commanded distance - a distance effect that is not one. The gain
    used is phase 4's, averaged over the two approach directions so that the
    lost motion cancels out of it, rather than one fitted here, which would
    remove a real distance effect along with the artificial one.

    Two errors come out of this, and they answer different questions:

      settling error - stop minus target, which is where the carriage ended
        up. The lost motion sits in this as a steady offset, the same for
        every move, because a move longer than the 20-odd microsteps of
        take-up ends in the same state whichever way it arrived.

      move error - travel minus commanded distance, which is how far the
        carriage went. A reversing move starts from the other side of the
        backlash and so travels about 2 um less than it was told to; a
        continuing move does not. This is the quantity phases 2 and 4
        measure, and it is where approach history belongs.
    """
    rows, previous = [], None
    for i in sorted(records):
        row, rec = plan[i], records[i]
        direction = 1 if row['target'] > row['start'] else -1
        if (row['phase'] == 'settling' and 'distance=' in row['note']
                and previous is not None):
            note = dict(kv.split('=') for kv in row['note'].split())
            st = settling_times(q, rec)
            commanded = (row['target'] - previous['target']) / 100.0
            travelled = (rec['settled'] - previous['settled']) / 100.0
            rows.append(dict(
                row=int(i),
                distance_mm=float(note['distance'].rstrip('m')),
                speed=note['speed'],
                speed_mm_s=row['speed'] / 100.0,
                direction=direction,
                history='continuation' if direction == previous['direction']
                else 'reversal',
                target_mm=row['target'] / 100.0,
                settled_mm=rec['settled'] / 100.0,
                from_mm=previous['target'] / 100.0,
                commanded_mm=commanded,
                travelled_mm=travelled,
                move_error_um=float((travelled - commanded
                                     * (1 + gain_um_per_mm / 1000.0))
                                    * 1000.0 * direction),
                raw_error_um=(rec['settled'] - rec['target']) * FULL_STEP_UM,
                settle_2um=st[2.0], settle_1um=st[1.0], settle_05um=st[0.5],
                noise_um=rec['noise'] * FULL_STEP_UM,
                t=rec['move_end'] / FS,
                move_start=int(rec['move_start']), move_end=int(rec['move_end']),
                dwell_end=int(rec['dwell_end'])))
        previous = dict(target=row['target'], settled=rec['settled'],
                        direction=direction)

    target = np.array([r['target_mm'] for r in rows])
    direction = np.array([r['direction'] for r in rows])
    residual = (np.array([r['raw_error_um'] for r in rows])
                - gain_um_per_mm * target)

    # What is left of the position error once the gain is out: a constant, a
    # direction term, and the screw. The direction term is fitted alongside
    # the screw so that the lost motion does not leak into the sine, and then
    # put back - it is the result, not a nuisance.
    angle = 2 * np.pi * (target % SCREW_MM) / SCREW_MM
    design = np.column_stack([np.ones(angle.size), direction,
                              np.cos(angle), np.sin(angle)])
    coef, *_ = np.linalg.lstsq(design, residual, rcond=None)
    screw = design[:, 2:] @ coef[2:]
    corrected = residual - coef[0] - screw

    for r, raw, fitted, value in zip(rows, residual, screw, corrected):
        r['residual_um'] = float(raw)
        r['screw_um'] = float(fitted)
        # Signed along the direction of travel, so that stopping short of the
        # target reads negative whichever way the carriage was going.
        r['error_um'] = float(value * r['direction'])

    model = dict(
        gain_um_per_mm=float(gain_um_per_mm),
        gain_source='phase 4, averaged over both approach directions',
        offset_um=float(coef[0]), direction_um=float(coef[1]),
        screw_cos_um=float(coef[2]), screw_sin_um=float(coef[3]),
        screw_amplitude_um=float(np.hypot(coef[2], coef[3])),
        screw_phase_deg=float(np.degrees(np.arctan2(coef[3], coef[2])) % 360),
        residual_rms_um=float(np.std(residual)),
        corrected_rms_um=float(np.std(corrected)),
        settling_mean_um=float(np.mean([r['error_um'] for r in rows])),
        settling_sd_um=float(np.std([r['error_um'] for r in rows], ddof=1)))
    return rows, model


def repeatability_gain(rep):
    """Gain error in µm per mm, from phase 4 with the direction averaged out.

    Each target in phase 4 is approached from both sides the same number of
    times, so the mean of the two directions is free of the lost motion and
    what is left varying with position is the gain.
    """
    target = np.array([p['target_mm'] for p in rep['per_target']])
    middle = np.array([(p['mean_forward'] + p['mean_reverse']) / 2
                       for p in rep['per_target']])
    slope, _ = np.polyfit(target, middle, 1)
    return float(slope)


def settling_traces(rows, q):
    '''Each stop's approach to its final value, for the first five seconds.'''
    out = []
    for r in rows:
        a = r['move_end']
        b = min(r['dwell_end'], a + int(SETTLE_WINDOW_S * FS))
        if b - a < 100:
            continue
        settled = r['settled_mm'] * 100.0
        out.append(dict(speed=r['speed'], history=r['history'],
                        direction=r['direction'],
                        distance_mm=r['distance_mm'],
                        trace=(q[a:b] - settled) * FULL_STEP_UM))
    return out


READ_FROM_S, READ_TO_S = 1.5, 2.5     # of the 2.5 s hold, after the creep


def repeatability_result(records, plan, q):
    """Phase 4 to ISO 230-2, from the plan rather than from the signal.

    Every move here is a row of the plan with a role in its label, so the
    moves are segmented by reading the plan and the recording is used only to
    say when each one happened. That matters because the quantity this phase
    measures is the difference between two approaches to the same point: if
    the segmentation came from the signal, a move too small to detect would
    drop a target rather than report it.

    The settled value is the mean over 1.5 to 2.5 s after the end of motion.
    Phase 3 showed the carriage creeping into place with half the creep gone
    by 0.3 s, so a window that opens at 1.5 s is past it; the standard
    deviation within the window is kept as a quality flag rather than used.

    Deviations are in machine coordinates, not signed along travel: ISO
    230-2 is written that way, and with the two directions reported
    separately the sign carries the direction already. One offset is removed
    from the whole phase - the mean of all 50 - because the campaign's zero
    is arbitrary. The gain is left in: it is one of the results.
    """
    rows, pending = [], None
    for i in sorted(records):
        row, rec = plan[i], records[i]
        if row['phase'] != 'repeatability':
            continue
        role = row['label'].split()[0]
        if role == 'run-up':
            pending = dict(row=int(i), from_mm=row['start'] / 100.0,
                           to_mm=row['target'] / 100.0,
                           move_start=int(rec['move_start']),
                           move_end=int(rec['move_end']))
            continue
        if role != 'target':
            continue
        note = row['note']
        a = rec['move_end'] + int(READ_FROM_S * FS)
        b = min(rec['move_end'] + int(READ_TO_S * FS), rec['dwell_end'])
        if b - a < 300:
            continue
        window = q[a:b]
        rows.append(dict(
            row=int(i), run_up=pending,
            target_mm=row['target'] / 100.0,
            direction='forward' if 'forward' in note else 'reverse',
            cycle=int(note.split('rep ')[1].split('/')[0]),
            settled_mm=float(np.mean(window)) / 100.0,
            window_sd_um=float(np.std(window)) * FULL_STEP_UM,
            move_end=int(rec['move_end']), dwell_end=int(rec['dwell_end']),
            move_start=int(rec['move_start']),
            t=rec['move_end'] / FS))
        pending = None

    offset = np.mean([(r['settled_mm'] - r['target_mm']) for r in rows]) * 1000
    for r in rows:
        r['deviation_um'] = float((r['settled_mm'] - r['target_mm']) * 1000
                                  - offset)

    # ---- how far each hold still had to creep when it began
    creep = {}
    for direction in ('forward', 'reverse'):
        first = [float(q[r['move_end']] / 100.0 - r['settled_mm']) * 1000
                 for r in rows if r['direction'] == direction]
        creep[direction + '_um'] = float(np.median(first))
    creep['half_life_s'] = 0.3        # measured in phase 3, same mechanism

    # ---- drift, before any statistic that assumes there is none
    cycle = np.array([r['cycle'] for r in rows], float)
    deviation = np.array([r['deviation_um'] for r in rows])
    trend = line_fit(cycle, deviation)
    drift = dict(per_cycle_um=trend['slope'], se_um=trend['se'],
                 over_phase_um=trend['slope'] * 4,
                 negligible=bool(abs(trend['slope'] * 4) < 0.1))

    # ---- ISO 230-2, per target and direction
    targets = sorted({r['target_mm'] for r in rows})
    per_target = []
    for target in targets:
        cell = {}
        for direction in ('forward', 'reverse'):
            values = np.array([r['deviation_um'] for r in rows
                               if r['target_mm'] == target
                               and r['direction'] == direction])
            cell[direction] = dict(n=int(values.size),
                                   mean=float(values.mean()),
                                   s=float(values.std(ddof=1)),
                                   values=values.tolist())
        up, down = cell['forward'], cell['reverse']
        bias = up['mean'] - down['mean']
        per_target.append(dict(
            target_mm=target,
            n_forward=up['n'], n_reverse=down['n'],
            mean_forward=up['mean'], mean_reverse=down['mean'],
            s_forward=up['s'], s_reverse=down['s'],
            forward=up['values'], reverse=down['values'],
            reversal=float(bias), mid=float((up['mean'] + down['mean']) / 2),
            R_forward=float(4 * up['s']), R_reverse=float(4 * down['s']),
            R=float(max(2 * up['s'] + 2 * down['s'] + abs(bias),
                        4 * up['s'], 4 * down['s']))))

    bias = np.array([p['reversal'] for p in per_target])
    mid = np.array([p['mid'] for p in per_target])
    position = np.array([p['target_mm'] for p in per_target])
    highs = [p[f'mean_{d}'] + 2 * p[f's_{d}']
             for p in per_target for d in ('forward', 'reverse')]
    lows = [p[f'mean_{d}'] - 2 * p[f's_{d}']
            for p in per_target for d in ('forward', 'reverse')]

    # The gain, from the direction-averaged deviation: the lost motion is
    # equal and opposite in the two directions, so it cancels in the mean and
    # what is left varying with position is the scale.
    gain = line_fit(position, mid)

    # What the gain line does not explain. The screw advances 2 mm a turn, so
    # the targets 5 mm apart sit alternately at a whole turn from zero and at
    # a half turn from it, where the once-per-revolution runout has the
    # opposite sign. If that is what the residual is, it separates by parity.
    residual = mid - (gain['slope'] * position + gain['intercept'])
    turns = position / 2.0
    whole = residual[np.isclose(turns % 1.0, 0)]
    halves = residual[~np.isclose(turns % 1.0, 0)]
    screw = dict(
        whole_turn_um=float(whole.mean()), half_turn_um=float(halves.mean()),
        amplitude_um=float((halves.mean() - whole.mean()) / 2),
        residual_rms_um=float(np.sqrt((residual ** 2).mean())),
        separates=bool(halves.min() > whole.max()),
        per_target={float(p): float(r) for p, r in zip(position, residual)})

    return dict(
        n=len(rows), offset_removed_um=float(offset),
        read_window_s=[READ_FROM_S, READ_TO_S],
        moves=rows, per_target=per_target, drift=drift,
        reversal_max=float(np.max(np.abs(bias))),
        reversal_mean=float(bias.mean()),
        reversal_sd=float(bias.std(ddof=1)),
        lost_motion_per_side_um=float(bias.mean() / 2),
        repeatability_forward=float(max(p['R_forward'] for p in per_target)),
        repeatability_reverse=float(max(p['R_reverse'] for p in per_target)),
        repeatability=float(max(p['R'] for p in per_target)),
        accuracy=float(max(highs) - min(lows)),
        gain_um_per_mm=gain['slope'], gain_se_um_per_mm=gain['se'],
        gain_pct=gain['slope'] / 1000 * 100,
        gain_se_pct=gain['se'] / 1000 * 100,
        gain_intercept_um=gain['intercept'], screw=screw, creep=creep,
        definitions=dict(
            deviation='settled position minus target, machine coordinates, '
                      'with the phase mean removed',
            B='x̄ forward minus x̄ reverse at one target',
            M='mean of the two directions at one target',
            R='max(2s↑ + 2s↓ + |B|, 4s↑, 4s↓) at one target',
            A='max(x̄ + 2s) − min(x̄ − 2s) over all targets and directions',
            gain='straight line through M against target position',
            screw='M minus that line, grouped by whether the target is a '
                  'whole or a half turn of the 2 mm screw from zero'))


def gain_pair_estimate(gain):
    """Phase 1's gain and lost motion, from its ±m pairs.

    Phase 1 visits +m and −m for m = 1..5 mm, the first always approached
    moving positive and the second moving negative. Half their difference is
    k·m + B/2: the scale appears in the slope and the lost motion in the
    intercept, each without the other. It is the only estimate of either that
    phase 1 can make without assuming the other.
    """
    cmd = np.array(gain['cmd']) / 100.0
    obs = np.array(gain['obs']) / 100.0
    magnitudes = sorted({abs(c) for c in cmd if c})
    half, scale = [], []
    for m in magnitudes:
        up = obs[np.isclose(cmd, m)]
        down = obs[np.isclose(cmd, -m)]
        if not up.size or not down.size:
            continue
        half.append((up.mean() - down.mean()) / 2 * 1000)
        scale.append(m)
    fit = line_fit(np.array(scale), np.array(half))
    # The slope is micrometres travelled per millimetre commanded, so the
    # nominal 1000 comes out of it before it is a gain error.
    return dict(gain_um_per_mm=fit['slope'] - 1000.0,
                gain_se_um_per_mm=fit['se'],
                gain_pct=(fit['slope'] - 1000.0) / 1000 * 100,
                gain_se_pct=fit['se'] / 1000 * 100,
                lost_motion_um=2 * fit['intercept'],
                lost_motion_se_um=2 * fit['se_intercept'],
                n_pairs=len(scale))


def campaign_summary(results):
    """The two numbers every phase measured, gathered in one place."""
    stair = [leg['lost_motion_um'] for leg in results['microstep']['per_leg']]
    settle = results['settling']
    move_error = {k: np.array([r['move_error_um'] for r in settle
                               if r['history'] == k])
                  for k in ('continuation', 'reversal')}
    difference = move_error['reversal'].mean() - move_error['continuation'].mean()
    se = float(np.sqrt(move_error['reversal'].var(ddof=1)
                       / move_error['reversal'].size
                       + move_error['continuation'].var(ddof=1)
                       / move_error['continuation'].size))
    rep = results['repeatability']
    pair = results['gain_pair']

    lost = [
        dict(phase='1', what='± pair fit', value=abs(pair['lost_motion_um']),
             error=pair['lost_motion_se_um'], pooled=True),
        dict(phase='2', what='three staircases',
             value=float(np.mean(stair)), error=float(np.std(stair, ddof=1)),
             pooled=True),
        dict(phase='3', what='reversal − continuation',
             value=abs(difference), error=se, pooled=True),
        dict(phase='4', what='mean B over five targets',
             value=abs(rep['reversal_mean']), error=rep['reversal_sd'],
             pooled=True),
    ]
    gain = [
        dict(phase='1', what='plain fit',
             value=(results['gain']['slope'] - 1) * 100,
             error=float('nan'), pooled=False),
        dict(phase='1', what='± pair fit', value=pair['gain_pct'],
             error=pair['gain_se_pct'], pooled=True),
        dict(phase='4', what='direction-averaged', value=rep['gain_pct'],
             error=rep['gain_se_pct'], pooled=True),
    ]
    for rows in (lost, gain):
        inside = [r['value'] for r in rows if r['pooled']]
        mean, sd = float(np.mean(inside)), float(np.std(inside, ddof=1))
        for r in rows:
            r['pooled_mean'], r['pooled_sd'] = mean, sd
    return dict(lost_motion=lost, gain=gain)


def _nan_to_none(a):
    '''NaN is not valid json; None is, and json.dump writes it as null.'''
    return [None if not np.isfinite(v) else float(v) for v in a]


def staircase_result(legs):
    '''Static microstep error and hysteresis, per staircase location.

    Slot k of the grid is commanded k/mres full steps out for the forward leg
    and the same positions back again for the reverse leg, so the two passes
    are compared slot by slot after reversing the second one. The commanded
    ramp itself is removed by a straight line fitted to the mean of the two
    passes: what is left is the error against a perfectly uniform microstep.
    '''
    out = []
    for leg in legs:
        n_up = leg['micro_steps']
        pos = leg['position']
        if pos.size < 2 * n_up:
            pos = np.r_[pos, np.full(2 * n_up + 1 - pos.size, np.nan)]
        # Slot n_up is the turnaround: it is the last forward slot and the
        # first reverse slot at once, so it belongs to neither pass and its
        # plateau straddles the direction change. Both passes stop short of it.
        last_read = int(np.flatnonzero(np.isfinite(pos))[-1])
        up = pos[:n_up]
        down = pos[n_up + 1:2 * n_up + 1][::-1]
        cmd = np.arange(n_up) / leg['mres']

        # Against the commanded ramp, which is where the microstep error is
        # defined: the position the drive was asked for is k microsteps out,
        # and the only freedom left is the constant that says where the leg
        # started. It is taken from the mean of the two passes, so neither
        # direction is privileged.
        cmd_um = cmd * FULL_STEP_UM
        raw_up = (up - pos[0]) * FULL_STEP_UM - cmd_um
        raw_dn = (down - pos[0]) * FULL_STEP_UM - cmd_um
        offset = np.nanmean(np.r_[raw_up, raw_dn])
        err_up_abs = raw_up - offset
        err_dn_abs = raw_dn - offset

        # Each microstep as its own displacement, in the direction of travel,
        # counted from the reversal each pass starts at. The turnaround slot
        # is left out of both: it straddles the direction change.
        step_up = np.diff(pos[:n_up]) * FULL_STEP_UM
        step_dn = -np.diff(pos[n_up + 1:2 * n_up + 1]) * FULL_STEP_UM
        steps = np.r_[step_up, step_dn]

        # Every commanded microstep, including the two either side of the
        # turnaround. Those two are left out of the pass-against-pass
        # comparison because the slot at the apex belongs to both passes at
        # once, but as displacements they are as measurable as any other: the
        # carriage dwells there for a full slot like everywhere else.
        step_all = np.diff(pos) * FULL_STEP_UM
        step_all[n_up:] *= -1
        half_step = 0.5 * FULL_STEP_UM / leg['mres']

        # How much travel the first microsteps after the reversal fail to
        # deliver, and how many of them it takes before they stop failing:
        # the lost motion, counted in the unit the drive commands in.
        nominal_step = FULL_STEP_UM / leg['mres']
        deficit = np.cumsum(np.nan_to_num(nominal_step - step_dn[:30]))
        take_up = int(np.argmax(deficit)) + 1
        # n95: the microstep by which 95 % of the lost motion has been made
        # up. The peak itself is where the deficit stops growing, which the
        # last few hundredths of a micrometre make noisy; n95 is the number
        # to quote as "how long a reversal costs".
        n95 = int(np.argmax(deficit >= 0.95 * deficit[take_up - 1])) + 1
        # The same for the forward pass. It is only non-zero where the leg was
        # entered against its approach, which is the whole point of measuring
        # it: it says whether a reversal after a continuous move costs the
        # same as a reversal inside the staircase.
        deficit_up = np.cumsum(np.nan_to_num(nominal_step - step_up[:30]))
        take_up_fwd = int(np.argmax(deficit_up)) + 1
        n95_fwd = int(np.argmax(
            deficit_up >= 0.95 * deficit_up[take_up_fwd - 1])) + 1

        both = np.nanmean(np.vstack([up, down]), axis=0)
        ok = np.isfinite(both)
        base = np.polyfit(cmd[ok], both[ok], 1)
        ref = np.polyval(base, cmd)
        err_up = (up - ref) * FULL_STEP_UM
        err_dn = (down - ref) * FULL_STEP_UM
        with np.errstate(invalid='ignore'):
            mean_err = np.nanmean(np.vstack([err_up, err_dn]), axis=0)
        hysteresis = err_up - err_dn
        out.append(dict(
            location_mm=leg['location'] / 100.0, n=int(n_up),
            period_s=leg['period_s'], slots=leg['slots'], mres=leg['mres'],
            expected_slots=leg['expected_slots'],
            missing=int((~np.isfinite(up)).sum() + (~np.isfinite(down)).sum()),
            cmd_full_steps=cmd.tolist(),
            err_up_um=np.where(np.isfinite(err_up), err_up, None).tolist(),
            err_dn_um=np.where(np.isfinite(err_dn), err_dn, None).tolist(),
            mean_err_um=np.where(np.isfinite(mean_err), mean_err, None).tolist(),
            hysteresis_um=np.where(np.isfinite(hysteresis), hysteresis,
                                   None).tolist(),
            phase=((cmd % 4.0) / 4.0).tolist(),
            gain_in_leg=float(base[0]),
            # The same two passes against the commanded ramp itself, with
            # only a constant taken out: no per-leg slope is fitted, so a leg
            # that ran long or short still shows it, and the pair keeps the
            # lost motion between them instead of sharing it out.
            cmd_index=list(range(n_up)),
            # The whole leg as measured, slot by slot in the order it
            # happened, relative to the slot it started from. The turnaround
            # slot is in here: it is excluded from the two passes because it
            # belongs to neither, but the sequence ran through it.
            position_um=_nan_to_none((pos - pos[0]) * FULL_STEP_UM),
            # Where the leg finished against where it was told to finish.
            # The sequence returns to its starting command, so anything left
            # here is lost motion that was never given back.
            closing_error_um=float(
                (pos[last_read] - pos[0]) * FULL_STEP_UM
                - (2 * n_up - last_read) * FULL_STEP_UM / leg['mres']),
            closing_slot=int(last_read),
            # Where the carriage stood against its command when the leg
            # began. The cumulative-error figure starts every leg at zero,
            # so this is the one number that figure removes: it belongs to
            # the approach move, not to the staircase.
            start_offset_um=float((pos[0] - leg['location']) * FULL_STEP_UM),
            err_up_abs_um=_nan_to_none(err_up_abs),
            err_dn_abs_um=_nan_to_none(err_dn_abs),
            offset_um=float(offset),
            step_nominal_um=FULL_STEP_UM / leg['mres'],
            step_up_um=_nan_to_none(step_up),
            step_dn_um=_nan_to_none(step_dn),
            step_median_um=float(np.nanmedian(steps)),
            step_sd_um=float(np.nanstd(steps)),
            step_min_um=float(np.nanmin(steps)),
            step_max_um=float(np.nanmax(steps)),
            step_stalled=int(np.nansum(steps < FULL_STEP_UM / leg['mres'] / 3)),
            # How many of the commanded microsteps actually show up as a
            # distinct movement. A microstep that moves less than half of
            # one is not resolvable from the step either side of it, and
            # those are concentrated right after a reversal.
            steps_commanded=2 * n_up,
            steps_measured=int(np.sum(np.isfinite(step_all))),
            steps_resolved=int(np.nansum(step_all > half_step)),
            # Two different failures, reported separately: a microstep that
            # moved the carriage less than half of one, and a microstep that
            # left no step the signal can be seen to take. The first is a
            # small step; the second is no visible step at all.
            steps_small=[int(k + 1) for k in
                         np.flatnonzero(step_all <= half_step)],
            steps_undetected=list(leg['undetected']),
            step_all_um=_nan_to_none(step_all),
            step_mean_um=float(np.nanmean(step_all)),
            step_sd_all_um=float(np.nanstd(step_all)),
            step_pv_um=float(np.nanmax(step_all) - np.nanmin(step_all)),
            reversal_take_up_steps=take_up, reversal_take_up_n95=n95,
            reversal_take_up_um=float(deficit[take_up - 1]),
            forward_take_up_steps=take_up_fwd, forward_take_up_n95=n95_fwd,
            forward_take_up_um=float(deficit_up[take_up_fwd - 1]),
            grid_shift=leg['grid_shift'], slots_read=leg['read'],
            static_pv_um=float(np.nanmax(mean_err) - np.nanmin(mean_err)),
            static_rms_um=float(np.sqrt(np.nanmean(mean_err ** 2))),
            hysteresis_mean_um=float(np.nanmean(hysteresis)),
            hysteresis_max_um=float(np.nanmax(np.abs(hysteresis)))))
    return out


HARMONICS = (1, 2, 4)      # cycles per electrical period: table, imbalance, detent
TRANSITION = 25            # microsteps after a reversal, excluded from the fits


def _line(k, y):
    '''Slope and intercept of y against k, over the finite points only.'''
    ok = np.isfinite(y) & np.isfinite(k)
    if ok.sum() < 10:
        return float('nan'), float('nan')
    slope, intercept = np.polyfit(k[ok], y[ok], 1)
    return float(slope), float(intercept)


def _harmonic_fit(k, y, per_period, harmonics=HARMONICS):
    '''Line plus harmonics, fitted together.

    Fitted together rather than detrended first: a leg covers two electrical
    periods, and over that span a straight line and a once-per-period sine
    are not independent. Taking the line out first would move part of the
    first harmonic into the trend and report the harmonic too small.
    '''
    ok = np.isfinite(y)
    cols = [np.ones(k.size), k.astype(float)]
    for h in harmonics:
        cols += [np.cos(2 * np.pi * h * k / per_period),
                 np.sin(2 * np.pi * h * k / per_period)]
    design = np.column_stack(cols)
    coef, *_ = np.linalg.lstsq(design[ok], y[ok], rcond=None)
    trend = design[:, :2] @ coef[:2]
    # String keys: this dict goes to json and comes back, and an int key
    # would not survive the round trip.
    amp = {str(h): float(np.hypot(coef[2 + 2 * n], coef[3 + 2 * n]))
           for n, h in enumerate(harmonics)}
    return dict(intercept=float(coef[0]), slope=float(coef[1]), amplitude=amp,
                detrended=y - trend, residual=y - design @ coef)


def microstep_analysis(res, legs, drift_rows, reversal_um):
    '''Drift, scale, step statistics and the folded microstep map.

    Three things are tangled in a staircase leg: the commanded ramp, a scale
    error that makes every microstep slightly the wrong size, and anything
    that moves the stage while the leg runs. The forward pass runs with time
    and the reverse pass against it, so a term that is the same in both is
    scale and a term that flips sign is drift. Nothing else in the campaign
    separates the two.
    '''
    per_leg, folded_rows = [], []
    for leg, raw in zip(res, legs):
        n, mres = leg['n'], leg['mres']
        nominal = leg['step_nominal_um']
        per_period = 4 * mres
        k = np.arange(n)
        ef = np.array([np.nan if v is None else v for v in leg['err_up_abs_um']])
        er = np.array([np.nan if v is None else v for v in leg['err_dn_abs_um']])

        # ---- 3a, 3b: scale and drift, from the two slopes
        #
        # The slope comes out of the joint fit below, not out of a straight
        # line through the errors. A leg spans two electrical periods and the
        # fit range, with the transition removed, spans about 1.6 of them -
        # too few for a line and a once-per-period sine to be independent. A
        # plain line fit takes part of that harmonic into its slope, and at
        # +0 mm that was enough to report the local step 2 % long and then
        # subtract the difference from the lost motion.
        fwd_ok = np.where(k >= TRANSITION, ef, np.nan)
        rev_ok = np.where(k <= n - 1 - TRANSITION, er, np.nan)
        fit_f = _harmonic_fit(k, fwd_ok, per_period)
        fit_r = _harmonic_fit(k, rev_ok, per_period)
        s_f, a_f = fit_f['slope'], fit_f['intercept']
        s_r, a_r = fit_r['slope'], fit_r['intercept']
        scale_err = (s_f + s_r) / 2                      # µm per microstep
        drift_um_s = (s_f - s_r) / 2 / leg['period_s']   # µm per second

        # ---- 3c: the hysteresis with that drift taken out
        #
        # Commanded microstep k is visited at time k.tau going out and
        # (2n - k).tau coming back, so with a drift d the measured errors are
        #     e_f(k) = true_f(k) + d.k.tau
        #     e_r(k) = true_r(k) + d.(2n - k).tau
        # and their difference at k = n/2, where the two visits are n slots
        # apart, is
        #     measured gap = true gap - d.n.tau.
        # The reverse pass is the later one, so an upward drift lifts it and
        # makes the measured gap too small; the correction is therefore
        # added, not subtracted. The same algebra gives the slopes used
        # above: s_f = s_true + d.tau and s_r = s_true - d.tau, so their mean
        # is the scale and half their difference is the drift.
        mid = n / 2
        separation_s = n * leg['period_s']
        gap = (a_f + s_f * mid) - (a_r + s_r * mid)
        drift_correction = drift_um_s * separation_s

        # The gap is the lost motion, and it is the number to quote. Both
        # passes run over the same ground, so the local scale error is in
        # both of them and cancels in the difference; the cumulative deficit
        # below has that scale in it once, which tilts its plateau. Three
        # separate things move this number away from the campaign's averaged
        # one and they are kept apart on purpose: the definition (fitted
        # lines at one k, with the transition and the harmonics out, against
        # a mean taken over every microstep including the ones either side of
        # a reversal), and then the drift.
        lost_gap = abs(gap + drift_correction)
        definition_shift = gap - leg['hysteresis_mean_um']

        # ---- 4a: the shape of the transition, from the cumulative deficit
        #
        # D(n) is referred to the nominal microstep, which is a fixed number.
        # It is not where the lost motion is read from - the gap above is -
        # because the local scale enters D(n) once and tilts its plateau: at
        # +0 mm the step runs 2.2 % long, which walks the plateau down by
        # 0.69 um over the window, and 1.93 + 0.69 is the 2.61 um the gap
        # reports. What D(n) is good for is the shape of the recovery, and
        # that is read at fractions of the gap value rather than at 95 % of
        # its own plateau, which is the part the tilt has eaten.
        local_step = nominal + scale_err
        reversals = {}
        for name, key in (('forward', 'step_up_um'), ('reverse', 'step_dn_um')):
            step = np.array([np.nan if v is None else v for v in leg[key]])
            deficit = np.cumsum(np.nan_to_num(nominal - step))
            against_local = np.cumsum(np.nan_to_num(local_step - step))
            plateau = deficit[39:60]
            lost = float(np.mean(plateau))

            # Recovery, measured against the gap rather than the plateau.
            # n50 and n80 are reached while the tilt is still small; the
            # running-mean rule does not use D(n) at all - it asks when the
            # steps themselves came back up to the local size.
            def reach(frac):
                hit = np.flatnonzero(deficit >= frac * lost_gap)
                return int(hit[0] + 1) if hit.size else 0

            window = 8
            recovered = 0
            for m in range(window, step.size):
                if np.nanmean(step[m - window:m]) >= 0.98 * local_step:
                    recovered = m
                    break
            # Kept for comparison, and not to be quoted: 95 % of a plateau
            # the tilt has already moved is the most tail-sensitive number
            # available here.
            n95 = int(np.argmax(deficit >= 0.95 * lost) + 1) if lost > 0 else 0
            # The deficit telescopes: D(n) = n x nominal - (pos_n - pos_0), so
            # its uncertainty is two position read-outs, not n of them. The
            # per-step scatter is real motion and cancels in the sum.
            read_noise = float(np.nanmedian(raw['slot_noise']) * FULL_STEP_UM)
            reversals[name] = dict(
                plateau_um=lost, n50=reach(0.5), n80=reach(0.8), n95=n95,
                n_steps_recovered=recovered,
                plateau_sd_um=float(np.std(plateau)),
                read_noise_um=float(np.sqrt(2) * read_noise),
                initial_step_ratio=float(np.nanmean(step[:3]) / nominal),
                plateau_vs_local_step_um=float(np.mean(against_local[39:60])),
                plateau_tilt_um=float(deficit[59] - deficit[39]),
                deficit_um=deficit.tolist())

        # ---- 4b: step size away from the reversal
        steps = {}
        for name, key in (('forward', 'step_up_um'), ('reverse', 'step_dn_um')):
            step = np.array([np.nan if v is None else v for v in leg[key]])
            tail = step[TRANSITION + 5:]
            tail = tail[np.isfinite(tail)]
            lag1 = float(np.corrcoef(tail[:-1], tail[1:])[0, 1])
            steps[name] = dict(
                n=int(tail.size), median_um=float(np.median(tail)),
                mean_um=float(tail.mean()), sd_um=float(tail.std()),
                min_um=float(tail.min()), max_um=float(tail.max()),
                stuck=int(np.sum(tail < 0.3)), jumps=int(np.sum(tail > 0.95)),
                lag1_autocorrelation=lag1)

        # ---- 4c: the line and the three harmonics, fitted together
        harmonic = {}
        for name, fit, keep in (('forward', fit_f, k >= TRANSITION),
                                ('reverse', fit_r, k <= n - 1 - TRANSITION)):
            harmonic[name] = dict(amplitude=fit['amplitude'],
                                  residual_rms_um=float(
                                      np.nanstd(fit['residual'])))
            for kk, value in zip(k[keep], fit['detrended'][keep]):
                if np.isfinite(value):
                    folded_rows.append((leg['location_mm'], name,
                                        int(kk % per_period), float(value)))

        per_leg.append(dict(
            location_mm=leg['location_mm'],
            slope_forward_um_per_step=s_f, slope_reverse_um_per_step=s_r,
            scale_error_um_per_step=scale_err,
            scale_error_pct=float(scale_err / nominal * 100),
            local_step_um=float(nominal + scale_err),
            drift_um_per_s=drift_um_s, drift_um_per_min=drift_um_s * 60,
            drift_correction_applied_um=float(drift_correction),
            separation_s=float(separation_s),
            # Three numbers for one quantity, kept apart so the shift from
            # one to the next can be attributed rather than lumped together.
            hysteresis_mean_um=float(leg['hysteresis_mean_um']),
            hysteresis_fitted_um=float(gap),
            definition_shift_um=float(definition_shift),
            hysteresis_drift_corrected_um=float(gap + drift_correction),
            lost_motion_um=float(lost_gap),
            drift_negligible=bool(abs(drift_correction) < 0.2),
            reversal=reversals, steps=steps, harmonic=harmonic))

    # ---- 4c, 4e: the folded map, and how alike the locations are
    per_period = 4 * res[0]['mres']
    angle = np.arange(per_period) * 360.0 / per_period
    folded = dict(angle_deg=angle.tolist(), per_period=per_period)

    def fold(rows):
        grid = [[] for _ in range(per_period)]
        for _, _, phase, value in rows:
            grid[phase].append(value)
        mean = np.array([np.mean(g) if g else np.nan for g in grid])
        mean -= np.nanmean(mean)
        return mean, np.array([np.std(g) if g else np.nan for g in grid]), \
            np.array([len(g) for g in grid])

    mean, sd, count = fold(folded_rows)
    fit = _harmonic_fit(np.arange(per_period), mean, per_period)
    folded.update(mean_um=mean.tolist(), sd_um=sd.tolist(),
                  n_per_angle=count.tolist(),
                  rms_um=float(np.nanstd(mean)),
                  pv_um=float(np.nanmax(mean) - np.nanmin(mean)),
                  amplitude=fit['amplitude'],
                  harmonic_curve_um=(fit['detrended']
                                     - fit['residual']).tolist())
    shapes = {}
    for name in ('forward', 'reverse'):
        m, sdev, _ = fold([r for r in folded_rows if r[1] == name])
        shapes[name] = m
        folded[name] = dict(mean_um=m.tolist(), sd_um=sdev.tolist(),
                            rms_um=float(np.nanstd(m)))
    # Same shape in both directions means the error is electrical - the drive
    # puts the carriage in the same wrong place whichever way it arrives.
    # A shape that changes with direction is friction shaping it.
    ok = np.isfinite(shapes['forward']) & np.isfinite(shapes['reverse'])
    folded['direction_correlation'] = float(
        np.corrcoef(shapes['forward'][ok], shapes['reverse'][ok])[0, 1])

    by_location = {}
    for leg in res:
        m, _, _ = fold([r for r in folded_rows if r[0] == leg['location_mm']])
        by_location[f"{leg['location_mm']:+.0f}"] = m
    pairs = {}
    names = list(by_location)
    for i, x in enumerate(names):
        for y in names[i + 1:]:
            ok = np.isfinite(by_location[x]) & np.isfinite(by_location[y])
            pairs[f'{x} vs {y}'] = float(
                np.corrcoef(by_location[x][ok], by_location[y][ok])[0, 1])
    folded['by_location_um'] = {k: v.tolist() for k, v in by_location.items()}
    folded['location_correlation'] = pairs

    # ---- 3d: the independent drift rate, from the reference returns
    #
    # The two checks that bracket the staircase are approached from opposite
    # sides, so the raw difference between them is mostly lost motion: they
    # have to be brought into one frame first. The reversal value measured in
    # phase 4 is what does it, and because that correction is larger than the
    # drift being looked for, the same-direction fit below is the one to
    # trust and this pair is the sanity check.
    checks = sorted(drift_rows, key=lambda r: r['t'])
    bracket = [r for r in checks if r['t'] < 1000][:2]
    cross = None
    if len(bracket) == 2:
        first, second = bracket
        same = first['direction'] == second['direction']
        e1 = second['error_um'] + (0.0 if same else reversal_um)
        cross = dict(
            t0=first['t'], t1=second['t'],
            error0_um=first['error_um'], error1_um=second['error_um'],
            approached_alike=same, reversal_applied_um=0.0 if same else reversal_um,
            rate_um_per_min=float((e1 - first['error_um'])
                                  / (second['t'] - first['t']) * 60))

    forward = [r for r in checks if r['direction'] > 0]
    same_direction = None
    if len(forward) >= 3:
        t = np.array([r['t'] for r in forward])
        e = np.array([r['error_um'] for r in forward])
        fit = line_fit(t / 60.0, e)
        same_direction = dict(n=len(forward), rate_um_per_min=fit['slope'],
                              se_um_per_min=fit['se'],
                              t=t.tolist(), error_um=e.tolist())

    return dict(per_leg=per_leg, folded=folded, drift_cross_check=cross,
                drift_same_direction=same_direction,
                transition_excluded=TRANSITION)


def noise_result(counts, align_info):
    '''Standstill behaviour, from the energised hold before the run starts.

    Short-term noise and slow drift are reported separately. Lumping them
    gives a number that is neither: over a 77 s hold the drift dominates, and
    quoting that as "encoder noise" would overstate it by most of an order of
    magnitude.
    '''
    first = int(align_info['seg_s'][0])
    rest = counts[:first] if first > 5000 else counts[:5000]
    t = np.arange(rest.size) / FS
    trend = np.polyfit(t, rest, 1)
    detrended = rest - np.polyval(trend, t)

    # Short-term: scatter within 1 s windows, so the slow term cannot enter.
    whole = (rest.size // 1000) * 1000
    per_second = rest[:whole].reshape(-1, 1000)
    short = float(np.mean(per_second.std(axis=1)))

    nper = 4096
    win = np.hanning(nper)
    blocks = detrended.size // nper
    psd = np.zeros(nper // 2 + 1)
    for k in range(blocks):
        psd += np.abs(np.fft.rfft(detrended[k * nper:(k + 1) * nper] * win)) ** 2
    psd /= max(blocks, 1)
    freq = np.fft.rfftfreq(nper, 1 / FS)
    band = freq > 1.0
    order = np.argsort(psd[band])[::-1]
    peaks, seen = [], []
    for k in order:
        f = float(freq[band][k])
        if any(abs(f - s) < 3 for s in seen):
            continue
        seen.append(f)
        peaks.append(round(f, 1))
        if len(peaks) == 5:
            break
    return dict(n=int(rest.size), seconds=float(rest.size / FS),
                short_term_counts=short, short_term_nm=short,
                detrended_std_counts=float(detrended.std()),
                ptp_counts=float(np.ptp(rest)),
                drift_counts_per_s=float(trend[0]),
                dominant_hz=peaks,
                psd_freq=freq.tolist(), psd=psd.tolist())


# --------------------------------------------------------------------------
# figures
# --------------------------------------------------------------------------

PHASE_LABEL = {
    'gain': '1 · gain',
    'staircase': '2 · staircase',
    'settling': '3 · settling',
    'repeatability': '4 · repeatability',
}


def phase_runs(plan, records):
    """Every planned row placed on the recording clock, merged into phase runs.

    Only part of the plan leaves a detectable segment, so bands drawn from the
    matched rows alone stop short: the staircase legs are quasi-static, the
    alignment skips them, and the longest stretch of phase 2 would fall outside
    every band. The matched rows are anchors between the plan's clock and the
    rig's - the two differ by minutes by the end of the run, but never by much
    between neighbouring anchors - and the rows in between are placed by
    interpolating on that anchor list, exact at every anchor and monotone in
    between. The runs that come out tile the run with no gaps, as planned.
    """
    anchor_plan, anchor_rec = [], []
    for i in sorted(records):
        anchor_plan += [plan[i]['t_start'], plan[i]['t_end']]
        anchor_rec += [records[i]['move_start'] / FS, records[i]['dwell_end'] / FS]
    edges = np.interp([r['t_start'] for r in plan] + [plan[-1]['t_end']],
                      anchor_plan, anchor_rec)

    runs = []
    for i, r in enumerate(plan):
        if runs and runs[-1][0] == r['phase']:
            runs[-1][2] = edges[i + 1]
        else:
            runs.append([r['phase'], edges[i], edges[i + 1]])
    # Butt each run against the next so rounding cannot leave a seam of
    # background between two phases that ran back to back.
    for a, b in zip(runs, runs[1:]):
        a[2] = b[1]
    return runs


def fig_overview(q, plan, records, out):
    t = np.arange(q.size) / FS
    fig, (ax, ax2) = plt.subplots(
        2, 1, figsize=(vs.FULL, 3.9), height_ratios=(2.2, 1), sharex=True,
        constrained_layout=True)

    runs = phase_runs(plan, records)
    labelled = set()
    for phase, lo, hi in runs:
        ax.axvspan(lo, hi, color=vs.PHASE_COLOR[phase], lw=0,
                   alpha=0.10 if phase != 'drift' else 0.22)
        # A drift check chops a phase into several runs; label the first one
        # only, or the overview reads as though the phase restarts each time.
        if phase != 'drift' and hi - lo > 120 and phase not in labelled:
            # Staggered: the gain phase is 214 s of a 2700 s axis, so its
            # label and the staircase's would otherwise run together.
            height = 0.965 if len(labelled) % 2 == 0 else 0.86
            labelled.add(phase)
            ax.text((lo + hi) / 2, height, PHASE_LABEL[phase], ha='center',
                    va='top', fontsize=8, fontweight='bold',
                    color=vs.PHASE_COLOR[phase],
                    transform=ax.get_xaxis_transform())

    ax.axhline(0, color=vs.INK_MUTED, lw=0.8, ls='--', alpha=0.8, zorder=0)
    ax.plot(t, q / 100.0, color=vs.INK, lw=0.5)

    drift = sorted(i for i in records if plan[i]['phase'] == 'drift')
    dt = np.array([records[i]['move_end'] / FS for i in drift])
    de = np.array([(records[i]['settled'] - records[i]['target']) * FULL_STEP_UM
                   for i in drift])
    ax.plot(dt, q[[records[i]['move_end'] for i in drift]] / 100.0, 'v', ms=5.5,
            color=vs.PHASE_COLOR['drift'], zorder=5)

    # The plan starts where its first move does: the recording opens on an
    # energised standstill hold and the operator jog that set the stage up,
    # neither of which is part of the sequence. They have their own figure.
    t0 = runs[0][1]
    t1 = max(rec['dwell_end'] for rec in records.values()) / FS
    pad = 0.004 * (t1 - t0)
    ax.set_xlim(t0 - pad, t1 + pad)
    ax.set_ylabel('position [mm]')
    vs.headline(fig, f'Stage identification run, as executed — '
                     f'{(t1 - t0) / 60:.0f} min, {len(plan)} commanded moves')
    # Headroom above the ±10.5 mm extremes so the phase labels sit in clear
    # space rather than across the repeatability peaks.
    ax.set_ylim(-11.5, 12.9)
    vs.panel(ax, 'a')
    vs.despine(ax)

    ax2.plot(dt, de, 'o-', color=vs.SERIES[0], ms=5, mec=vs.SURFACE, mew=1.2)
    for x, y in zip(dt, de):
        ax2.annotate(f'{y:+.1f}', (x, y), textcoords='offset points',
                     xytext=(0, 8), ha='center', fontsize=7,
                     color=vs.INK_SOFT)
    ax2.axhline(0, color=vs.INK_MUTED, lw=0.8)
    ax2.set_ylabel('reference return\nerror [µm]')
    ax2.set_xlabel('time within the recording [s]')
    vs.panel(ax2, 'b')
    vs.despine(ax2)
    return vs.save(fig, out, 'stage_01_overview')


def fig_gain_moves(res, out):
    """Executed minus commanded, one column per move, in execution order.

    Columns rather than points because the quantity is one number per move
    and the moves are a sequence, not a continuum: the eye can follow the
    first pass across and find the second one repeating it column for column,
    which is the whole claim.
    """
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(vs.FULL, 2.8),
                               constrained_layout=True)
        cmd = np.array(res['cmd_um']) / 1000.0
        err = np.array(res['err_um'])
        direction = np.array(res['direction'])
        occurrence = np.array(res['occurrence'])
        x = np.arange(cmd.size)

        for d, colour in ((1, vs.FORWARD), (-1, vs.REVERSE)):
            for v, filled in ((1, True), (2, False)):
                k = (direction == d) & (occurrence == v)
                if k.any():
                    ax.bar(x[k], err[k], width=0.72, zorder=2,
                           color=colour if filled else vs.SURFACE,
                           edgecolor=colour, linewidth=0.7)

        # Where the second pass starts: the first move back to the target the
        # phase opened on. Not the first repeated leg - the move that gets
        # there is a 6 mm one that happens only once, and the divider would
        # land a column late.
        target = np.array(res['target_mm'])
        repeat = int(np.flatnonzero(target[1:] == target[0])[0] + 1)
        ax.axvline(repeat - 0.5, color=vs.INK_MUTED, lw=0.6, ls=':', zorder=1)
        ax.text(repeat - 0.35, 5.2, 'second pass', fontsize=7,
                color='#52514e', va='top')

        ax.axhline(0, color=vs.INK_MUTED, lw=0.6, zorder=3)
        ax.grid(False)
        ax.set_xlim(-0.8, cmd.size - 0.2)
        ax.set_ylim(-5.5, 5.5)
        ax.set_xticks(x)
        ax.set_xticklabels([f'{c:+.0f}' for c in cmd], fontsize=6)
        ax.tick_params(axis='x', length=0, pad=1.5)
        ax.set_yticks(np.arange(-4, 5, 2))
        ax.set_xlabel('commanded move [mm], in the order executed')
        ax.set_ylabel('executed − commanded [µm]')
        vs.headline(fig, 'Move length is right to 1 part in 10⁵: the ±4 µm '
                         'is in where it stops')
        vs.despine(ax)

        fig.legend(handles=[
            Patch(facecolor=vs.FORWARD, edgecolor=vs.FORWARD, label='move +ve'),
            Patch(facecolor=vs.REVERSE, edgecolor=vs.REVERSE, label='move −ve'),
            Patch(facecolor=vs.INK_MUTED, edgecolor=vs.INK_MUTED,
                  label='1st time'),
            Patch(facecolor=vs.SURFACE, edgecolor=vs.INK_MUTED,
                  label='2nd time')],
            loc='outside lower center', ncol=4, handletextpad=0.25,
            columnspacing=1.6, borderpad=0.1,
            handlelength=1.1, handleheight=0.8)
        return vs.save(fig, out, 'stage_02c_gain_moves')


def _decay_tau(t, y):
    """Time constant of a single exponential, by least squares on log |y|.

    Only the part of the curve that is above the standstill noise is fitted;
    past that the log of the magnitude is the log of the noise and would drag
    the slope towards zero.
    """
    k = (np.abs(y) > 0.02) & (t > 0) & (t < 1.0)
    if k.sum() < 20:
        return None
    slope, _ = np.polyfit(t[k], np.log(np.abs(y[k])), 1)
    return None if slope >= 0 else -1.0 / slope


# Both phase 1 figures are drawn at the size they are meant to be printed -
# one journal column, 8 x 5 cm - rather than at a comfortable screen size and
# shrunk later, because type that is set at the wrong size is the one thing
# scaling cannot fix. The extra dpi is what keeps them sharp on the page.
COLUMN = (3.15, 1.97)
COLUMN_DPI = 400
COLUMN_RC = {
    'font.size': 6.5, 'axes.labelsize': 6.5, 'axes.titlesize': 7.4,
    'xtick.labelsize': 6.0, 'ytick.labelsize': 6.0, 'legend.fontsize': 6.0,
    'axes.linewidth': 0.6, 'xtick.major.width': 0.6, 'ytick.major.width': 0.6,
    'xtick.major.size': 2.5, 'ytick.major.size': 2.5,
}


def fig_gain_raw(res, moves, out):
    """The error at each stop before anything is fitted away, and the fit.

    This is the step in front of the residual figure. The zero here is the
    hold the run started from - the position the operator's ZERO was issued
    at - rather than the offset the campaign fits against the commanded
    targets, because the question "how far from where it was told to go did
    it stop" is asked in the frame the command was given in. The fitted line
    is what the next figure removes: its height is that zero being a few
    micrometres off, and its tilt is the gain error.
    """
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(vs.FULL, 2.6),
                               constrained_layout=True)
        cmd = np.array(res['cmd']) / 100.0
        # Measured minus commanded, in the frame the command was given in.
        raw = (np.array(res['obs']) - np.array(res['cmd'])) * FULL_STEP_UM             - moves['hold_um']
        direction = np.array(res['direction'])

        slope_um_per_mm = (res['slope'] - 1) * 1000.0
        offset_um = res['intercept'] * FULL_STEP_UM - moves['hold_um']
        span = np.array([cmd.min() - 0.5, cmd.max() + 0.5])
        for d, colour, name in ((1, vs.FORWARD, '+ve'), (-1, vs.REVERSE, '−ve')):
            k = direction == d
            ax.plot(cmd[k], raw[k], 'o', color=colour, ms=3.4, mew=0,
                    ls='none', zorder=3, label=f'approached {name}')
        ax.plot(span, slope_um_per_mm * span + offset_um, color=vs.SERIES[2],
                lw=1.2, zorder=2, label='fitted line')

        ax.grid(True, axis='y', color=vs.GRID, lw=0.5)
        ax.set_axisbelow(True)
        ax.set_xlim(*span)
        ax.set_xticks(np.arange(-5, 6))
        ax.set_ylim(0, 7.5)
        ax.set_yticks([0, 2.5, 5.0, 7.5])
        ax.set_xlabel('commanded position [mm]')
        ax.set_ylabel('measured − commanded [µm]')
        vs.headline(fig, 'Phase 1: raw error at each stop, with the fitted '
                         'line')
        vs.despine(ax)

        fig.legend(loc='outside lower center', ncol=3, handletextpad=0.25,
                   columnspacing=1.6, borderpad=0.1)
        return vs.save(fig, out, 'stage_02a_raw_error')


def fig_gain(res, out):
    """Where the stage stops: the residual after the gain fit, nothing else.

    No shading and no fitted curves. The structure in the points - the
    alternation with odd and even millimetres, and the two long pairs at
    +/-1 mm - is left for the reader to see and for the text to explain, so
    the figure is not making a claim the eye cannot check.
    """
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(vs.FULL, 2.8),
                               constrained_layout=True)
        cmd = np.array(res['cmd']) / 100.0
        resid = np.array(res['resid_um'])
        direction = np.array(res['direction'])
        visit = np.array(res['visit'])

        # Ordered +1, -1, +2, -2 ... rather than along the axis, so the two
        # ends of each move sit side by side and the reader is comparing the
        # pair that was measured together. The return to zero is left out:
        # it is one move with no counterpart, and it is still in every number
        # quoted for the phase.
        order = [m * s for m in range(1, 6) for s in (1, -1)]
        slot = {x: j for j, x in enumerate(order)}
        keep = np.array([c in slot for c in cmd])
        pos = np.array([slot.get(c, -1) for c in cmd], float)

        for x in order:
            for d in (1, -1):
                k = (cmd == x) & (direction == d)
                if k.sum() > 1:
                    ax.plot([slot[x], slot[x]], [resid[k].min(), resid[k].max()],
                            color=vs.INK_MUTED, lw=0.6, alpha=0.6, zorder=3)

        for d, marker, colour in ((1, '^', vs.FORWARD), (-1, 'v', vs.REVERSE)):
            for v, filled in ((1, True), (2, False)):
                k = keep & (direction == d) & (visit == v)
                if k.any():
                    ax.plot(pos[k], resid[k], marker, ms=4.2, color=colour,
                            ls='none', mfc=colour if filled else 'none',
                            mew=0.9, zorder=4)

        ax.axhline(0, color=vs.INK_MUTED, lw=0.6, zorder=2)
        ax.grid(True, axis='y', color=vs.GRID, lw=0.5, zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_ylim(-3.5, 3.5)
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels([f'{x:+.0f}' for x in order])
        ax.set_yticks(np.arange(-3, 4))
        ax.set_xlabel('commanded position [mm], paired by size')
        ax.set_ylabel('residual [µm]')
        vs.headline(fig, 'Residual after gain calibration: ±3 µm over 10 mm')
        vs.despine(ax)

        mark = dict(ls='none', ms=4.2, mew=0.9)
        fig.legend(handles=[
            Line2D([], [], marker='^', color=vs.FORWARD,
                   label='approached +ve', **mark),
            Line2D([], [], marker='v', color=vs.REVERSE,
                   label='approached −ve', **mark),
            Line2D([], [], marker='o', color=vs.INK_MUTED, label='1st visit',
                   **mark),
            Line2D([], [], marker='o', color=vs.INK_MUTED, mfc='none',
                   label='2nd visit', **mark)],
            loc='outside lower center', ncol=4, handletextpad=0.25,
            columnspacing=1.6, borderpad=0.1)
        return vs.save(fig, out, 'stage_02_gain')


def gain_dwell_result(res, records, q):
    """The dwells of phase 1 reduced to a median and an envelope per direction.

    Twenty-one traces on one axis is a blob. The median says where the stage
    is and the envelope says how much that varied, which is the whole of what
    the figure has to carry.
    """
    span = int(np.median([records[i]['dwell_end'] - records[i]['move_end']
                          for i in res['row']]))
    out = dict(span_s=span / FS, t=((np.arange(span) + 1) / FS).tolist(),
               read_from_s=0.6 * span / FS, by_direction={})
    for want in (1, -1):
        block = []
        for i, d in zip(res['row'], res['direction']):
            if d != want:
                continue
            rec = records[i]
            trace = ((q[rec['move_end']:rec['dwell_end']] - rec['settled'])
                     * FULL_STEP_UM)
            padded = np.full(span, np.nan)
            n = min(span, trace.size)
            padded[:n] = trace[:n]
            block.append(padded)
        block = np.array(block)
        med = np.nanmedian(block, axis=0)
        out['by_direction'][want] = dict(
            n=len(block), median=med.tolist(),
            lo=np.nanmin(block, axis=0).tolist(),
            hi=np.nanmax(block, axis=0).tolist(),
            creep_um=float(med[0]), tau_s=_decay_tau(np.array(out['t']), med))
    return out


def fig_gain_dwells(dwell, out):
    """How the stage gets there: the dwells, as a median and an envelope.

    The time axis is linear because the claim is "flat by a second", and that
    is a claim about the end of the dwell rather than about the first 50 ms.
    """
    with plt.rc_context():
        fig, ax = plt.subplots(figsize=(vs.FULL, 2.6),
                               constrained_layout=True)
        t = np.array(dwell['t'])

        ax.axvspan(dwell['read_from_s'], dwell['span_s'], color=vs.INK_MUTED,
                   alpha=0.10, lw=0, zorder=0)
        ax.axhline(0, color=vs.INK_MUTED, lw=0.6, zorder=1)
        for d, colour, name in ((1, vs.FORWARD, '+ve'), (-1, vs.REVERSE, '−ve')):
            part = dwell['by_direction'][d]
            ax.fill_between(t, part['lo'], part['hi'], color=colour,
                            alpha=0.16, lw=0, zorder=2)
            ax.plot(t, part['median'], color=colour, lw=0.9, zorder=3,
                    label=f"approached {name} (n = {part['n']})")

        ax.grid(False)
        ax.set_xlim(0, dwell['span_s'])
        ax.set_xticks(np.arange(0, dwell['span_s'], 0.5))
        ax.set_ylim(-0.25, 0.25)
        ax.set_xlabel('time from the end of commanded motion [s]')
        ax.set_ylabel('position − settled value [µm]')
        vs.headline(fig, 'The stage creeps < 0.2 µm into position and '
                         'settles within ~1 s')
        ax.annotate('read window', (0.5 * (dwell['read_from_s']
                                           + dwell['span_s']), -0.235),
                    ha='center', va='bottom', fontsize=7,
                    color='#52514e')
        vs.despine(ax)

        fig.legend(loc='outside lower center', ncol=2, handletextpad=0.5,
                   columnspacing=1.6, borderpad=0.1)
        return vs.save(fig, out, 'stage_02b_gain_dwells')


def _electrical_axis(ax, mres):
    """Top axis in electrical degrees: 4 full steps = 4 x mres microsteps."""
    per_period = 4 * mres
    twin = ax.secondary_xaxis(
        'top', functions=(lambda i: i * 360.0 / per_period,
                          lambda deg: deg * per_period / 360.0))
    twin.set_xticks([0, 180, 360, 540, 720])
    twin.set_xlabel('electrical angle (deg)', fontsize=8)
    twin.tick_params(labelsize=7.5, colors=vs.INK_MUTED)
    for side in twin.spines.values():
        side.set_visible(False)


def fig_cumulative_error(res, out):
    """The error carried through the whole stepping sequence, out and back.

    Phase 2 drawn as one continuous run rather than as two passes against
    commanded position: 128 microsteps out, the reversal, 128 back. Plotted
    this way the lost motion is the step at the turnaround, the screw and
    microstep error is the wobble along each arm, and whether the leg came
    home is the value at the right-hand end.
    """
    fig, axes = plt.subplots(1, len(res), figsize=(vs.FULL, 2.9),
                             sharey=True,
                             constrained_layout=True, squeeze=False)
    closes, traces = [], []
    for col, leg in enumerate(res):
        ax = axes[0][col]
        n, mres = leg['n'], leg['mres']
        # Measured minus commanded, slot by slot, from the leg's own start
        # through the whole sequence in the order it ran - including the
        # turnaround slot, which the two passes exclude but the carriage went
        # through. This way round a carriage that is behind its command on
        # the way out reads negative, which is what "behind" ought to look
        # like, and it is the same convention as the folded map.
        measured = np.array([np.nan if v is None else v
                             for v in leg['position_um']])
        k = np.arange(measured.size)
        commanded = np.where(k <= n, k, 2 * n - k) * leg['step_nominal_um']
        error = measured - commanded

        ax.axvline(n, color=vs.INK_MUTED, lw=0.8, ls=':', zorder=1)
        ax.axhline(0, color=vs.INK_MUTED, lw=0.7, zorder=1)
        ax.plot(k[:n + 1], error[:n + 1], color=vs.FORWARD, lw=1.2, zorder=3,
                label='out')
        ax.plot(k[n:], error[n:], color=vs.REVERSE, lw=1.2, zorder=3,
                label='back')
        # The two values worth reading off: where the carriage stood against
        # its command when the direction changed, and where it stood when the
        # leg finished. Both labels hang below their point, which is the one
        # direction clear of the trace at both places - the reversal is the
        # bottom of the out arm, and the back arm arrives at the end from
        # above - so neither crosses a line.
        last = int(np.flatnonzero(np.isfinite(error))[-1])
        closes.append(error[last])
        for x, y, text, ha, dx in (
                (n, error[n], 'at the reversal', 'center', 0),
                (k[last], error[last], 'at the end', 'right', -5)):
            ax.plot(x, y, 'o', color=vs.INK, ms=4.5, mec=vs.SURFACE, mew=0.8,
                    zorder=5)
            # The end label hangs back from its point rather than sitting
            # over it: at the right-hand edge of the panel a centred one
            # would run off the axes.
            ax.annotate(f'{text}\n{y:+.2f} µm', (x, y),
                        textcoords='offset points', xytext=(dx, -13),
                        ha=ha, va='top', fontsize=7.5,
                        color=vs.INK_SOFT, linespacing=1.4,
                        # A label can still land on the zero reference where
                        # a leg turns round near it, so it carries the
                        # background with it.
                        bbox=dict(facecolor=vs.SURFACE, edgecolor='none',
                                  pad=1.4, alpha=0.9),
                        arrowprops=dict(arrowstyle='-', lw=0.6,
                                        color=vs.INK_MUTED, shrinkB=3))
        for full_step in range(0, 2 * n + 1, 4 * mres):
            ax.axvline(full_step, color=vs.GRID, lw=0.5, zorder=0)

        ax.grid(False)
        ax.set_xlim(0, 2 * n)
        ax.set_xticks(np.arange(0, 2 * n + 1, 64))
        ax.set_xlabel('microstep through the leg')
        vs.corner(ax, f"{leg['location_mm']:+.0f} mm", 'upper left')
        vs.panel(ax, 'abc'[col])
        if col == 0:
            ax.set_ylabel('measured − commanded, cumulative [µm]')
            # Low and left, where no leg goes: the location label has the
            # top corner.
            ax.legend(loc='lower left', ncol=2)
        vs.despine(ax)
        traces.append(error)

    # Headroom above the data for the labels, shared so the three panels stay
    # comparable.
    stack = np.concatenate(traces)
    lo, hi = np.nanmin(stack), np.nanmax(stack)
    axes[0][0].set_ylim(lo - 0.42 * (hi - lo), hi + 0.12 * (hi - lo))

    vs.headline(fig, f'Cumulative error: about 2 µm lost at each turnaround, '
                     f'and up to {max(abs(c) for c in closes):.1f} µm of it '
                     f'never recovered')
    return vs.save(fig, out, 'stage_03_cumulative_error')


def fig_microsteps(res, out):
    """Every commanded microstep as its own displacement, all 256 of them.

    One panel per location, each with the distribution of its own step sizes
    beside it. The microsteps that fail to clear half of one are marked and
    numbered, because where they are is the whole of what this figure has to
    say: they are the first microsteps after a direction reversal and
    nowhere else.
    """
    fig = plt.figure(figsize=(vs.FULL, 3.0), constrained_layout=True)
    gs = fig.add_gridspec(1, 3 * len(res), width_ratios=[6, 1.5, 0.7] * len(res))
    nominal = res[0]['step_nominal_um']
    half = nominal / 2
    axes = []

    for col, leg in enumerate(res):
        ax = fig.add_subplot(gs[0, 3 * col], sharey=axes[0] if axes else None)
        axes.append(ax)
        n = leg['n']
        step = np.array([np.nan if v is None else v for v in leg['step_all_um']])
        k = np.arange(1, step.size + 1)
        small = step <= half

        # The band is the leg's own mean and spread, not the nominal: the
        # question of whether the steps are the right size is the dashed
        # line, and the question of how alike they are is the band.
        ax.axhspan(leg['step_mean_um'] - leg['step_sd_all_um'],
                   leg['step_mean_um'] + leg['step_sd_all_um'],
                   color=vs.INK_MUTED, alpha=0.10, lw=0, zorder=0)
        ax.axhline(nominal, color=vs.INK, lw=0.9, ls='--', zorder=2)
        ax.axhline(half, color=vs.INK_MUTED, lw=0.7, ls=':', zorder=2)
        for name, keep, colour in (('out', k <= n, vs.FORWARD),
                                   ('back', k > n, vs.REVERSE)):
            ax.plot(k[keep], step[keep], 'o', ms=2.2, mec='none', color=colour,
                    alpha=0.85, zorder=3, label=name)
        # Two marks, because there are two different failures and one of
        # them is not a failure to move: a ring is a microstep that moved the
        # carriage less than half of one, a cross is a microstep that left no
        # step the signal can be seen to take. Most are both; microstep 133
        # at 0 mm is small but perfectly visible, and marking it as missing
        # would be wrong.
        gone = np.isin(k, leg['steps_undetected'])
        ax.plot(k[small], step[small], 'o', ms=7, mfc='none',
                mec=vs.STATUS_BAD, mew=1.1, zorder=4)
        ax.plot(k[gone], step[gone], 'x', ms=5.5, color=vs.STATUS_BAD,
                mew=1.3, zorder=5)
        for index in np.unique(np.r_[k[small], k[gone]]):
            # Labelled to the side rather than below: two of these land on
            # top of each other at -5 mm, where microsteps 1 and 2 are both
            # short.
            ax.annotate(f'{index}', (index, step[index - 1]),
                        textcoords='offset points', xytext=(10, 0),
                        ha='left', va='center', fontsize=7,
                        color=vs.STATUS_BAD)
        ax.axvline(n + 0.5, color=vs.INK_MUTED, lw=0.8, ls=':', zorder=1)

        ax.grid(False)
        ax.set_xlim(0, 2 * n + 1)
        ax.set_xticks([1, 64, 128, 192, 256])
        ax.set_ylim(0, 1.05)
        ax.set_xlabel('commanded microstep')
        vs.corner(ax, f"{leg['location_mm']:+.0f} mm", 'upper left')
        vs.panel(ax, 'abc'[col])
        if col == 0:
            ax.set_ylabel('measured displacement [µm]')
        else:
            ax.tick_params(labelleft=False)
        vs.despine(ax)

        # The distribution beside each panel: the bulk is tight and the
        # failures are the tail, which a time series alone understates.
        axh = fig.add_subplot(gs[0, 3 * col + 1], sharey=ax)
        bins = np.arange(0, 1.05, 0.025)
        axh.hist(step[k <= n], bins=bins, orientation='horizontal',
                 color=vs.FORWARD, alpha=0.55, lw=0)
        axh.hist(step[k > n], bins=bins, orientation='horizontal',
                 color=vs.REVERSE, alpha=0.55, lw=0)
        axh.axhline(nominal, color=vs.INK, lw=0.9, ls='--', zorder=2)
        axh.grid(False)
        axh.tick_params(labelleft=False, labelbottom=False, length=0)
        axh.set_xlabel('count', fontsize=8)
        vs.despine(axh, keep=())

    resolved = sum(leg['steps_resolved'] for leg in res)
    total = sum(leg['steps_measured'] for leg in res)
    gone_total = sum(len(leg['steps_undetected']) for leg in res)
    fig.legend(handles=[
        Line2D([], [], marker='o', ls='none', ms=4, color=vs.FORWARD,
               label='out'),
        Line2D([], [], marker='o', ls='none', ms=4, color=vs.REVERSE,
               label='back'),
        Line2D([], [], marker='o', ls='none', ms=7, mfc='none',
               mec=vs.STATUS_BAD, mew=1.1,
               label='moved less than half a microstep'),
        Line2D([], [], marker='x', ls='none', ms=5.5, color=vs.STATUS_BAD,
               mew=1.3, label='left no step the signal can be seen to take')],
        loc='outside lower left', ncol=4, fontsize=8.5, handletextpad=0.3,
        columnspacing=1.6)
    vs.headline(fig, f'{resolved} of {total} microsteps clear half a '
                     f'microstep and {total - gone_total} leave a visible '
                     f'step; the exceptions all follow a reversal')
    return vs.save(fig, out, 'stage_03b_step_size')


def fig_recovery(res, micro, out):
    """What the microsteps after a reversal fail to deliver, added up."""
    fig, axes = plt.subplots(1, len(res), figsize=(vs.FULL, 2.7),
                             sharey=True,
                             constrained_layout=True, squeeze=False)
    for col, (leg, m) in enumerate(zip(res, micro['per_leg'])):
        ax = axes[0][col]
        for name, colour, label in (('forward', vs.FORWARD, 'out'),
                                    ('reverse', vs.REVERSE, 'back')):
            part = m['reversal'][name]
            deficit = np.array(part['deficit_um'])
            ax.plot(np.arange(1, deficit.size + 1), deficit, color=colour,
                    lw=1.2, zorder=3, label=label)
            if part['plateau_um'] > 0.5:
                ax.hlines(part['plateau_um'], 40, 60, color=colour, lw=2.4,
                          alpha=0.6, zorder=4)
        ax.axhline(m['lost_motion_um'], color=vs.INK, lw=0.9, ls='--',
                   zorder=2)
        ax.axhline(0, color=vs.INK_MUTED, lw=0.7, zorder=1)
        ax.axvspan(40, 60, color=vs.INK_MUTED, alpha=0.08, lw=0, zorder=0)
        ax.grid(False)
        ax.set_ylim(-0.8, 3.2)
        ax.set_xlim(0, 127)
        ax.set_xlabel('microstep, counted from the reversal')
        vs.corner(ax, f"{leg['location_mm']:+.0f} mm", 'lower right')
        vs.panel(ax, 'abc'[col])
        if col == 0:
            ax.set_ylabel('cumulative shortfall D(n) [µm]')
            ax.legend(loc='lower left', ncol=2)
        vs.despine(ax)

    vs.headline(fig, 'Recovery after a reversal: half the lost motion in 6 '
                     'microsteps, the rest over tens')
    return vs.save(fig, out, 'stage_03d_recovery')


def fig_microstep_map(micro, out):
    """The microstep error folded onto one electrical period.

    The tilt and the lost-motion gap are both gone - the line is fitted out
    with the harmonics in one go, since over the two periods a leg covers a
    straight line and a once-per-period sine are not independent - so what is
    left is the shape the drive repeats every four full steps, averaged over
    every place it was measured.
    """
    folded = micro['folded']
    angle = np.array(folded['angle_deg'])
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(vs.FULL, 2.9), sharey=True,
                                  constrained_layout=True)

    mean = np.array(folded['mean_um'])
    sd = np.array(folded['sd_um'])
    ax.fill_between(angle, mean - sd, mean + sd, color=vs.SERIES[0],
                    alpha=0.18, lw=0, zorder=2, label='± 1 sd across passes')
    ax.plot(angle, mean, color=vs.SERIES[0], lw=1.4, zorder=3,
            label=f"mean of {int(np.max(folded['n_per_angle']))} passes")
    ax.plot(angle, folded['harmonic_curve_um'], color=vs.INK, lw=1.0, ls='--',
            zorder=4, label='harmonics 1, 2 and 4 per period')
    ax.axhline(0, color=vs.INK_MUTED, lw=0.7, zorder=1)
    for full_step in range(1, 4):
        ax.axvline(full_step * 90, color=vs.GRID, lw=0.6, zorder=0)
    ax.set_xlim(0, 360)
    ax.set_xticks(np.arange(0, 361, 90))
    ax.set_xlabel('electrical angle [deg], one period = 4 full steps')
    ax.set_ylabel('microstep error [µm]')
    vs.panel(ax, 'a')
    ax.legend(loc='upper right', fontsize=7.5)
    vs.despine(ax)

    for name, colour in (('forward', vs.FORWARD), ('reverse', vs.REVERSE)):
        part = folded[name]
        ax2.plot(angle, part['mean_um'], color=colour, lw=1.3)
        # Labelled where the two separate, which is the left-hand end; the
        # rms of each is a number for the caption.
        ax2.annotate(name, (angle[1], part['mean_um'][1]),
                     textcoords='offset points', xytext=(4, 0), ha='left',
                     va='center', fontsize=7.5, color=colour)
    ax2.axhline(0, color=vs.INK_MUTED, lw=0.7, zorder=1)
    for full_step in range(1, 4):
        ax2.axvline(full_step * 90, color=vs.GRID, lw=0.6, zorder=0)
    ax2.set_xlim(0, 360)
    ax2.set_xticks(np.arange(0, 361, 90))
    ax2.set_xlabel('electrical angle [deg]')
    vs.panel(ax2, 'b')
    vs.despine(ax2)

    vs.headline(fig, f"Microstep error over one electrical period: "
                     f"{folded['pv_um']:.2f} µm peak to peak, the same shape "
                     f"everywhere it is measured")
    return vs.save(fig, out, 'stage_03c_microstep_map')


SPEEDS = ('slow', 'moderate', 'fast')
SPEED_COLOUR = {'slow': '#9ec5f4', 'moderate': '#3987e5', 'fast': '#0d366b'}
HISTORY_COLOUR = {'continuation': '#86857f', 'reversal': '#eb6834'}
DIRECTION_COLOUR = {1: '#2a78d6', -1: '#eb6834'}


def _step(ax, x, y, n, **kwargs):
    """A numbered marker: ①..⑤ as a circled digit the body font will have."""
    return ax.annotate(str(n), (x, y), ha='center', va='center', fontsize=7,
                       zorder=6,
                       bbox=dict(boxstyle='circle,pad=0.22', fc='white',
                                 ec='black', lw=0.6), **kwargs)


def _anatomy_move(rows):
    """One move to open up: slow enough to see the ramp, typical in its error.

    Slow and 1 mm so the ramp lasts seconds rather than a tenth of one, and
    then the one whose settling error is nearest the median of the phase, so
    that the number the panel ends on is the number the next figure is about.
    """
    median = np.median([r['error_um'] for r in rows])
    slow = [r for r in rows if r['distance_mm'] == 1.0 and r['speed'] == 'slow']
    return min(slow or rows, key=lambda r: abs(r['error_um'] - median))


def fig_settling_sequence(rows, drift, q, model, out):
    """How a settling point is made: the phase, one move, and its error."""
    fig = plt.figure(figsize=(vs.FULL, 5.2), constrained_layout=True)
    gs = fig.add_gridspec(3, 1, height_ratios=(1.25, 0.95, 1.2))
    ax = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    ax3 = fig.add_subplot(gs[2], sharex=ax2)

    # ---- (a) the phase
    t0 = min(r['move_start'] for r in rows) / FS - 5
    t1 = max(r['dwell_end'] for r in rows) / FS + 5
    a, b = int(t0 * FS), int(t1 * FS)
    ax.plot((np.arange(a, b) / FS - t0) / 60, q[a:b] / 100.0,
            color='#bdbcb8', lw=0.4, zorder=2)
    ax.plot([(r['move_end'] / FS - t0) / 60 for r in rows],
            [r['settled_mm'] for r in rows], 'o', ms=1.8, mec='none',
            color='#52514e', zorder=3)
    inside = [d for d in drift if t0 < d['t'] < t1]
    ax.plot([(d['t'] - t0) / 60 for d in inside], [0] * len(inside), 'v',
            ms=4, color='black', zorder=4)
    ax.annotate('reference returns', ((inside[0]['t'] - t0) / 60 + 0.3, 0.4),
                fontsize=7, color='#52514e')

    move = _anatomy_move(rows)
    z0 = move['move_start'] / FS - 1.5
    # Stop at the end of the dwell: a step further and the next move drops
    # into the panel and takes the micrometre scale with it.
    z1 = min(move['dwell_end'], move['move_end'] + int(5 * FS)) / FS
    ax.axvspan((z0 - t0) / 60, (z1 - t0) / 60, color='black', alpha=0.12,
               lw=0, zorder=0)
    ax.annotate('(b), (c)', (((z0 + z1) / 2 - t0) / 60, 11.8), ha='center',
                va='bottom', fontsize=7)
    ax.set_xlim(0, (t1 - t0) / 60)
    ax.set_ylim(-11.5, 11.5)
    ax.set_yticks([-10, 0, 10])
    ax.set_ylabel('position [mm]')
    ax.set_xlabel('time within phase 3 [min]')
    vs.panel(ax, 'a')

    # ---- (b) one move, at the scale the move happens on
    za, zb = int(z0 * FS), int(z1 * FS)
    zt = np.arange(za, zb) / FS - move['move_end'] / FS
    ramp_start = move['move_start'] / FS - move['move_end'] / FS
    commanded = np.clip(
        move['from_mm'] + (zt - ramp_start) * move['speed_mm_s']
        * move['direction'],
        min(move['from_mm'], move['target_mm']),
        max(move['from_mm'], move['target_mm']))

    ax2.plot(zt, commanded, color='#86857f', lw=1.0, ls='--', zorder=2)
    ax2.plot(zt, q[za:zb] / 100.0, color='#2a78d6', lw=1.0, zorder=3)
    read_from = 0.6 * (z1 - move['move_end'] / FS)
    ax2.axvspan(ramp_start / 2, zt[-1], color='black', alpha=0.07, lw=0,
                zorder=0)
    ax2.annotate('zoomed below', (zt[-1] - 0.2, move['target_mm']),
                 textcoords='offset points', xytext=(0, 14), ha='right',
                 fontsize=7, color='#52514e')
    _step(ax2, ramp_start * 0.55, np.interp(ramp_start * 0.55, zt, commanded),
          1, textcoords='offset points', xytext=(-13, 7))
    ax2.annotate('move', (ramp_start * 0.55, np.interp(ramp_start * 0.55, zt,
                                                       commanded)),
                 textcoords='offset points', xytext=(-19, 7), ha='right',
                 va='center', fontsize=7.5)
    ax2.annotate('commanded', (ramp_start - 0.2, move['from_mm']),
                 ha='right', va='center', fontsize=7, color='#52514e')
    ax2.annotate('measured', (zt[-1] - 0.2, move['target_mm']),
                 textcoords='offset points', xytext=(0, -9), ha='right',
                 va='top', fontsize=7, color='#2a78d6')
    ax2.set_ylabel('position [mm]')
    ax2.tick_params(labelbottom=False)
    span = abs(move['commanded_mm'])
    ax2.set_ylim(min(move['from_mm'], move['target_mm']) - 0.45 * span,
                 max(move['from_mm'], move['target_mm']) + 0.45 * span)
    vs.panel(ax2, 'b')

    # ---- (c) the same move, at the scale the error happens on
    #
    # The constant the campaign takes off this move - its share of the gain,
    # the screw and the offset - is taken off the whole trace, not just the
    # settled value, so the number this panel ends on is the number the move
    # contributes to the next figure.
    # Evaluated at the commanded position at each instant, not once at the
    # target: a 1 mm move at this speed crosses half a screw revolution in
    # four seconds, so a constant screw term would leave that revolution in
    # the trace and it would read as a tracking error. At the stop the two
    # agree by construction, so the number the panel ends on is still the
    # number the move contributes to the next figure.
    angle = 2 * np.pi * (commanded % SCREW_MM) / SCREW_MM
    correction = (model['gain_um_per_mm'] * commanded + model['offset_um']
                  + model['screw_cos_um'] * np.cos(angle)
                  + model['screw_sin_um'] * np.sin(angle))
    error = (((q[za:zb] / 100.0 - commanded) * 1000.0 - correction)
             * move['direction'])
    ax3.axhline(0, color='black', lw=0.7, ls='--', zorder=2)
    # Above the line and over the ramp, where the trace is well below it.
    ax3.annotate('target', (ramp_start / 2, 0), textcoords='offset points',
                 xytext=(0, 3), fontsize=7, color='#52514e', va='bottom',
                 ha='center')
    ax3.axvspan(read_from, zt[-1], color='black', alpha=0.07, lw=0, zorder=0)
    ax3.plot(zt, error, color='#2a78d6', lw=0.8, zorder=3)

    settled = move['error_um']
    arrow_at = read_from + 0.55
    ax3.annotate('', (arrow_at, settled), (arrow_at, 0),
                 arrowprops=dict(arrowstyle='<->', lw=0.7, color='black'))
    # The label sits above the target line: the gap between that line and
    # the trace is a micrometre, which at this scale is not a line of type.
    _step(ax3, arrow_at, 0.0, 5, textcoords='offset points', xytext=(0, 13))
    ax3.annotate(f'settling error {settled:+.2f} µm', (arrow_at - 0.2, 0),
                 textcoords='offset points', xytext=(0, 10), ha='right',
                 va='center', fontsize=7.5)
    # At the far end of the band, so it clears the arrow and its label.
    _step(ax3, zt[-1] - 0.12, 0.93, 4, xycoords=('data', 'axes fraction'))
    ax3.annotate('read window', (zt[-1] - 0.3, 0.93),
                 xycoords=('data', 'axes fraction'),
                 textcoords='offset points', xytext=(0, 0), ha='right',
                 va='center', fontsize=7.5)
    _step(ax3, 1.0, np.interp(1.0, zt, error), 3,
          textcoords='offset points', xytext=(0, -11))
    ax3.annotate('creep', (1.0, np.interp(1.0, zt, error)),
                 textcoords='offset points', xytext=(0, -18), ha='center',
                 va='top', fontsize=7.5)
    ax3.annotate('tracking error while moving',
                 (ramp_start * 0.5, np.interp(ramp_start * 0.5, zt, error)),
                 textcoords='offset points', xytext=(0, -8), ha='center',
                 va='top', fontsize=7, color='#52514e')

    low, high = float(np.min(error)), float(np.max(error))
    pad = 0.25 * (high - low)
    ax3.set_ylim(low - pad, high + pad)
    ax3.set_xlim(zt[0], zt[-1])
    ax3.set_ylabel('error, signed along travel [µm]\n(negative = short)')
    ax3.set_xlabel('time from the end of commanded motion [s]')
    vs.panel(ax3, 'c')

    # ---- the instant both panels share
    for axis, y in ((ax2, 0.5), (ax3, 0.5)):
        axis.axvline(0, color='black', lw=0.6, ls=':', zorder=1)
    _step(ax3, 0, 0.08, 2, xycoords=('data', 'axes fraction'))
    ax3.annotate('t = 0, end of commanded motion', (0, 0.08),
                 xycoords=('data', 'axes fraction'),
                 textcoords='offset points', xytext=(9, 0), ha='left',
                 va='center', fontsize=7.5)

    vs.headline(fig, '90 moves: 3 distances × 3 speeds × 2 directions × 5 '
                     'repeats, in shuffled order')
    return vs.save(fig, out, 'stage_04_sequence')


def fig_settling_error(rows, delta_um, out):
    """Where the stage stops, how far it goes, and why those differ.

    The two measurements disagree about the reversal and the sketch is there
    to say why they are both right: the carriage always comes to rest the
    same shortfall behind its target, so where it stops knows nothing about
    history, while a move that reverses starts from the far side of that
    shortfall and so travels twice it less than it was told to.
    """
    fig = plt.figure(figsize=(vs.FULL, 4.7), constrained_layout=True)
    gs = fig.add_gridspec(2, 2, height_ratios=(1.9, 1), width_ratios=(1.7, 1))
    ax = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, :])
    distances = sorted({r['distance_mm'] for r in rows})
    rng = np.random.default_rng(4)

    # ---- (a) where it stops
    for g, distance in enumerate(distances):
        for s, speed in enumerate(SPEEDS):
            x = g * 4 + s
            cell = [r['error_um'] for r in rows
                    if r['distance_mm'] == distance and r['speed'] == speed]
            ax.plot(x + rng.uniform(-0.3, 0.3, len(cell)), cell, 'o', ms=2.4,
                    mec='none', color=SPEED_COLOUR[speed], alpha=0.9, zorder=3)
            ax.hlines(np.mean(cell), x - 0.42, x + 0.42, color='black',
                      lw=1.0, zorder=4)
            ax.annotate(speed[0], (x, -0.03),
                        xycoords=('data', 'axes fraction'), ha='center',
                        va='top', fontsize=7, color='#52514e')
        ax.annotate(f'{distance:g} mm', (g * 4 + 1, -0.13),
                    xycoords=('data', 'axes fraction'), ha='center', va='top',
                    fontsize=8)

    ax.axhline(0, color='#86857f', lw=0.6, zorder=1)
    ax.axhline(-delta_um, color='black', lw=0.7, ls='--', zorder=2)
    # In the margin the data leaves free at the right, not over the points.
    ax.annotate('−δ, steady-state\nshortfall (phase 4)',
                (len(distances) * 4 - 0.6, -delta_um), fontsize=7,
                ha='left', va='center', color='#52514e', linespacing=1.4)
    ax.set_xlim(-1, len(distances) * 4 + 1.9)
    ax.set_xticks([])
    ax.set_ylim(-2.3, 0.9)
    ax.set_ylabel('settling error, signed along travel [µm]\n'
                  '(negative = short)')
    vs.panel(ax, 'a')

    # ---- (b) how far it goes
    for h, (name, colour) in enumerate(HISTORY_COLOUR.items()):
        cell = [r['move_error_um'] for r in rows if r['history'] == name]
        ax2.plot(h + rng.uniform(-0.22, 0.22, len(cell)), cell, 'o', ms=2.4,
                 mec='none', color=colour, alpha=0.9, zorder=3)
        ax2.hlines(np.mean(cell), h - 0.34, h + 0.34, color='black', lw=1.0,
                   zorder=4)
        ax2.annotate(f'{name}\nn = {len(cell)}', (h, -0.03),
                     xycoords=('data', 'axes fraction'), ha='center',
                     va='top', fontsize=7.5, color=colour, linespacing=1.4)

    ax2.axhline(0, color='black', lw=0.7, ls='--', zorder=2)
    ax2.axhline(-2 * delta_um, color='black', lw=0.7, ls='--', zorder=2)
    ax2.annotate('0', (0.62, 0), fontsize=7, ha='left', va='bottom',
                 color='#52514e')
    ax2.annotate('−2δ', (0.62, -2 * delta_um), fontsize=7, ha='left',
                 va='bottom', color='#52514e')
    ax2.set_xlim(-0.6, 1.6)
    ax2.set_xticks([])
    span = [r['move_error_um'] for r in rows]
    ax2.set_ylim(min(span) - 0.3, max(span) + 0.6)
    ax2.set_ylabel('move error, signed along travel [µm]\n'
                   '(negative = travelled short)')
    vs.panel(ax2, 'b')

    # ---- (c) why both are true
    d = 1.0                      # delta, drawn far larger than it is
    a_pos, b_pos = 0.0, 6.0
    ax3.axhline(0, color='#86857f', lw=0.6, zorder=1)
    # Both tick labels are pushed away from the middle of the sketch, where
    # the two starting dots and their labels are.
    for x, label, ha, dx in ((a_pos, 'previous target', 'right', -0.12),
                             (b_pos, 'target', 'left', 0.12)):
        ax3.plot([x, x], [-0.45, 0.45], color='#86857f', lw=0.6, ls=':',
                 zorder=1)
        ax3.annotate(label, (x + dx, 0.5), ha=ha, va='bottom', fontsize=7,
                     color='#52514e')

    ax3.plot(a_pos - d, 0, 'o', ms=4, color=HISTORY_COLOUR['continuation'],
             zorder=4)
    ax3.plot(a_pos + d, 0, 'o', ms=4, color=HISTORY_COLOUR['reversal'],
             zorder=4)
    ax3.plot(b_pos - d, 0, 'o', ms=4, color='black', zorder=4)
    ax3.annotate('stopped here\nafter a +ve move', (a_pos - d, -0.2),
                 ha='center', va='top', fontsize=7,
                 color=HISTORY_COLOUR['continuation'], linespacing=1.4)
    ax3.annotate('stopped here\nafter a −ve move', (a_pos + d, 0.25),
                 ha='center', va='bottom', fontsize=7,
                 color=HISTORY_COLOUR['reversal'], linespacing=1.4)
    ax3.annotate('stops δ short,\nwhichever it was', (b_pos - d, 0.25),
                 ha='center', va='bottom', fontsize=7, linespacing=1.4)

    for y, start, colour, text in (
            (-0.75, a_pos - d, HISTORY_COLOUR['continuation'],
             'continuing: travels exactly what it was told to'),
            (-1.15, a_pos + d, HISTORY_COLOUR['reversal'],
             'reversing: travels 2δ less')):
        ax3.annotate('', (b_pos - d, y), (start, y),
                     arrowprops=dict(arrowstyle='->', lw=0.8, color=colour))
        ax3.annotate(text, ((start + b_pos - d) / 2, y - 0.08), ha='center',
                     va='top', fontsize=7, color=colour)

    ax3.annotate('', (b_pos, 0.95), (a_pos, 0.95),
                 arrowprops=dict(arrowstyle='<->', lw=0.7, color='black'))
    ax3.annotate('commanded distance', ((a_pos + b_pos) / 2, 1.0),
                 ha='center', va='bottom', fontsize=7)
    ax3.annotate('2δ', (a_pos, -0.32), ha='center', va='top', fontsize=7)
    ax3.annotate('', (a_pos + d, -0.22), (a_pos - d, -0.22),
                 arrowprops=dict(arrowstyle='<->', lw=0.7, color='black'))

    ax3.set_xlim(a_pos - 3.0, b_pos + 1.6)
    ax3.set_ylim(-1.75, 1.5)
    ax3.axis('off')
    ax3.annotate('position →', (b_pos + 1.1, 0), ha='right', va='bottom',
                 fontsize=7, color='#52514e')
    ax3.annotate('δ drawn far larger than it is', (a_pos - 2.3, -1.6),
                 ha='left', va='bottom', fontsize=7, color='#86857f')
    vs.panel(ax3, 'c', dx=0.0, dy=0.88)

    vs.headline(fig, 'Every stop lands δ short; only the travel knows about '
                     'the reversal')
    return vs.save(fig, out, 'stage_04a_settling_error')


def fig_settling_traces(traces, out):
    """The last five seconds of every stop, by direction of travel."""
    fig, ax = plt.subplots(figsize=(vs.FULL, 2.7), constrained_layout=True)
    ax.axhline(0, color='#86857f', lw=0.6, zorder=1)

    # Each trace is drawn to its own length; the medians are taken over a
    # padded array so that a short dwell shortens the sample behind the
    # median rather than cutting every trace back to the shortest one.
    longest = max(len(p['trace']) for p in traces)
    stack = {1: np.full((0, longest), np.nan), -1: np.full((0, longest), np.nan)}
    for direction in stack:
        block = [p['trace'] for p in traces if p['direction'] == direction]
        padded = np.full((len(block), longest), np.nan)
        for n, trace in enumerate(block):
            padded[n, :len(trace)] = trace
        stack[direction] = padded

    t = (np.arange(longest) + 1) / FS
    for part in traces:
        ax.plot(t[:len(part['trace'])], part['trace'],
                color=DIRECTION_COLOUR[part['direction']], lw=0.3, alpha=0.22,
                zorder=2)

    half = {}
    for direction, colour in DIRECTION_COLOUR.items():
        median = np.nanmedian(stack[direction], axis=0)
        ax.plot(t, median, color=colour, lw=1.3, zorder=4)
        crossing = np.flatnonzero(np.abs(median) <= abs(median[0]) / 2)
        half[direction] = t[crossing[0]] if crossing.size else np.nan
        ax.annotate('moving ' + ('+ve' if direction > 0 else '−ve'),
                    (0.0012, median[2]), textcoords='offset points',
                    xytext=(0, -9 if direction > 0 else 9), fontsize=7.5,
                    color=colour, ha='left',
                    va='top' if direction > 0 else 'bottom')

    reach = float(np.nanmax(list(half.values())))
    ax.axvline(reach, color='black', lw=0.6, ls=':', zorder=3)
    ax.annotate(f'half the creep gone by {reach:.1f} s', (reach, 0.385),
                textcoords='offset points', xytext=(-4, 0), fontsize=7,
                ha='right', va='top')

    ax.set_xscale('log')
    ax.set_xlim(0.001, SETTLE_WINDOW_S)
    ax.set_ylim(-0.4, 0.4)
    ax.set_xlabel('time from the end of commanded motion [s]')
    ax.set_ylabel('position − settled value [µm]')
    vs.headline(fig, f'Every stop creeps a fraction of a micrometre into '
                     f'place, half of it gone by {reach:.1f} s')
    return vs.save(fig, out, 'stage_04b_settling_traces')


APPROACH_COLOUR = {'forward': '#2a78d6', 'reverse': '#eb6834'}


def _both_ways(rep):
    """One target visited from both sides in the same cycle, to open up."""
    for cycle in (3, 2, 4, 1, 5):
        pair = [r for r in rep['moves']
                if r['target_mm'] == 0.0 and r['cycle'] == cycle]
        if len(pair) == 2:
            return sorted(pair, key=lambda r: r['move_start'])
    return rep['moves'][:2]


def fig_repeat_phase(rep, q, out):
    """Phase 4 end to end: the zig-zag, its stops, and the pair drawn in detail.

    The reference returns are deliberately not marked. They are a property of
    the campaign rather than of this phase - every phase takes them, and what
    they establish, that nothing drifts, is reported as a number instead.
    """
    fig, ax = plt.subplots(figsize=(vs.FULL, 1.7), constrained_layout=True)
    moves = rep['moves']
    t0 = min(r['move_start'] for r in moves) / FS - 8
    t1 = max(r['dwell_end'] for r in moves) / FS + 5
    a, b = int(t0 * FS), int(t1 * FS)

    pair = _both_ways(rep)
    z0 = pair[0]['run_up']['move_end'] / FS - 1.5
    z1 = pair[1]['move_end'] / FS + 2.0
    ax.axvspan((z0 - t0) / 60, (z1 - t0) / 60, color='black', alpha=0.12,
               lw=0, zorder=0)
    ax.annotate('one target, both approaches',
                (((z0 + z1) / 2 - t0) / 60, 11.3), ha='center', va='bottom',
                fontsize=7, color=vs.INK_SOFT)

    ax.plot((np.arange(a, b) / FS - t0) / 60, q[a:b] / 100.0,
            color='#bdbcb8', lw=0.4, zorder=2)
    ax.plot([(r['move_end'] / FS - t0) / 60 for r in moves],
            [r['settled_mm'] for r in moves], 'o', ms=1.8, mec='none',
            color='#52514e', zorder=3)

    ax.set_xlim(0, (t1 - t0) / 60)
    ax.set_ylim(-11.5, 13.8)
    ax.set_yticks([-10, 0, 10])
    ax.set_ylabel('position [mm]')
    ax.set_xlabel('time within phase 4 [min]')
    vs.despine(ax)
    vs.headline(fig, 'Five targets, each reached once from each side per '
                     'cycle, five cycles over')
    return vs.save(fig, out, 'stage_05_phase')


def fig_repeat_sequence(rep, q, out):
    """How a reversal value is measured: the approach, then the two stops.

    The phase itself - five targets, both directions, five cycles - is a
    sentence, and the plan it came from is drawn in the sequence figure, so
    neither is repeated here. What cannot be inferred from either is the
    run-up that establishes the approach direction and the window the settled
    value is read from, and those are the two panels.
    """
    # Extra width and inter-panel padding keep (b) clear of its top y tick.
    fig, (ax2, ax3) = plt.subplots(1, 2,
                                   figsize=(vs.FULL * 1.045, 2.3),
                                   constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.08, wspace=0.08)

    pair = _both_ways(rep)
    z0 = pair[0]['run_up']['move_end'] / FS - 1.5
    z1 = pair[1]['move_end'] / FS + 2.0

    # ---- (a) the same target from both sides, at millimetre scale
    za, zb = int(z0 * FS), int(z1 * FS)
    origin = pair[0]['move_end'] / FS
    zt = np.arange(za, zb) / FS - origin
    # The commanded profile across both approaches: each move ramps between
    # the positions the plan gives it, over the time the move took.
    target = pair[0]['target_mm']
    marks, levels = [zt[0]], [pair[0]['run_up']['from_mm']]
    for move in pair:
        for leg, destination in ((move['run_up'], move['run_up']['to_mm']),
                                 (move, move['target_mm'])):
            marks += [leg['move_start'] / FS - origin,
                      leg['move_end'] / FS - origin]
            levels += [levels[-1], destination]
    marks.append(zt[-1])
    levels.append(levels[-1])
    ax2.plot(zt, q[za:zb] / 100.0, color='#2a78d6', lw=1.2, zorder=3)
    ax2.plot(marks, levels, color='#52514e', lw=0.9, zorder=4,
             dashes=(3.5, 2.5))
    ax2.axhline(target, color='#86857f', lw=0.5, ls=':', zorder=1)
    for move in pair:
        ax2.plot(move['move_end'] / FS - origin, move['settled_mm'], 'o',
                 ms=4, mec='none', color=APPROACH_COLOUR[move['direction']],
                 zorder=4)
    # Both approaches are built the same way; the marks go on the second,
    # where the overshoot is above the target and there is room to label it.
    second = pair[1]
    hook = second['run_up']['to_mm']
    # The overshoot is held between the end of the run-up and the start of the
    # approach: the one place in the move where its size can be drawn.
    held = (second['run_up']['move_end']
            + second['move_start']) / 2 / FS - origin
    ax2.annotate('', (held, hook), (held, target),
                 arrowprops=dict(arrowstyle='<->', lw=0.7, color='black'))
    _step(ax2, held, hook, 1, textcoords='offset points', xytext=(0, 11))
    ax2.annotate('run-up, 0.5 mm past', (held, hook),
                 textcoords='offset points', xytext=(0, 18), ha='center',
                 va='bottom', fontsize=7.5)
    ramp = (second['move_start'] + second['move_end']) / 2 / FS - origin
    _step(ax2, ramp, (hook + target) / 2, 2,
          textcoords='offset points', xytext=(13, 1))
    ax2.annotate('approach', (ramp, (hook + target) / 2),
                 textcoords='offset points', xytext=(19, 1), ha='left',
                 va='center', fontsize=7.5)
    # Both labels on the first overshoot plateau, where the dashed command
    # lies over the measured, the two can be told apart, and nothing else is.
    plateau = (pair[0]['run_up']['move_end']
               + pair[0]['move_start']) / 2 / FS - origin
    ax2.annotate('measured', (plateau, pair[0]['run_up']['to_mm']),
                 textcoords='offset points', xytext=(0, 5), ha='center',
                 va='bottom', fontsize=7, color='#2a78d6')
    ax2.annotate('commanded', (plateau, pair[0]['run_up']['to_mm']),
                 textcoords='offset points', xytext=(0, -5), ha='center',
                 va='top', fontsize=7, color='#52514e')
    ax2.set_xlim(zt[0], zt[-1])
    span = abs(hook - target)
    ax2.set_ylim(target - 1.45 * span, target + 1.45 * span)
    ax2.set_xlabel('time from the first stop [s]')
    ax2.set_ylabel('position [mm]')
    # This is a sequence schematic: keep the physical units but suppress the
    # particular numerical tick values of the selected example.
    ax2.tick_params(axis='both', which='both', labelbottom=False,
                    labelleft=False)
    vs.panel(ax2, 'a')

    # ---- (b) the two holds, at micrometre scale, which is the same data
    # as the two coloured stops of (a) magnified a thousand times
    settled, shown, last = {}, [], 0.0
    for move in pair:
        end = move['move_end']
        stop = min(move['dwell_end'], end + int(3.0 * FS))
        t = (np.arange(end, stop) - end) / FS
        deviation = (q[end:stop] / 100.0 - move['target_mm']) * 1000
        shown.append(deviation)
        last = max(last, float(t[-1]))
        ax3.plot(t, deviation, color=APPROACH_COLOUR[move['direction']],
                 lw=0.95, zorder=4)
        settled[move['direction']] = move['deviation_um'] + rep['offset_removed_um']
        ax3.annotate('approached from ' + ('below' if move['direction']
                                           == 'forward' else 'above'),
                     (t[0], deviation[:200].mean()),
                     textcoords='offset points', xytext=(2, 7), ha='left',
                     va='bottom', fontsize=7.5,
                     color=APPROACH_COLOUR[move['direction']])

    ax3.axvspan(READ_FROM_S, READ_TO_S, color='black', alpha=0.07, lw=0,
                zorder=0)
    read_mid = (READ_FROM_S + READ_TO_S) / 2
    _step(ax3, READ_FROM_S + 0.12, 0.94, 3,
          xycoords=('data', 'axes fraction'))
    ax3.annotate('read window', (READ_FROM_S + 0.12, 0.94),
                 xycoords=('data', 'axes fraction'),
                 textcoords='offset points', xytext=(9, 0), ha='left',
                 va='center', fontsize=7.5)

    # Draw B inside the window from which both settled levels are read, and
    # label it back along the gap, which is now most of the panel.
    ax3.annotate('', (read_mid, settled['forward']),
                 (read_mid, settled['reverse']),
                 arrowprops=dict(arrowstyle='<->', lw=0.9, color='black'),
                 zorder=6)
    middle = float(np.mean(list(settled.values())))
    _step(ax3, read_mid, middle, 4, textcoords='offset points',
          xytext=(-12, 0))
    ax3.annotate('reversal value B', (read_mid, middle),
                 textcoords='offset points', xytext=(-18, 0), ha='right',
                 va='center', fontsize=7.5)

    # The pair is all that is drawn, so the axes are its own extent: the two
    # levels and the gap between them fill the panel, with headroom only
    # where the read-window mark needs it.
    ax3.set_xlim(0, last)
    low = min(float(np.nanmin(v)) for v in shown)
    high = max(float(np.nanmax(v)) for v in shown)
    ax3.set_ylim(low - 0.20 * (high - low), high + 0.38 * (high - low))
    ax3.set_xlabel('time from the end of commanded motion [s]')
    ax3.set_ylabel('deviation from target [µm]')
    vs.corner(ax3, 'the two stops of (a)', 'lower right',
              fontsize=7, color=vs.INK_MUTED)
    vs.panel(ax3, 'b', dx=-0.045)

    vs.headline(fig, '5 targets × 2 directions × 5 cycles, each approached '
                     'over a 0.5 mm run-up')
    return vs.save(fig, out, 'stage_05_sequence')


def fig_repeatability(rep, out):
    """The ISO 230-2 result, target by target, with the numbers beside it."""
    fig = plt.figure(figsize=(vs.FULL, 2.5), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=(1.35, 1.0))
    ax = fig.add_subplot(gs[0])
    axt = fig.add_subplot(gs[1])
    per = rep['per_target']
    x = np.array([p['target_mm'] for p in per])

    # One dashed envelope per target, enclosing both approaches: it says
    # which marks belong to the same commanded position, which is otherwise
    # only implied by them being near each other.
    for p in per:
        values = np.array(p['forward'] + p['reverse'])
        lo = min(values.min(), p['mean_forward'] - 2 * p['s_forward'],
                 p['mean_reverse'] - 2 * p['s_reverse'])
        hi = max(values.max(), p['mean_forward'] + 2 * p['s_forward'],
                 p['mean_reverse'] + 2 * p['s_reverse'])
        pad = 0.18
        ax.add_patch(Rectangle(
            (p['target_mm'] - 0.72, lo - pad), 1.44, (hi - lo) + 2 * pad,
            fill=False, edgecolor=vs.INK_MUTED, lw=0.6, ls=(0, (3, 2)),
            zorder=1))

    for direction, marker, dx in (('forward', '^', -0.3),
                                  ('reverse', 'v', +0.3)):
        colour = APPROACH_COLOUR[direction]
        mean = np.array([p[f'mean_{direction}'] for p in per])
        sd = np.array([p[f's_{direction}'] for p in per])
        for p, offset in zip(per, x + dx):
            ax.plot(np.full(len(p[direction]), offset), p[direction], 'o',
                    ms=2, mec='none', color=colour, alpha=0.3, zorder=2)
        ax.errorbar(x + dx, mean, yerr=2 * sd, fmt=marker, ms=5,
                    color=colour, mec='none', elinewidth=0.9, capsize=2.5,
                    capthick=0.9, zorder=4)

    span = np.array([x.min() - 1.4, x.max() + 1.4])
    ax.plot(span, rep['gain_um_per_mm'] * span + rep['gain_intercept_um'],
            color='black', lw=1.2, dashes=(5, 2.5), zorder=5)
    gain_x = 2.6
    ax.annotate(f"gain {rep['gain_pct']:+.3f} %",
                (gain_x, rep['gain_um_per_mm'] * gain_x
                 + rep['gain_intercept_um']),
                textcoords='offset points', xytext=(0, 9), ha='center',
                va='bottom', fontsize=7)

    ax.axhline(0, color='#86857f', lw=0.6, zorder=1)
    ax.set_xticks(x)
    ax.set_xlim(x.min() - 1.25, x.max() + 1.25)
    ax.set_xlabel('target position [mm]')
    ax.set_ylabel('deviation from target [µm]\n(phase mean removed)')
    ax.annotate('▲ approached from below', (0.015, 0.96),
                xycoords='axes fraction', fontsize=7.5, va='top',
                color=APPROACH_COLOUR['forward'])
    ax.annotate('▼ approached from above', (0.015, 0.88),
                xycoords='axes fraction', fontsize=7.5, va='top',
                color=APPROACH_COLOUR['reverse'])
    vs.panel(ax, 'a')

    # ---- (b) the same five targets as numbers
    axt.set_xlim(0, 1)
    axt.set_ylim(0, 1)
    axt.set_axis_off()
    head = ('x [mm]', 'x̄↑', 'x̄↓', 'min', 'max', 'mean')
    body = []
    for p in per:
        values = np.array(p['forward'] + p['reverse'])
        body.append(tuple(
            cell.replace('-', '\u2212') for cell in (
                f"{p['target_mm']:+.0f}",
                f"{p['mean_forward']:+.2f}",
                f"{p['mean_reverse']:+.2f}",
                f"{values.min():+.2f}",
                f"{values.max():+.2f}",
                f"{values.mean():+.2f}")))

    right = (0.13, 0.31, 0.48, 0.66, 0.83, 1.0)
    top, step = 0.78, 0.105
    axt.annotate('deviation from target [µm], 10 visits per position',
                 (0, 0.95), fontsize=7, color=vs.INK_SOFT, va='bottom')
    for line in (top + 0.055, top - 0.045):
        axt.plot([0, 1], [line, line], color=vs.INK_SOFT, lw=0.7,
                 clip_on=False)
    for column, label in zip(right, head):
        axt.annotate(label, (column, top), ha='right', va='center',
                     fontsize=7.5, fontweight='bold')
    for row, cells in enumerate(body, start=1):
        y = top - row * step
        for column, cell in zip(right, cells):
            axt.annotate(cell, (column, y), ha='right', va='center',
                         fontsize=7.5)
    axt.plot([0, 1], [top - len(body) * step - 0.05] * 2, color=vs.INK_SOFT,
             lw=0.7, clip_on=False)
    vs.panel(axt, 'b', dx=0.0)

    vs.headline(fig, f"From either direction the stage returns to each "
                     f"target within {rep['repeatability_forward']:.2f} µm")
    return vs.save(fig, out, 'stage_05a_repeatability')


def fig_campaign(summary, out):
    """The two quantities every phase measured, gathered in one view."""
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(vs.FULL, 2.4),
                                  constrained_layout=True)

    for axis, rows, unit, label in (
            (ax, summary['lost_motion'], 'µm',
             'lost motion between directions [µm]'),
            (ax2, summary['gain'], '%', 'gain [% of nominal]')):
        y = np.arange(len(rows))[::-1]
        value = np.array([r['value'] for r in rows])
        error = np.array([r['error'] for r in rows])
        finite = np.isfinite(error)
        mean, sd = rows[0]['pooled_mean'], rows[0]['pooled_sd']
        axis.axvspan(mean - sd, mean + sd, color='black', alpha=0.07, lw=0,
                     zorder=0)
        axis.axvline(mean, color='#52514e', lw=0.7, zorder=1)
        axis.annotate(f"{mean:.2f} ± {sd:.2f} {unit}" if unit == 'µm'
                      else f"{mean:+.3f} ± {sd:.3f} {unit}",
                      (mean, 1.0), xycoords=('data', 'axes fraction'),
                      textcoords='offset points', xytext=(0, 2), ha='center',
                      va='bottom', fontsize=7, color='#52514e')
        axis.errorbar(value, y, xerr=np.where(finite, error, 0), fmt='o',
                      ms=4, color='#2a78d6', mec='none', elinewidth=0.9,
                      capsize=2.5, capthick=0.9, zorder=3)
        axis.set_yticks(y)
        axis.set_yticklabels([f"phase {r['phase']}\n{r['what']}"
                              for r in rows], fontsize=7.5)
        axis.set_ylim(-0.7, len(rows) - 0.25)
        axis.set_xlabel(label)
        axis.tick_params(axis='y', length=0)

    ax2.axvline(0, color='#86857f', lw=0.6, zorder=1)
    vs.panel(ax, 'a')
    vs.panel(ax2, 'b')
    lost, gain = summary['lost_motion'][0], summary['gain'][1]
    vs.headline(fig, f"Four phases, two quantities: lost motion "
                     f"{lost['pooled_mean']:.2f} ± {lost['pooled_sd']:.2f} µm "
                     f"and gain {gain['pooled_mean']:+.4f} ± "
                     f"{gain['pooled_sd']:.4f} %")
    return vs.save(fig, out, 'stage_05b_campaign')


def fig_noise(res, counts, align_info, out):
    fig, (ax, ax2) = plt.subplots(1, 2, figsize=(vs.FULL, 2.6),
                                  constrained_layout=True)
    first = align_info['seg_s'][0]
    rest = counts[:first]
    t = np.arange(rest.size) / FS
    ax.plot(t, (rest - np.median(rest)) / 1000.0, color=vs.SERIES[0], lw=0.5)
    ax.set_xlabel('time within the recording [s]')
    ax.set_ylabel('position [µm]')
    vs.panel(ax, 'a')
    vs.despine(ax)

    freq = np.array(res['psd_freq'])
    psd = np.array(res['psd'])
    keep = (freq > 1.0) & (freq < 500)
    ax2.semilogy(freq[keep], psd[keep] / psd[keep].max(), color=vs.SERIES[0],
                 lw=0.7)
    for f in res['dominant_hz'][:3]:
        ax2.annotate(f'{f:g} Hz', (f, 1.0), textcoords='offset points',
                     xytext=(0, -2), ha='center', va='top', fontsize=7,
                     color=vs.INK_SOFT)
    ax2.set_xlabel('frequency [Hz]')
    ax2.set_ylabel('relative power')
    vs.panel(ax2, 'b')
    vs.despine(ax2)
    vs.headline(fig, 'The encoder at standstill, energised, before the run')
    return vs.save(fig, out, 'stage_06_standstill')


# --------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--cache', type=Path, default=HERE / 'cache')
    ap.add_argument('--plan', type=Path,
                    default=STAGE_DIR / 'sequence' / 'sequence_plan.csv')
    ap.add_argument('--out', type=Path, default=HERE / 'plots')
    args = ap.parse_args()

    vs.use_thesis()
    data = np.load(args.cache / 'IDdata.npz')
    counts = data['counter'].astype(np.float64)
    plan = load_plan(args.plan)
    print(f'{counts.size} samples ({counts.size / FS:.1f} s), {len(plan)} planned moves')

    info = build_alignment(counts, plan)
    gain, rms, worst = fit_scale(counts, info['seg_s'], info['seg_e'], plan,
                                 info['shift'])
    matched = len(info['pairs'])
    print(f'aligned {matched}/{len(plan)} moves; {len(info["skipped"])} planned '
          f'moves left no segment, {len(info["extra"])} segments match no move')
    print(f'scale {gain:,.1f} counts/full step (fit rms {rms:.3f} full steps)')

    row_end = {i: info['seg_e'][j] for i, j in info['pairs']}
    offset = np.median([counts[row_end[i]] - gain * plan[i]['target']
                        for i in sorted(row_end)[:20]])
    q = (counts - offset) / gain

    records = settled_positions(q, info, plan)
    legs = staircase_legs(q, counts, info, plan)
    results = dict(
        recording=dict(samples=int(counts.size), seconds=float(counts.size / FS),
                       sample_hz=FS, channel=str(data['channel']),
                       start=str(data['starttime']), end=str(data['endtime'])),
        alignment=dict(planned=len(plan), matched=matched,
                       skipped=[dict(row=int(i), phase=plan[i]['phase'],
                                     label=plan[i]['label']) for i in info['skipped']],
                       extra=[dict(t=float(info['seg_s'][j] / FS),
                                   displacement_full_steps=float(
                                       (counts[info['seg_e'][j]]
                                        - counts[info['seg_s'][j]]) / gain))
                              for j in info['extra']],
                       scale_fit_rms_full_steps=rms,
                       scale_fit_max_full_steps=worst),
        gain=gain_result(records, plan, gain),
        settling=[],
        repeatability=repeatability_result(records, plan, q),
        staircase=staircase_result(legs),
        standstill=noise_result(counts, info),
        nominal=dict(counts_per_rev=NOMINAL_COUNTS_PER_REV,
                     full_step_um=FULL_STEP_UM))

    drift = sorted(i for i in records if plan[i]['phase'] == 'drift')
    results['gain_pair'] = gain_pair_estimate(results['gain'])
    results['settling'], results['settling_model'] = settling_result(
        records, plan, q, repeatability_gain(results['repeatability']))
    results['gain_dwells'] = gain_dwell_result(results['gain'], records, q)
    results['gain_moves'] = gain_moves_result(results['gain'], records, plan, q)

    results['drift'] = [dict(t=records[i]['move_end'] / FS,
                             error_um=(records[i]['settled'] - records[i]['target'])
                             * FULL_STEP_UM,
                             # Which way the stage came into the reference.
                             # Two checks approached from opposite sides
                             # differ by the lost motion before they differ
                             # by anything that happened in between.
                             direction=1 if plan[i]['target'] > plan[i]['start']
                             else -1) for i in drift]

    # The reference returns that bracket the staircase give an independent
    # drift rate, but only the ones approached from the same side are
    # comparable; the first two both come up to zero from -5 mm.
    if results['staircase']:
        results['microstep'] = microstep_analysis(
            results['staircase'], legs, results['drift'],
            results['repeatability']['reversal_mean'])

    # Last, because it reads from every phase.
    results['campaign'] = campaign_summary(results)

    print(f"gain: {results['gain']['counts_per_full_step']:,.1f} counts/full step "
          f"= {results['gain']['counts_per_rev']:,.0f} counts/rev")
    print(f"settling moves: {len(results['settling'])}, "
          f"repeatability targets: {len(results['repeatability']['per_target'])}, "
          f"staircase legs: {len(results['staircase'])}")

    transients = []
    for i, rec in records.items():
        if plan[i]['phase'] != 'settling' or 'distance=' not in plan[i]['note']:
            continue
        a = rec['move_end'] - int(0.25 * FS)
        b = min(rec['dwell_end'], rec['move_end'] + int(1.5 * FS))
        if a < 0 or b - a < 100:
            continue
        transients.append((q[a:b] - rec['settled']) * FULL_STEP_UM)

    print('figures:')
    fig_overview(q, plan, records, args.out)
    fig_gain_raw(results['gain'], results['gain_moves'], args.out)
    fig_gain(results['gain'], args.out)
    fig_gain_dwells(results['gain_dwells'], args.out)
    fig_gain_moves(results['gain_moves'], args.out)
    if results['staircase']:
        fig_cumulative_error(results['staircase'], args.out)
        fig_microsteps(results['staircase'], args.out)
        fig_recovery(results['staircase'], results['microstep'], args.out)
        fig_microstep_map(results['microstep'], args.out)
    if results['settling']:
        fig_settling_sequence(results['settling'], results['drift'], q,
                              results['settling_model'], args.out)
        fig_settling_error(results['settling'],
                           abs(results['repeatability']['reversal_mean']) / 2,
                           args.out)
        fig_settling_traces(settling_traces(results['settling'], q), args.out)
    if results['repeatability']['per_target']:
        fig_repeat_phase(results['repeatability'], q, args.out)
        fig_repeat_sequence(results['repeatability'], q, args.out)
        fig_repeatability(results['repeatability'], args.out)
        fig_campaign(results['campaign'], args.out)
    fig_noise(results['standstill'], counts, info, args.out)

    slim = json.loads(json.dumps(results))
    slim['standstill'].pop('psd', None)
    slim['standstill'].pop('psd_freq', None)
    # The dwell traces are 2 k samples per direction; the json keeps the
    # summary of them, the figure gets the arrays in memory.
    slim['gain_dwells'].pop('t', None)
    for part in slim['gain_dwells']['by_direction'].values():
        for key in ('median', 'lo', 'hi'):
            part.pop(key, None)
    out = HERE / 'stage_results.json'
    out.write_text(json.dumps(slim, indent=1), encoding='utf-8')
    print(f'wrote {out.name}')


if __name__ == '__main__':
    main()
