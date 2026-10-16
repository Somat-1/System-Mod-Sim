#!/usr/bin/env python3
'''Render the electrical-alignment acceptance check from the motor campaign.

The figure uses actions 13 and 14 (64 native 1/16 microsteps outward and 64
back). Encoder direction is inferred before any metric is calculated; a
negative command/encoder covariance causes the encoder sign to be flipped.
Settled positions are medians over the final quarter of each inter-step dwell.
'''

import argparse
import json
import re
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
MOTOR_DIR = HERE.parent
DEFAULT_RECORDING = (MOTOR_DIR / 'data' /
                     '2026-10-13_motor_identification' /
                     'motor_identification_run.csv')
DEFAULT_OUT = MOTOR_DIR / 'analysis' / 'motor_encoder_report' / 'plots'

BLUE = '#2a78d6'
ORANGE = '#e96b36'
PURPLE = '#6653bd'
INK = '#161616'
MUTED = '#66645f'
GRID = '#deddd8'
GREEN = '#16784a'
RED = '#bd3a32'


def style():
    mpl.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'font.family': 'DejaVu Sans', 'font.size': 8.5,
        'axes.titlesize': 10, 'axes.titleweight': 'bold',
        'axes.labelsize': 8.5, 'axes.edgecolor': GRID,
        'axes.linewidth': 0.8, 'axes.grid': True,
        'axes.axisbelow': True, 'grid.color': GRID, 'grid.linewidth': 0.55,
        'xtick.color': MUTED, 'ytick.color': MUTED,
        'xtick.labelsize': 7.5, 'ytick.labelsize': 7.5,
        'legend.frameon': False, 'legend.fontsize': 7.3,
        'pdf.fonttype': 42, 'ps.fonttype': 42,
    })


def load_alignment(path):
    '''Read only alignment actions 13--14 and the startup driver readback.'''
    driver = {}
    rows = []
    driver_pattern = re.compile(r'([a-z_]+)=([^,]+)')
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('# DRIVER'):
                driver = {key: value.strip() for key, value
                          in driver_pattern.findall(line)}
            elif line.startswith('D,'):
                fields = line.rstrip().split(',')
                action = int(fields[3])
                if action in (13, 14):
                    rows.append((int(fields[1]), int(fields[2]), action,
                                 int(fields[4]), int(fields[5]), int(fields[6]),
                                 int(fields[7]), int(fields[8]), int(fields[9])))
            elif line.startswith('# ACTION_END,index=14'):
                break
    if not rows:
        raise ValueError(f'no alignment telemetry (actions 13--14) in {path}')
    return np.asarray(rows, dtype=np.int64), driver


def settled_positions(data, cpr, full_steps):
    '''Return settled command/encoder pairs following every marked step edge.'''
    event_rows = np.flatnonzero(data[:, 3] != 0)
    if event_rows.size == 0:
        raise ValueError('alignment telemetry contains no marked step edges')

    baseline_count = float(data[event_rows[0], 5])
    command_u256, encoder_counts, direction = [], [], []
    for number, row in enumerate(event_rows):
        stop = event_rows[number + 1] if number + 1 < event_rows.size else len(data)
        # Use the final quarter of the hold so transient motion is excluded.
        start = row + max(1, int((stop - row) * 0.75))
        window = data[start:stop, 5]
        if window.size == 0:
            window = data[row:stop, 5]
        command_u256.append(data[row, 4])
        encoder_counts.append(float(np.median(window)) - baseline_count)
        direction.append(data[row, 3])

    command_deg = np.asarray(command_u256) * 360.0 / (full_steps * 256.0)
    encoder_deg = np.asarray(encoder_counts) * 360.0 / cpr
    direction = np.asarray(direction, dtype=int)
    covariance = float(np.cov(command_deg, encoder_deg, bias=True)[0, 1])
    encoder_sign = 1 if covariance >= 0 else -1
    return (event_rows, command_deg, encoder_deg * encoder_sign, direction,
            encoder_sign, baseline_count)


def analyse(data, driver, cpr, full_steps, expected_mres, setpoint_ma):
    event_rows, command_deg, encoder_deg, direction, sign, baseline = \
        settled_positions(data, cpr, full_steps)
    step_deg = 360.0 / (full_steps * expected_mres)
    slope, intercept = np.polyfit(command_deg, encoder_deg, 1)
    error_deg = encoder_deg - command_deg

    outward = direction == 1
    returning = direction == -1
    out_map = {round(c, 9): e for c, e in
               zip(command_deg[outward], encoder_deg[outward])}
    ret_map = {round(c, 9): e for c, e in
               zip(command_deg[returning], encoder_deg[returning])}
    shared = sorted(set(out_map) & set(ret_map))
    opening = np.asarray([out_map[key] - ret_map[key] for key in shared])

    mode_ok = (driver.get('spreadcycle') == '1' and
               driver.get('intpol') == '0' and
               driver.get('stealth') == '0' and
               driver.get('coolstep') == '0' and
               driver.get('ihold') == driver.get('irun'))
    mres_values = np.unique(data[:, 8])
    current_values = np.unique(data[:, 7])
    positive_events = int(np.sum(data[:, 3] == 1))
    negative_events = int(np.sum(data[:, 3] == -1))
    direction_ok = sign == 1
    scale_ok = abs(slope - 1.0) <= 0.005
    return_error = float(abs(encoder_deg[-1] - command_deg[-1]))
    index_values = np.unique(data[:, 6])
    counts_per_microstep = cpr / (full_steps * expected_mres)

    checks = [
        ('Mode readback lines',
         f"spread=1; intpol=0; IHOLD=IRUN={driver.get('irun', '?')}",
         'SpreadCycle; intpol=0; IHOLD=IRUN confirmed', mode_ok),
        ('mres column',
         ', '.join(str(value) for value in mres_values),
         f'{expected_mres} throughout',
         bool(np.array_equal(mres_values, [expected_mres]))),
        ('current_ma',
         ', '.join(f'{value} mA' for value in current_values),
         f'close to {setpoint_ma} mA setpoint',
         bool(np.all(np.abs(current_values - setpoint_ma) <= 10))),
        ('Step events', f'+1: {positive_events}; −1: {negative_events}',
         '64 outward; 64 back', positive_events == 64 and negative_events == 64),
        ('Direction sign',
         'same sign; no flip' if sign == 1 else 'opposite; encoder flipped',
         'encoder and command same sign', direction_ok),
        ('Fitted scale', f'{slope:.4f} ({100 * (slope - 1):+.3f}%)',
         '1.000 ± ~0.5%', scale_ok),
        ('Return error', f'{return_error:.3f}° ({return_error / step_deg:.2f} µstep)',
         f'< ~1 microstep ({step_deg:.3f}°)', return_error < step_deg),
        ('index_count',
         f'{index_values[0]} → {index_values[-1]} (Δ{index_values[-1] - index_values[0]:+d})',
         'no change, or expected Z change', index_values.size == 1),
        ('Quantisation', f'{counts_per_microstep:.2f} counts/microstep',
         '≥ ~10 counts per microstep', counts_per_microstep >= 10),
    ]
    return dict(event_rows=event_rows, command_deg=command_deg,
                encoder_deg=encoder_deg, direction=direction, sign=sign,
                baseline=baseline, slope=float(slope), intercept=float(intercept),
                error_deg=error_deg, loop_mean_abs=float(np.mean(np.abs(opening))),
                loop_max_abs=float(np.max(np.abs(opening))),
                return_error=return_error, step_deg=step_deg,
                counts_per_microstep=counts_per_microstep, checks=checks)


def despine(ax):
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)


def render(data, metrics, cpr, full_steps, source, dense=False):
    fig = plt.figure(figsize=(13.2, 9.0))
    grid = fig.add_gridspec(2, 3, height_ratios=(1.0, 0.79),
                            left=0.065, right=0.985, bottom=0.045, top=0.875,
                            wspace=0.30, hspace=0.34)
    axes = [fig.add_subplot(grid[0, column]) for column in range(3)]
    table_ax = fig.add_subplot(grid[1, :])

    failed = sum(not row[3] for row in metrics['checks'])
    verdict = 'PASS' if failed == 0 else f'FAIL — {failed} of {len(metrics["checks"])} checks'
    verdict_colour = GREEN if failed == 0 else RED
    fig.suptitle('Electrical alignment validation', x=0.065, y=0.968,
                 ha='left', fontsize=16, fontweight='bold', color=INK)
    fig.text(0.065, 0.925,
             'Actions 13–14 · 64 microsteps outward and 64 return · '
             f'encoder CPR = {cpr:,}', ha='left', color=MUTED, fontsize=9)
    fig.text(0.985, 0.948, verdict, ha='right', va='center', color='white',
             fontsize=9, fontweight='bold',
             bbox=dict(boxstyle='round,pad=0.38', facecolor=verdict_colour,
                       edgecolor='none'))

    # (a) Command and encoder against time.
    ax = axes[0]
    time_s = (data[:, 1] - data[0, 1]) / 1e6
    command_raw = data[:, 4] * 360.0 / (full_steps * 256.0)
    encoder_raw = ((data[:, 5] - metrics['baseline']) * metrics['sign'] *
                   360.0 / cpr)
    ax.step(time_s, command_raw, where='post', color=BLUE, linewidth=1.25,
            label='command')
    ax.plot(time_s, encoder_raw, color=ORANGE, linewidth=1.0, label='encoder')
    events = metrics['event_rows']
    event_sign = data[events, 3]
    for sign, colour, y, label in ((1, BLUE, 0.015, 'event = +1'),
                                   (-1, ORANGE, 0.045, 'event = −1')):
        chosen = events[event_sign == sign]
        ax.plot(time_s[chosen], np.full(chosen.size, y), linestyle='none',
                marker='|', markersize=4.0, markeredgewidth=0.55, color=colour,
                transform=ax.get_xaxis_transform(), clip_on=False, label=label)
    ax.set_xlabel('time from first alignment edge [s]')
    ax.set_ylabel('relative shaft angle [deg]')
    ax.set_title('(a) Command and encoder vs time', loc='left')
    ax.legend(loc='upper right', ncol=2)
    ax.text(0.025, 0.95,
            f"peak command {np.max(command_raw):.2f}°\n"
            f"direction {'unchanged' if metrics['sign'] == 1 else 'sign-corrected'}",
            transform=ax.transAxes, ha='left', va='top', color=MUTED,
            fontsize=7.3)
    despine(ax)

    # (b) Settled encoder positions against command.
    ax = axes[1]
    command = metrics['command_deg']
    encoder = metrics['encoder_deg']
    direction = metrics['direction']
    if dense:
        ax.scatter(command_raw, encoder_raw, s=3.0, color=MUTED, alpha=0.14,
                   linewidths=0, rasterized=True, label='recorded samples')
    for sign, colour, label in ((1, BLUE, 'outward'), (-1, ORANGE, 'return')):
        keep = direction == sign
        ax.plot(command[keep], encoder[keep], color=colour, linewidth=1.0,
                marker='o', markersize=2.4, markeredgewidth=0, label=label)
    limit = max(np.max(command), np.max(encoder)) * 1.02
    line = np.linspace(0, limit, 200)
    ax.plot(line, line, color=MUTED, linestyle=':', linewidth=1.0,
            label='ideal 1:1')
    ax.plot(line, metrics['slope'] * line + metrics['intercept'], color=PURPLE,
            linestyle='--', linewidth=1.1, label='least-squares fit')
    ax.set_xlim(-0.12, limit)
    ax.set_ylim(-0.20, limit)
    ax.set_aspect('equal', adjustable='box')
    ax.set_xlabel('commanded angle [deg]')
    ax.set_ylabel('encoder angle [deg]')
    ax.set_title('(b) Encoder vs command', loc='left')
    ax.legend(loc='lower right')
    ax.text(0.04, 0.96,
            f"fit: encoder = {metrics['slope']:.4f} command "
            f"{metrics['intercept']:+.3f}°\n"
            f"loop |opening|: mean {metrics['loop_mean_abs']:.3f}°, "
            f"max {metrics['loop_max_abs']:.3f}°",
            transform=ax.transAxes, ha='left', va='top', fontsize=7.2,
            color=MUTED)
    despine(ax)

    # (c) Direct position error; no fit is removed, so drift remains visible.
    ax = axes[2]
    if dense:
        ax.scatter(command_raw, encoder_raw - command_raw, s=3.0, color=MUTED,
                   alpha=0.14, linewidths=0, rasterized=True,
                   label='recorded samples')
    for sign, colour, label in ((1, BLUE, 'outward'), (-1, ORANGE, 'return')):
        keep = direction == sign
        ax.plot(command[keep], metrics['error_deg'][keep], color=colour,
                linewidth=1.05, marker='o', markersize=2.4,
                markeredgewidth=0, label=label)
    ax.axhline(0, color=MUTED, linewidth=0.8)
    ax.axhline(metrics['step_deg'], color=MUTED, linewidth=0.7, linestyle=':')
    ax.axhline(-metrics['step_deg'], color=MUTED, linewidth=0.7, linestyle=':',
               label='±1 microstep')
    ax.set_xlim(-0.12, limit)
    ax.set_xlabel('commanded angle [deg]')
    ax.set_ylabel(r'$\theta_{enc}-\theta_{cmd}$ [deg]')
    ax.set_title('(c) Error vs commanded position', loc='left')
    ax.legend(loc='lower left')
    ax.text(0.98, 0.96,
            f"range {np.min(metrics['error_deg']):+.3f}° to "
            f"{np.max(metrics['error_deg']):+.3f}°",
            transform=ax.transAxes, ha='right', va='top', fontsize=7.2,
            color=MUTED)
    despine(ax)

    # Pass/fail table.
    table_ax.axis('off')
    table_ax.text(0.0, 0.985, 'Alignment acceptance table', transform=table_ax.transAxes,
                  ha='left', va='top', fontsize=10, fontweight='bold', color=INK)
    table_ax.text(1.0, 0.985,
                  'Settled value = median of final quarter of each 120 ms dwell',
                  transform=table_ax.transAxes, ha='right', va='top',
                  fontsize=7.2, color=MUTED)
    rows = [[check, observed, expected, 'PASS' if passed else 'FAIL']
            for check, observed, expected, passed in metrics['checks']]
    table = table_ax.table(cellText=rows,
                           colLabels=['Check', 'Observed', 'Expected', 'Status'],
                           colWidths=[0.19, 0.27, 0.44, 0.10],
                           cellLoc='left', colLoc='left',
                           bbox=[0.0, 0.0, 1.0, 0.90])
    table.auto_set_font_size(False)
    table.set_fontsize(7.15)
    for (row, column), cell in table.get_celld().items():
        cell.set_edgecolor('#d7d6d2')
        cell.set_linewidth(0.45)
        cell.PAD = 0.16
        if row == 0:
            cell.set_facecolor('#ecebe7')
            cell.get_text().set_fontweight('bold')
            cell.get_text().set_color(INK)
        else:
            cell.set_facecolor('#fafafa' if row % 2 else 'white')
            if column == 3:
                passed = metrics['checks'][row - 1][3]
                cell.set_facecolor('#e1f1e8' if passed else '#f7e2df')
                cell.get_text().set_color(GREEN if passed else RED)
                cell.get_text().set_fontweight('bold')

    fig.text(0.985, 0.012, f'Source: {source.name}', ha='right', va='bottom',
             fontsize=6.8, color=MUTED)
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recording', type=Path, default=DEFAULT_RECORDING)
    parser.add_argument('--out', type=Path, default=DEFAULT_OUT)
    parser.add_argument('--cpr', type=int, default=20_000)
    parser.add_argument('--full-steps', type=int, default=200)
    parser.add_argument('--mres', type=int, default=16)
    parser.add_argument('--current-ma', type=int, default=360)
    parser.add_argument('--dense', action='store_true',
                        help='show every recorded sample behind settled points')
    parser.add_argument('--interactive', action='store_true',
                        help='open the finished figure in an interactive window')
    args = parser.parse_args()

    style()
    data, driver = load_alignment(args.recording)
    metrics = analyse(data, driver, args.cpr, args.full_steps,
                      args.mres, args.current_ma)
    args.out.mkdir(parents=True, exist_ok=True)
    fig = render(data, metrics, args.cpr, args.full_steps, args.recording,
                 dense=args.dense or args.interactive)

    png = args.out / 'alignment_validation.png'
    pdf = args.out / 'alignment_validation.pdf'
    fig.savefig(png, dpi=240)
    fig.savefig(pdf)
    machine = {
        'source': str(args.recording), 'cpr': args.cpr,
        'encoder_sign_applied': metrics['sign'],
        'fit_slope': metrics['slope'],
        'fit_intercept_deg': metrics['intercept'],
        'loop_mean_abs_deg': metrics['loop_mean_abs'],
        'loop_max_abs_deg': metrics['loop_max_abs'],
        'return_error_deg': metrics['return_error'],
        'counts_per_microstep': metrics['counts_per_microstep'],
        'checks': [{'check': row[0], 'observed': row[1],
                    'expected': row[2], 'pass': bool(row[3])}
                   for row in metrics['checks']],
    }
    metrics_path = args.out / 'alignment_validation_metrics.json'
    metrics_path.write_text(json.dumps(machine, indent=2), encoding='utf-8')
    print(f'wrote {png}')
    print(f'wrote {pdf}')
    print(f'wrote {metrics_path}')
    if args.interactive:
        print('opening interactive Matplotlib window')
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
