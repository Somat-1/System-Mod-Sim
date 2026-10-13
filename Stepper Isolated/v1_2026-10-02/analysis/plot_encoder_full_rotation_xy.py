#!/usr/bin/env python3
'''XY plot of encoder measurement versus reconstructed stepper command.'''

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

HERE = Path(__file__).resolve().parent
DEFAULT_CSV = HERE.parent / 'recordings' / 'encoder_full_rotation_20261009_162013.csv'
DEFAULT_PNG = HERE / 'plots' / 'encoder_full_rotation_xy.png'


def load_csv(path):
    '''Return comment metadata and numeric sample columns.'''
    meta, lines = {}, []
    with path.open(encoding='utf-8') as stream:
        for line in stream:
            if line.startswith('#'):
                key, sep, value = line[1:].partition(':')
                if sep:
                    meta[key.strip()] = value.strip()
            elif line.strip():
                lines.append(line)
    rows = list(csv.DictReader(lines))
    if not rows:
        raise ValueError(f'No samples found in {path}')
    t_us = np.array([float(row['t_us']) for row in rows])
    count = np.array([float(row['count']) for row in rows])
    return meta, t_us, count


def make_plot(csv_path, output_path, encoder_cpr, full_steps, microsteps, duration,
              deviation_gain):
    meta, t_us, count = load_csv(csv_path)
    elapsed = (t_us - t_us[0]) * 1e-6
    recorded_duration = float(meta.get('clip_duration_s', elapsed[-1]))
    duration = recorded_duration if duration is None else duration

    commands = full_steps * microsteps
    step_angle = 360.0 / commands
    command_index = np.rint(np.clip(elapsed / duration, 0, 1) * commands)
    commanded = command_index * step_angle
    measured = (count - count[0]) * 360.0 / encoder_cpr
    error = measured - commanded
    mean_error = float(np.mean(error))
    rms_error = float(np.sqrt(np.mean(error**2)))
    endpoint_error = float(error[-1])
    max_positive_index = int(np.argmax(error))
    max_negative_index = int(np.argmin(error))
    displayed = commanded + deviation_gain * error

    plt.rcParams.update({'font.size': 10, 'axes.titlesize': 14})
    fig, ax = plt.subplots(figsize=(8.2, 7.2), constrained_layout=True)
    low = min(0.0, float(displayed.min()))
    high = max(360.0, float(displayed.max()))
    margin = max(4.0, 0.02 * (high - low))
    limits = (low - margin, high + margin)

    ax.plot(commanded, displayed, color='#0072B2', linewidth=1.15,
            label=f'Encoder measurement (deviation x{deviation_gain:g})', zorder=3)
    ax.plot(limits, limits, '--', color='#222222', linewidth=1.4,
            label='Ideal (measured = commanded)', zorder=2)
    ax.set(xlim=limits, ylim=limits, xlabel='Commanded angle [deg]',
           ylabel='Displayed encoder angle [deg]',
           title=f'Encoder Full-Rotation Tracking ({deviation_gain:g}x deviation)')
    ax.set_aspect('equal', adjustable='box')
    ax.grid(True, color='#D7D7D7', linewidth=0.7, alpha=0.8)

    for index, color, description in (
        (max_positive_index, '#D55E00', 'Maximum positive deviation'),
        (max_negative_index, '#009E73', 'Maximum negative deviation'),
    ):
        x_value = commanded[index]
        actual_error = error[index]
        displayed_value = displayed[index]
        ax.plot([x_value, x_value], [x_value, displayed_value], color=color,
                linewidth=2.3, marker='o', markersize=5, zorder=5,
                label=(f'{description}: {actual_error:+.4f} deg '
                       f'at {x_value:.4f} deg command'))

    ax.legend(loc='upper left')

    details = (
        f'{full_steps} full steps/rev, 1/{microsteps} microstepping\n'
        f'{commands:,} commands at {step_angle:.4f} deg/microstep\n'
        f'Duration: {recorded_duration:.6f} s; encoder: {encoder_cpr:,.0f} counts/rev\n'
        f'Encoder deviations magnified {deviation_gain:g}x for visibility only\n'
        f'Error: mean {mean_error:+.3f} deg, RMS {rms_error:.3f} deg, '
        f'endpoint {endpoint_error:+.3f} deg'
    )
    ax.text(0.98, 0.03, details, transform=ax.transAxes, ha='right', va='bottom',
            fontsize=8.5, bbox={'boxstyle': 'round,pad=0.45', 'facecolor': 'white',
                               'alpha': 0.9, 'edgecolor': '#AAAAAA'})

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=240, bbox_inches='tight')
    fig.savefig(output_path.with_suffix('.svg'), bbox_inches='tight')
    plt.close(fig)
    print(f'Samples: {len(elapsed):,}')
    print(f'Measured endpoint: {measured[-1]:.3f} deg')
    print(f'Error: mean={mean_error:+.4f}, RMS={rms_error:.4f}, '
          f'endpoint={endpoint_error:+.4f} deg')
    print(f'Wrote: {output_path}')
    print(f'''Wrote: {output_path.with_suffix('.svg')}''')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, default=DEFAULT_CSV)
    parser.add_argument('--output', type=Path, default=DEFAULT_PNG)
    parser.add_argument('--encoder-cpr', type=float, default=20_000.0)
    parser.add_argument('--full-steps', type=int, default=200)
    parser.add_argument('--microsteps', type=int, default=16)
    parser.add_argument('--duration', type=float)
    parser.add_argument('--deviation-gain', type=float, default=15.0,
                        help='Visual magnification applied about the ideal line')
    args = parser.parse_args()
    make_plot(args.input.resolve(), args.output.resolve(), args.encoder_cpr,
              args.full_steps, args.microsteps, args.duration,
              args.deviation_gain)


if __name__ == '__main__':
    main()
