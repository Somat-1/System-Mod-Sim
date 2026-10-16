#!/usr/bin/env python3
'''Open the complete motor-identification recording in Matplotlib.

Every one of the 1.9 million telemetry rows is used for the command, encoder,
and raw error traces.  No settling selection, fit, filtering, eccentricity
correction, or calibration is applied.  The lower panels show the recorded
driver settings, index count, action number, and all marked STEP events.
'''

import argparse
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


HERE = Path(__file__).resolve().parent
DEFAULT_RECORDING = (HERE.parent / 'data' /
                     '2026-10-13_motor_identification' /
                     'motor_identification_run.csv')

BLUE = '#2a78d6'
ORANGE = '#e96b36'
GREEN = '#1b9e77'
PURPLE = '#6653bd'
GREY = '#65635f'
GRID = '#deddd8'


def load(path):
    table = np.loadtxt(path, delimiter=',', comments='#', skiprows=3,
                       usecols=range(1, 10), dtype=np.int64)
    names = ('seq', 't_us', 'action', 'event', 'command_u256', 'encoder',
             'index_count', 'current_ma', 'mres')
    return {name: table[:, column] for column, name in enumerate(names)}


def change_trace(time, values):
    '''Keep first, last, and both sides of each change for a faithful step line.'''
    change = np.flatnonzero(np.diff(values) != 0)
    indices = np.unique(np.r_[0, change, change + 1, len(values) - 1])
    return time[indices], values[indices]


def padded_limits(*arrays, fraction=0.07):
    '''Finite min/max with a useful margin, including constant-valued windows.'''
    finite = [np.asarray(array)[np.isfinite(array)] for array in arrays]
    finite = [array for array in finite if array.size]
    if not finite:
        return -1.0, 1.0
    low = min(float(np.min(array)) for array in finite)
    high = max(float(np.max(array)) for array in finite)
    span = high - low
    pad = max(span * fraction, max(abs(low), abs(high), 1.0) * 0.015)
    return low - pad, high + pad


class VisibleWindowAutoscaler:
    '''Keep y limits tied to the samples visible after an x zoom or pan.'''

    def __init__(self, fig, axes, twins, time, series):
        self.fig = fig
        self.axes = axes
        self.twins = twins
        self.time = time
        self.series = series
        self.busy = False
        self.last_xlim = None

    def __call__(self, _axis=None, force=False):
        if self.busy:
            return
        xlim = tuple(sorted(self.axes[0].get_xlim()))
        if not force and self.last_xlim == xlim:
            return
        self.last_xlim = xlim
        first = max(0, int(np.searchsorted(self.time, xlim[0], side='left')) - 1)
        last = min(len(self.time),
                   int(np.searchsorted(self.time, xlim[1], side='right')) + 1)
        if last <= first:
            return
        sl = slice(first, last)
        self.busy = True
        try:
            self.axes[0].set_ylim(*padded_limits(
                self.series['command'][sl], self.series['encoder'][sl]))
            self.axes[1].set_ylim(*padded_limits(self.series['error'][sl]))
            self.axes[2].set_ylim(*padded_limits(self.series['current'][sl]))
            self.twins[0].set_ylim(*padded_limits(self.series['mres'][sl]))
            self.axes[3].set_ylim(*padded_limits(self.series['index'][sl]))
            self.twins[1].set_ylim(*padded_limits(self.series['action'][sl]))
            self.fig.canvas.draw_idle()
        finally:
            self.busy = False

    def key_press(self, event):
        if event.key and event.key.lower() == 'y':
            self(event.inaxes, force=True)


def make_figure(data, cpr, full_steps):
    mpl.rcParams.update({
        'figure.facecolor': 'white', 'axes.facecolor': 'white',
        'font.family': 'DejaVu Sans', 'font.size': 9,
        'axes.grid': True, 'axes.axisbelow': True,
        'grid.color': GRID, 'grid.linewidth': 0.55,
        'axes.edgecolor': GRID, 'axes.linewidth': 0.8,
        'xtick.color': GREY, 'ytick.color': GREY,
        'legend.frameon': False, 'legend.fontsize': 8,
        'agg.path.chunksize': 20_000,
        'path.simplify': True, 'path.simplify_threshold': 0.35,
    })

    time_min = ((data['t_us'] - data['t_us'][0]) / 60e6).astype(np.float32)
    command_deg = (data['command_u256'] *
                   (360.0 / (full_steps * 256.0))).astype(np.float32)
    encoder_deg = (data['encoder'] * (360.0 / cpr)).astype(np.float32)
    error_deg = (encoder_deg - command_deg).astype(np.float32)

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(15.5, 9.2),
                             gridspec_kw={'height_ratios': [2.0, 1.35, 0.9, 0.9]})
    fig.subplots_adjust(left=0.065, right=0.925, top=0.89, bottom=0.075,
                        hspace=0.16)
    fig.suptitle('Complete motor-identification measurement', x=0.065, y=0.965,
                 ha='left', fontsize=16, fontweight='bold')
    fig.text(0.065, 0.925,
             f'{len(time_min):,} recorded samples · '
             f'{time_min[-1]:.2f} min · no filtering or point reduction',
             ha='left', fontsize=9, color=GREY)

    # All recorded samples: command and encoder angle.
    ax = axes[0]
    ax.plot(time_min, command_deg, color=BLUE, linewidth=0.55,
            label='command_u256 → angle')
    ax.plot(time_min, encoder_deg, color=ORANGE, linewidth=0.50, alpha=0.90,
            label='encoder counts → angle')
    ax.set_ylabel('shaft angle [deg]')
    ax.set_title('(a) Raw command and encoder', loc='left', fontweight='bold')
    ax.legend(loc='upper right', ncol=2)

    # Direct raw error; no scale or static-map fit removed.
    ax = axes[1]
    ax.plot(time_min, error_deg, color=PURPLE, linewidth=0.48)
    ax.axhline(0, color=GREY, linewidth=0.7)
    ax.set_ylabel(r'$\theta_{enc}-\theta_{cmd}$ [deg]')
    ax.set_title('(b) Raw encoder-minus-command error', loc='left',
                 fontweight='bold')

    # Settings are piecewise constant, so only exact change vertices are needed.
    ax = axes[2]
    tc, current = change_trace(time_min, data['current_ma'])
    tm, mres = change_trace(time_min, data['mres'])
    ax.step(tc, current, where='post', color=GREEN, linewidth=1.0,
            label='current_ma')
    ax.set_ylabel('current [mA]', color=GREEN)
    ax.tick_params(axis='y', colors=GREEN)
    twin = ax.twinx()
    twin.step(tm, mres, where='post', color=ORANGE, linewidth=1.0,
              label='MRES denominator')
    twin.set_ylabel('MRES denominator', color=ORANGE)
    twin.tick_params(axis='y', colors=ORANGE)
    ax.set_title('(c) Recorded current and microstep resolution', loc='left',
                 fontweight='bold')
    lines = ax.get_lines() + twin.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc='upper right',
              ncol=2)

    # Recorded sequence metadata and every event-marked STEP edge.
    ax = axes[3]
    ti, index = change_trace(time_min, data['index_count'])
    ta, action = change_trace(time_min, data['action'])
    ax.step(ti, index, where='post', color=GREY, linewidth=0.9,
            label='index_count')
    ax.set_ylabel('index count', color=GREY)
    action_ax = ax.twinx()
    action_ax.step(ta, action, where='post', color=PURPLE, linewidth=0.9,
                   label='action')
    action_ax.set_ylabel('action number', color=PURPLE)
    action_ax.tick_params(axis='y', colors=PURPLE)
    event = data['event'] != 0
    plus = event & (data['event'] > 0)
    minus = event & (data['event'] < 0)
    ax.plot(time_min[plus], np.full(np.sum(plus), 0.02), linestyle='none',
            marker='|', markersize=3.5, markeredgewidth=0.45, color=BLUE,
            transform=ax.get_xaxis_transform(), label='event +1')
    ax.plot(time_min[minus], np.full(np.sum(minus), 0.055), linestyle='none',
            marker='|', markersize=3.5, markeredgewidth=0.45, color=ORANGE,
            transform=ax.get_xaxis_transform(), label='event −1')
    ax.set_title('(d) Index, action, and all marked STEP events', loc='left',
                 fontweight='bold')
    lines = ax.get_lines() + action_ax.get_lines()
    ax.legend(lines, [line.get_label() for line in lines], loc='upper left',
              ncol=4)
    ax.set_xlabel('elapsed recording time [min]')
    ax.set_xlim(float(time_min[0]), float(time_min[-1]))

    for axis in axes:
        axis.spines['top'].set_visible(False)
        axis.spines['right'].set_visible(False)
    for axis in (twin, action_ax):
        axis.spines['top'].set_visible(False)

    fig.text(0.925, 0.018,
             'Horizontal zoom/pan auto-rescales every y-axis · press Y to refresh · '
             'toolbar: home, back/forward, configure, save',
             ha='right', fontsize=7.5, color=GREY)

    autoscaler = VisibleWindowAutoscaler(
        fig, axes, (twin, action_ax), time_min,
        dict(command=command_deg, encoder=encoder_deg, error=error_deg,
             current=data['current_ma'], mres=data['mres'],
             index=data['index_count'], action=data['action']))
    for axis in axes:
        axis.callbacks.connect('xlim_changed', autoscaler)
    fig.canvas.mpl_connect('key_press_event', autoscaler.key_press)
    # Keep the callback owner alive for the lifetime of the interactive figure.
    fig._visible_window_autoscaler = autoscaler
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recording', type=Path, default=DEFAULT_RECORDING)
    parser.add_argument('--cpr', type=int, default=20_000)
    parser.add_argument('--full-steps', type=int, default=200)
    parser.add_argument('--save', type=Path,
                        help='optional image path for a non-interactive copy')
    parser.add_argument('--no-show', action='store_true',
                        help='render without opening the GUI')
    args = parser.parse_args()

    print(f'loading {args.recording}')
    data = load(args.recording)
    print(f'loaded {len(data["seq"]):,} rows')
    fig = make_figure(data, args.cpr, args.full_steps)
    if args.save:
        args.save.parent.mkdir(parents=True, exist_ok=True)
        fig.savefig(args.save, dpi=180)
        print(f'wrote {args.save}')
    if not args.no_show:
        print('opening complete-recording Matplotlib window')
        plt.show()
    plt.close(fig)


if __name__ == '__main__':
    main()
