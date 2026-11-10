#!/usr/bin/env python3
"""Short illustrative Maxwell-slip animation; no experimental response data."""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE / 'rendered' / 'sequence_details'
TEAL, GOLD, GREY = '#17647a', '#b47722', '#abb4ba'
FPH, FPS = 16, 20


def save_gif(frames, path):
    """Use a shared palette for stable colours and efficient frame differences."""
    for width in (960, 850, 740, 640, 540):
        size = (width, round(frames[0].height * width / frames[0].width))
        resized = [frame.convert('RGB').resize(size, Image.Resampling.LANCZOS)
                   for frame in frames]
        # Reserve both cycle colours, including their antialiased shades;
        # quantizing only the mostly-white first frame can discard gold.
        colours = [(v, v, v) for v in np.linspace(0, 255, 32).astype(int)]
        for base in ((23, 100, 122), (180, 119, 34)):
            colours.extend(tuple(round(c + (255-c)*t) for c in base)
                           for t in np.linspace(0, 1, 16))
        palette = Image.new('P', (1, 1))
        palette.putpalette([int(c) for rgb in colours for c in rgb] * 4)
        indexed = [frame.quantize(palette=palette, dither=Image.Dither.NONE)
                   for frame in resized]
        indexed[0].save(path, save_all=True, append_images=indexed[1:],
                        duration=round(1000/FPS), loop=0, optimize=True)
        if path.stat().st_size < 5_000_000:
            return
    raise RuntimeError('GIF exceeds the 5 MB limit')


def simulate(beta):
    limits = np.array([.05, .15, .4, 1., 2.5, 6.])
    weights = np.full(6, beta / (2 * limits.sum()))
    state = limits.copy()  # saturated positive by the initial approach
    strokes = [.195, 1.992, 20., 20., 1.992, .195]
    records = []
    for k, stroke in enumerate(strokes):
        u, lag = [0.], [-weights @ state]
        for target in [stroke, 0., stroke, 0.]:
            for point in np.linspace(u[-1], target, FPH + 1)[1:]:
                state = np.clip(state + point - u[-1], -limits, limits)
                u.append(point)
                lag.append(-weights @ state)
        u, lag = np.array(u), np.array(lag)
        # Directional loss on each half, then mean the two halves of a cycle.
        edges = np.arange(0, 4 * FPH + 1, FPH)
        half_loss = stroke - np.abs(np.diff((u + lag)[edges]))
        assert np.all(half_loss >= -1e-10)
        assert np.all(half_loss <= beta + 1e-10)
        assert np.isclose(lag[2 * FPH], lag[4 * FPH])
        records.append(dict(s=stroke, u=u, lag=lag,
                            loss=half_loss.reshape(2, 2).mean(axis=1)))
    return records


def main():
    beta = json.loads((HERE / 'f0b_config.json').read_text())['beta_0_um']
    records = simulate(beta)
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
                         'mathtext.fontset': 'dejavuserif', 'font.size': 10,
                         'axes.labelsize': 10, 'axes.titlesize': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})
    fig = plt.figure(figsize=(10, 6.7), dpi=110, facecolor='white')
    gs = fig.add_gridspec(2, 2, left=.09, right=.97, bottom=.20, top=.83,
                          height_ratios=[.85, 1.25], hspace=.54, wspace=.32)
    cmd = fig.add_subplot(gs[0, :])
    loop = fig.add_subplot(gs[1, 0])
    loss = fig.add_subplot(gs[1, 1])
    fig.text(.09, .965, 'L · two cycles reveal friction history', fontsize=16)
    heading = fig.text(.09, .907, '', fontsize=12, color=TEAL)
    fig.text(.09, .085, 'Illustrative Maxwell-slip response · placeholder pre-sliding scales · no F0b measurements',
             fontsize=9, color='#59636b')
    fig.text(.09, .047, f'Saturation calibrated to β₀ = {beta:.3f} µm; animation pace is not experiment time.',
             fontsize=9, color='#59636b')
    for ax in (cmd, loop, loss):
        ax.set_axisbelow(True)
        ax.grid(color='#e2e6e8', linewidth=.6)
    cmd.set(title='(a) Commanded displacement', xlabel='half-cycle index · four halves per stroke',
            ylabel='u − x₀ [µm]', xlim=(0, 24), ylim=(-.018, 36))
    cmd.set_yscale('symlog', linthresh=.04, linscale=.65)
    cmd.set_yticks([0, .195, 1.992, 20], ['0', '0.195', '1.992', '20'])
    cmd.set_xticks(np.arange(0, 25, 4))
    for k, r in enumerate(records):
        cmd.plot(k * 4 + np.arange(len(r['u'])) / FPH, r['u'], color='#d4dadd', lw=.8)
    cmd.axvline(12, ls=':', color=GREY)
    cmd.text(.25, 1.02, 'increasing strokes', transform=cmd.transAxes, ha='center', fontsize=9)
    cmd.text(.75, 1.02, 'decreasing strokes', transform=cmd.transAxes, ha='center', fontsize=9)
    progress = cmd.axvspan(0, 4, color=TEAL, alpha=.045)
    trail, = cmd.plot([], [], color=TEAL, lw=1.7)
    dot, = cmd.plot([], [], 'o', color=GOLD, ms=6)
    loop.set(title='(b) Lag loop · current stroke', xlabel='position within stroke, (u − x₀) / s',
             ylabel='lag e = x − u [µm]', xlim=(-.04, 1.04), ylim=(-1.08, 1.08))
    for val in [-beta / 2, beta / 2]:
        loop.axhline(val, color=GREY, ls=':', lw=.8)
    ghost, = loop.plot([], [], color=GREY, lw=4, alpha=.38, label='increasing-pass loop')
    first, = loop.plot([], [], '--', color=GOLD, lw=2, label='cycle 1')
    second, = loop.plot([], [], color=TEAL, lw=1.6, label='cycle 2')
    live, = loop.plot([], [], 'o', color=TEAL, ms=6)
    loop.legend(handles=[first, second, ghost], loc='upper right', frameon=True, fontsize=8)
    status = loop.text(.03, .04, '', transform=loop.transAxes, fontsize=9,
                       bbox=dict(facecolor='white', edgecolor='none', alpha=.85))
    loss.set(title='(c) Lost motion per cycle', xlabel='stroke s [µm]',
             ylabel=r'mean half-cycle loss $\ell$ [µm]', xlim=(.12, 35), ylim=(-.08, 1.97))
    loss.set_xscale('log')
    loss.axhline(beta, color=GREY, ls=':', lw=.9)
    loss.text(.14, beta + .04, f'β₀ = {beta:.3f} µm', fontsize=9, color='#59636b')
    loss.legend(handles=[Line2D([], [], color=GOLD, marker='s', ls='none', label='cycle 1'),
                         Line2D([], [], color=TEAL, marker='o', ls='none', label='cycle 2')],
                loc='upper left', bbox_to_anchor=(0, .88), fontsize=9, frameon=False)
    loss.text(.03, .51, 'filled: increasing\nopen: decreasing', transform=loss.transAxes, fontsize=8)
    marks = []
    frames = []
    for k, r in enumerate(records):
        progress.remove()
        progress = cmd.axvspan(4*k, 4*k+4, color=TEAL, alpha=.045)
        if k >= 3:
            old = records[5-k]
            ghost.set_data(old['u'] / old['s'], old['lag'])
        else:
            ghost.set_data([], [])
        for j in range(0, 4 * FPH + 1):
            cycle = 1 if j <= 2 * FPH else 2
            regime = 'small stroke' if r['s'] < 1 else 'pre-sliding transition' if r['s'] < 10 else 'sliding saturation'
            heading.set_text(f"{r['s']:g} µm · {regime} · {'increasing' if k < 3 else 'decreasing'} · cycle {cycle}/2")
            trail.set_data(k*4 + np.arange(j+1)/FPH, r['u'][:j+1])
            dot.set_data([k*4+j/FPH], [r['u'][j]])
            end = min(j, 2*FPH)
            first.set_data(r['u'][:end+1]/r['s'], r['lag'][:end+1])
            second.set_data(r['u'][2*FPH:j+1]/r['s'], r['lag'][2*FPH:j+1])
            live.set_data([r['u'][j]/r['s']], [r['lag'][j]])
            status.set_text('First cycle inherits the starting state' if cycle == 1 else
                            ('Closed loop; history may shift its offset' if k >= 3 else 'Second cycle closes the repeat loop'))
            if j in (2*FPH, 4*FPH):
                c = j // (2*FPH) - 1
                marks.append(loss.plot(r['s'], r['loss'][c], marker='s' if c == 0 else 'o',
                                      ms=6 if c == 0 else 9, mfc=(GOLD if c == 0 else TEAL) if k < 3 else 'white',
                                      mec=GOLD if c == 0 else TEAL, mew=1.4, ls='none')[0])
            fig.canvas.draw()
            frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
            frames.append(frame)
            if j in (2*FPH, 4*FPH):
                frames.extend([frame.copy() for _ in range(8)])
    path = OUT / 'f0b_L_friction_loops.gif'
    save_gif(frames, path)
    print(f'{path}\n{len(frames)} frames; {len(frames)/FPS:.1f} s; {path.stat().st_size/1e6:.1f} MB')
    plt.close(fig)


if __name__ == '__main__':
    main()
