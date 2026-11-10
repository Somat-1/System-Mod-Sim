#!/usr/bin/env python3
"""Illustrative S slow-ramp animation (spring-block mechanism + lag); no experimental data.

S1 (wipe from below, ramp down: a reversal) and S4 (wipe from above, ramp
down: a continuation) cover the same 20 um at 1/256. A quasi-static
Maxwell-slip model in displacement units, u - x = sum(gamma_i z_i) with
|z_i| <= a_i, is solved exactly from slip event to slip event. Slip limits
are placeholders; the saturated lag is calibrated to beta0/2.
"""
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent
OUT = HERE / 'rendered' / 'sequence_details'
OCHRE, DARK, GREY, INK = '#c68e37', '#2b3640', '#abb4ba', '#59636b'
FPS, SUB = 20, 4                       # frames per second; frames per microstep in 1b
USTEP = 2000 / (200 * 256)             # 39.06 nm in um
STEPS, RATE = 512, 6.4                 # 20 um ramp at 6.4 usteps/s
LIMITS = np.array([.02, .08, .25, .8, 2.5])          # slip limits a_i [um]
GAMMA0 = np.array([7.7, 2.6, .96, .26, .03])         # relative stiffness k_i / K
WIPE_UM = 6.0                          # schematic wipe length shown in (a)

# Panel (a) world-to-screen: one scale for lag and travel, so the cartoon is
# proportional; the view pans to keep the carriage in frame.
SCALE, L0, BLOCK_W = 2.6, 3.5, 2.4     # units per um, spring rest length, block width
BLOCK_SCREEN, PAN_EDGE = 10.0, 7.5      # block left edge at ramp start; pan threshold


def save_gif(frames, path):
    """Shared palette (greys + ochre + slate) for stable colours and small diffs."""
    for width in (1000, 900, 800, 700, 600):
        size = (width, round(frames[0].height * width / frames[0].width))
        resized = [f.convert('RGB').resize(size, Image.Resampling.LANCZOS) for f in frames]
        colours = [(v, v, v) for v in np.linspace(0, 255, 32).astype(int)]
        for base in ((198, 142, 55), (43, 54, 64)):
            colours.extend(tuple(round(c + (255 - c) * t) for c in base)
                           for t in np.linspace(0, 1, 16))
        palette = Image.new('P', (1, 1))
        palette.putpalette([int(c) for rgb in colours for c in rgb] * 4)
        indexed = [f.quantize(palette=palette, dither=Image.Dither.NONE) for f in resized]
        indexed[0].save(path, save_all=True, append_images=indexed[1:],
                        duration=round(1000 / FPS), loop=0, optimize=True)
        if path.stat().st_size < 5_000_000:
            return
    raise RuntimeError('GIF exceeds the 5 MB limit')


def advance(z, du, gamma):
    """Move the command by du; return new z and carriage displacement dx (exact)."""
    d, remaining, dx_total = np.sign(du), du, 0.0
    while abs(remaining) > 1e-15:
        slipping = (np.abs(z) >= LIMITS - 1e-12) & (np.sign(z) == d)
        k = 1 + gamma[~slipping].sum()
        room = np.where(slipping, np.inf, LIMITS - d * z)
        dx = remaining / k
        if abs(dx) > room.min():                       # next element reaches its limit
            dx = d * room.min()
        z = np.where(slipping, z, np.clip(z + dx, -LIMITS, LIMITS))
        remaining -= dx * k
        dx_total += dx
    return z, dx_total


def ramp(sign0, gamma, direction=-1):
    """Ramp from a saturated state z = sign0 * a; SUB records per microstep."""
    z = sign0 * LIMITS
    u, x = 0.0, -gamma @ z
    rec = [(u, x, z)]
    for _ in range(STEPS * SUB):
        z, dx = advance(z, direction * USTEP / SUB, gamma)
        u += direction * USTEP / SUB
        x += dx
        assert abs((u - x) - gamma @ z) < 1e-9
        rec.append((u, x, z))
    u, x = np.array([r[0] for r in rec]), np.array([r[1] for r in rec])
    z = np.array([r[2] for r in rec])
    travel = direction * u                              # commanded travel, >= 0
    e = direction * ((x - u) - (x[0] - u[0]))           # lag in ramp direction, from 0
    per_step = e[::SUB]
    slope = np.diff(per_step) / USTEP                   # de/ds per microstep
    return dict(u=u, x=x, z=z, s=travel, e=e, slope=slope, dir=direction, sign0=sign0)


def regime(slope):
    if slope < -.8:
        return 'stuck', 'stuck: steps go into the spring'
    if slope < -.1:
        return 'pre-sliding', 'pre-sliding: bristles bending'
    return 'sliding', 'sliding: block follows'


def spring_xy(x0, x1, y, n=9, amp=.22):
    xs = np.linspace(x0, x1, 2 * n + 2)
    ys = y + amp * np.r_[0, [(-1) ** k for k in range(2 * n)], 0]
    return xs, ys


def main():
    beta = json.loads((HERE / 'f0b_config.json').read_text())['beta_0_um']
    gamma = GAMMA0 * (beta / 2) / (GAMMA0 @ LIMITS)
    S1, S4 = ramp(+1, gamma), ramp(-1, gamma)
    assert abs(S1['e'][-1] + beta) < .02 and abs(S4['e']).max() < 1e-9
    OUT.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family': 'serif', 'font.serif': ['DejaVu Serif'],
                         'mathtext.fontset': 'dejavuserif', 'font.size': 10,
                         'axes.labelsize': 10, 'axes.titlesize': 11,
                         'axes.spines.top': False, 'axes.spines.right': False})

    # ------------------------------------------------------------------ figure
    fig = plt.figure(figsize=(10, 6.7), dpi=110, facecolor='white')
    gs = fig.add_gridspec(2, 2, left=.07, right=.97, bottom=.17, top=.84,
                          height_ratios=[1.25, .85], hspace=.5, wspace=.22,
                          width_ratios=[1, 1.05])
    mech = fig.add_subplot(gs[0, 0])
    lag = fig.add_subplot(gs[0, 1])
    stair = fig.add_subplot(gs[1, :])
    fig.text(.07, .955, 'S · slow ramp after a reversal vs in continuation', fontsize=16)
    heading = fig.text(.07, .898, '', fontsize=12, color=OCHRE)
    fig.text(.07, .072, 'Illustrative quasi-static Maxwell-slip model · placeholder slip limits '
             '· no F0b measurements', fontsize=9, color=INK)
    fig.text(.07, .036, f'Saturated lag Σγa = β₀/2 = {beta/2 + 1e-9:.3f} µm · (a) lag and travel share one '
             'scale, bristles not to scale · animation pace ≠ experiment time',
             fontsize=9, color=INK)

    # (a) mechanism
    mech.set(xlim=(0, 16), ylim=(-.9, 4.3), title='(a) Mechanism')
    mech.title.set_ha('left'); mech.title.set_x(0)
    mech.axis('off')
    ground_y, block_y0, block_y1, spring_y, rail_y = 0., .9, 2.3, 1.6, 3.4
    mech.add_patch(Rectangle((0, -.45), 16, .45, facecolor='white', edgecolor=GREY,
                             hatch='////', lw=0))
    mech.plot([0, 16], [ground_y] * 2, color=INK, lw=1)
    mech.plot([0, 16], [rail_y] * 2, color=GREY, lw=2.5, solid_capstyle='butt')
    ticks = [mech.plot([], [], color=INK, lw=.8)[0] for _ in range(12)]
    mech.text(15.9, -.7, 'ground ticks 1 µm apart', ha='right', fontsize=8, color=INK)
    block = mech.add_patch(Rectangle((0, block_y0), BLOCK_W, block_y1 - block_y0,
                                     facecolor='#eef1f2', edgecolor=DARK, lw=1.3, zorder=4))
    block_label = mech.text(0, (block_y0 + block_y1) / 2, 'carriage\nx', ha='center',
                            va='center', fontsize=8.5, color=DARK, zorder=5)
    handle, = mech.plot([], [], color=DARK, lw=4, solid_capstyle='butt', zorder=5)
    handle_label = mech.text(0, rail_y + .3, 'command u', ha='center', fontsize=8.5, color=DARK)
    spring, = mech.plot([], [], lw=1.6, zorder=3)
    spring_label = mech.text(16, 4.15, '', ha='right', fontsize=9, color=INK)
    bristles = [mech.plot([], [], lw=1.4, zorder=3)[0] for _ in LIMITS]
    feet = [mech.plot([], [], 'o', ms=3.5, zorder=3)[0] for _ in LIMITS]
    status = mech.text(0, 4.15, '', fontsize=10, fontweight='bold', color=DARK)

    # (b) lag vs commanded travel
    lag.set(title='(b) Lag vs commanded travel', xlabel='commanded travel since ramp start [µm]',
            ylabel='lag e (ramp direction) [µm]', xlim=(0, 20), ylim=(-2.05, .45))
    lag.title.set_ha('left'); lag.title.set_x(0)
    lag.set_axisbelow(True)
    lag.grid(color='#e2e6e8', linewidth=.6)
    lag.plot([0, 2.05], [0, -2.05], color=GREY, ls=':', lw=1)
    lag.text(1.95, -1.98, 'slope −1: carriage stuck', fontsize=8, color=INK, ha='left', va='bottom')
    lag.axhline(-beta, color=GREY, ls=':', lw=1)
    lag.text(19.7, -beta - .07, f'lag flipped by β₀ = {beta:.3f} µm', fontsize=8, color=INK,
             ha='right', va='top')
    s1_line, = lag.plot([], [], color=OCHRE, lw=2, label='S1 · reversal')
    s4_line, = lag.plot([], [], '--', color=OCHRE, lw=2, label='S4 · continuation')
    d_line, = lag.plot([], [], color=DARK, lw=2.6, label='Δ = S1 − S4')
    live, = lag.plot([], [], 'o', color=DARK, ms=6, zorder=6)
    lag.legend(handles=[s1_line, s4_line, d_line], loc='upper right', ncol=3, fontsize=8,
               frameon=True, handlelength=1.8, columnspacing=1)
    inset = lag.inset_axes([.42, .33, .55, .32])
    inset.set_xscale('log')
    inset.set(xlim=(.039, 20), ylim=(-2.05, .3))
    inset.xaxis.set_major_locator(matplotlib.ticker.FixedLocator([.04, .2, 1, 5, 20]))
    inset.xaxis.set_major_formatter(matplotlib.ticker.FixedFormatter(['0.04', '0.2', '1', '5', '20']))
    inset.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    inset.set_yticks([0, -1, -2])
    inset.tick_params(labelsize=7)
    inset.set_title('same, log travel', fontsize=8, loc='left', pad=2)
    inset.axhline(-beta, color=GREY, ls=':', lw=.8)
    in_s1, = inset.plot([], [], color=OCHRE, lw=1.6)
    in_s4, = inset.plot([], [], '--', color=OCHRE, lw=1.6)
    in_d, = inset.plot([], [], color=DARK, lw=2)
    # Regime boundaries from the S1 slope (used by the phase-1d annotations).
    steps_s = np.arange(1, STEPS + 1) * USTEP
    reg = np.array([regime(v)[0] for v in S1['slope']])
    b_pre = steps_s[np.argmax(reg != 'stuck')]
    b_slide = steps_s[np.argmax(reg == 'sliding')]
    spans = [inset.axvspan(lo, hi, color=c, alpha=.12, lw=0, visible=False)
             for lo, hi, c in ((.039, b_pre, DARK), (b_pre, b_slide, OCHRE), (b_slide, 20, GREY))]
    notes = [lag.text(*xy, t, fontsize=8.5, color=c, visible=False, ha=ha)
             for xy, t, c, ha in (((.35, -.35), 'deadband\n(slope ≈ −1)', DARK, 'left'),
                                  ((2.6, -1.25), 'pre-sliding\n(bristles bend)', OCHRE, 'left'),
                                  ((9.0, -1.62), 'saturated (lag flipped)', INK, 'left'))]

    # (c) first 30 microsteps on a time axis
    n30 = 30
    stair.set(title='(c) First 30 microsteps of the current ramp',
              xlabel='time since ramp start [s] (6.4 µsteps/s)',
              ylabel='travel in ramp direction [µm]', xlim=(0, n30 / RATE), ylim=(-.05, 1.3))
    stair.title.set_ha('left'); stair.title.set_x(0)
    stair.set_axisbelow(True)
    stair.grid(color='#e2e6e8', linewidth=.6)
    t30 = np.arange(n30 + 1) / RATE
    cmd_line, = stair.plot([], [], color=GREY, lw=1.2, drawstyle='steps-post', label='command u')
    ghost_line, = stair.plot([], [], color=OCHRE, lw=1.2, alpha=.35, drawstyle='steps-post',
                             label='carriage x · S1 (earlier)')
    car_line, = stair.plot([], [], color=OCHRE, lw=2, drawstyle='steps-post', label='carriage x')
    stair.legend(handles=[cmd_line, car_line], loc='upper left', fontsize=8, frameon=False)
    stair.text(.995, .04, f'1 µstep = {USTEP*1000:.1f} nm at 1/256', transform=stair.transAxes,
               ha='right', fontsize=8, color=INK)

    # ------------------------------------------------------------- drawing
    view = {'cam': 0.0}

    def draw_mech(u, x, z, d_motion, slipping, cam):
        """Draw one mechanism state; world positions in um relative to x0."""
        bx = BLOCK_SCREEN + SCALE * x - cam
        length = L0 + SCALE * (x - u)
        hx = bx - length
        block.set_x(bx)
        block_label.set_x(bx + BLOCK_W / 2)
        handle.set_data([hx, hx], [spring_y - .45, rail_y])
        handle_label.set_x(hx)
        load = (x - u) / (beta / 2)
        colour = GREY if abs(load) < .25 else OCHRE
        spring.set_data(*spring_xy(hx, bx, spring_y))
        spring.set_color(colour)
        spring_label.set_text('drive spring: ' + ('relaxed' if abs(load) < .25 else
                                                  'pulling ↓' if load > 0 else 'pushing ↑'))
        for i, (line, foot) in enumerate(zip(bristles, feet)):
            tip = bx + .25 + i * (BLOCK_W - .5) / (len(LIMITS) - 1)
            base = tip - .55 * z[i] / LIMITS[i]
            c = OCHRE if slipping[i] else GREY
            line.set_data([base, tip], [ground_y, block_y0])
            line.set_color(c)
            foot.set_data([base], [ground_y])
            foot.set_color(c)
        first = np.floor(cam / SCALE) - 4
        for k, tick in enumerate(ticks):
            gx = BLOCK_SCREEN + SCALE * (first + k) - cam
            tick.set_data([gx, gx], [ground_y, ground_y - .22])

    def slipping_of(z, d):
        return (np.abs(z) >= LIMITS - 1e-9) & (np.sign(z) == d)

    def ramp_cam(r, i):
        """Static view until the carriage nears the left edge, then pan with it."""
        bx = BLOCK_SCREEN + SCALE * r['x'][i]
        return min(0.0, bx - PAN_EDGE) if bx < PAN_EDGE else 0.0

    frames = []

    def grab(n=1):
        fig.canvas.draw()
        frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
        frames.extend([frame] + [frame.copy() for _ in range(n - 1)])

    def show_ramp_curves(r, i, line, in_line):
        line.set_data(r['s'][:i + 1], r['e'][:i + 1])
        keep = r['s'][:i + 1] > 0
        in_line.set_data(r['s'][:i + 1][keep], r['e'][:i + 1][keep])
        live.set_data([r['s'][i]], [r['e'][i]])
        k = min(i // SUB, n30)
        cmd_line.set_data(t30[:k + 1], r['s'][::SUB][:k + 1])
        car_line.set_data(t30[:k + 1], (r['dir'] * (r['x'] - r['x'][0]))[::SUB][:k + 1])

    def wipe(r, name, n=20):
        """Schematic saturated approach to x0 from the side opposite the ramp."""
        sign0 = r['sign0']
        z = sign0 * LIMITS
        for line in (cmd_line, car_line):
            line.set_data([], [])
        live.set_data([], [])
        for j in range(n + 1):
            u = -sign0 * WIPE_UM * (1 - j / n)
            x = u - gamma @ z
            cam = SCALE * (x - r['x'][0])            # follow the block, stop at x0
            draw_mech(u, x, z, sign0, slipping_of(z, sign0), cam)
            side = 'below (last motion up)' if sign0 > 0 else 'above (last motion down)'
            status.set_text(f'memory wipe from {side}')
            heading.set_text(f'{name} · preparation · wipe to x₀ from {side.split(" (")[0]}')
            grab()

    def run(r, name, kind):
        idx = list(range(1, n30 * SUB + 1)) + list(range((n30 + 2) * SUB, STEPS * SUB + 1, 8 * SUB))
        for i in idx:
            step = int(np.ceil(i / SUB))
            reg_key, reg_text = regime(r['slope'][step - 1])
            draw_mech(r['u'][i], r['x'][i], r['z'][i], r['dir'],
                      slipping_of(r['z'][i], r['dir']), ramp_cam(r, i))
            status.set_text(reg_text)
            heading.set_text(f'{name} · {kind} · step {step}/{STEPS} · carriage: {reg_key}')
            show_ramp_curves(r, i, s1_line if name == 'S1' else s4_line,
                             in_s1 if name == 'S1' else in_s4)
            grab()

    # --------------------------------------------------------------- phases
    # 1: S1 reversal
    wipe(S1, 'S1')
    run(S1, 'S1', 'reversal')
    for sp in spans:
        sp.set_visible(True)
    for note in notes:
        note.set_visible(True)
    heading.set_text('S1 · reversal · deadband → pre-sliding → saturated')
    grab(28)
    for sp in spans:
        sp.set_visible(False)
    for note in notes:
        note.set_visible(False)
    # 2: S4 continuation, with S1 kept for comparison
    ghost_line.set_data(t30, (S1['dir'] * (S1['x'] - S1['x'][0]))[::SUB][:n30 + 1])
    stair.legend(handles=[cmd_line, car_line, ghost_line], loc='upper left', fontsize=8,
                 frameon=False)
    wipe(S4, 'S4')
    run(S4, 'S4', 'continuation')
    grab(10)
    # 3: overlay and difference
    live.set_data([], [])
    delta = S1['e'] - S4['e']
    for j in range(1, 21):
        i = round(j / 20 * STEPS * SUB)
        d_line.set_data(S1['s'][:i + 1], delta[:i + 1])
        keep = S1['s'][:i + 1] > 0
        in_d.set_data(S1['s'][:i + 1][keep], delta[:i + 1][keep])
        heading.set_text('Δ = S1 − S4 over the same commanded positions')
        grab()
    heading.set_text('Δ = S1 − S4: same positions, so geometry and microstep error cancel; '
                     'Δ is the friction transient')
    status.set_text('')
    grab(40)

    path = OUT / 'f0b_S_slow_ramps.gif'
    title_frames = title_card(beta, fig.get_size_inches(), fig.dpi)
    save_gif(title_frames + frames, path)
    total = len(title_frames) + len(frames)
    print(f'{path}\n{total} frames; {total/FPS:.1f} s; {path.stat().st_size/1e6:.2f} MB\n'
          f'S1 regimes: stuck < {b_pre:.3f} um, pre-sliding < {b_slide:.3f} um, '
          f'final lag {S1["e"][-1]:+.3f} um; initial slope {S1["slope"][0]:+.3f}')
    plt.close(fig)


def title_card(beta, size, dpi, n=24):
    fig = plt.figure(figsize=size, dpi=dpi, facecolor='white')
    fig.text(.5, .62, 'Slow ramp after a reversal vs in continuation', ha='center', fontsize=20)
    fig.text(.5, .5, 'S1: wipe from below, then ramp down 20 µm  →  reversal\n'
             'S4: wipe from above, then ramp down 20 µm  →  continuation',
             ha='center', va='center', fontsize=12.5, color=OCHRE, linespacing=1.8)
    fig.text(.5, .34, f'512 microsteps at 1/256 · 0.25 µm/s · illustrative friction model '
             f'(saturation β₀/2 = {beta/2:.3f} µm)', ha='center', fontsize=10, color=INK)
    fig.canvas.draw()
    frame = Image.fromarray(np.asarray(fig.canvas.buffer_rgba())[:, :, :3].copy())
    plt.close(fig)
    return [frame] + [frame.copy() for _ in range(n - 1)]


if __name__ == '__main__':
    main()
