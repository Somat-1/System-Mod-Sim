#!/usr/bin/env python3
'''Shared figure style, so every plot in the campaign reads as one set.

Colours are the validated categorical palette: hues are assigned by slot in a
fixed order and never cycled. Scatter-type figures use at most the first three
slots, which are the ones that clear the colour-vision separation floors on
every pair rather than only on adjacent pairs.
'''

from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt

SURFACE = '#fcfcfb'
INK = '#0b0b0b'
INK_SOFT = '#52514e'
INK_MUTED = '#86857f'
GRID = '#e4e3df'

# Fixed categorical order; index by meaning, never by loop position.
SERIES = ['#2a78d6', '#eb6834', '#1baf7a', '#eda100',
          '#e87ba4', '#008300', '#6250d6', '#e34948']

# One hue, light -> dark, for magnitude.
SEQUENTIAL = ['#cde2fb', '#9ec5f4', '#6da7ec', '#3987e5',
              '#256abf', '#184f95', '#0d366b']

FORWARD, REVERSE = SERIES[0], SERIES[1]
STATUS_BAD = '#e34948'

PHASE_COLOR = {
    'gain': SERIES[0],
    'staircase': SERIES[1],
    'settling': SERIES[2],
    'repeatability': SERIES[3],
    'drift': INK_MUTED,
}


def use():
    '''Apply the campaign style to the global matplotlib state.'''
    mpl.rcParams.update({
        'figure.facecolor': SURFACE,
        'axes.facecolor': SURFACE,
        'savefig.facecolor': SURFACE,
        'axes.edgecolor': GRID,
        'axes.labelcolor': INK_SOFT,
        'axes.titlecolor': INK,
        'axes.titlesize': 11,
        'axes.titleweight': 'bold',
        'axes.titlelocation': 'left',
        'axes.labelsize': 9,
        'axes.linewidth': 0.8,
        'axes.grid': True,
        'axes.axisbelow': True,
        'grid.color': GRID,
        'grid.linewidth': 0.6,
        'xtick.color': INK_MUTED,
        'ytick.color': INK_MUTED,
        'xtick.labelsize': 8,
        'ytick.labelsize': 8,
        'xtick.direction': 'out',
        'ytick.direction': 'out',
        'legend.frameon': False,
        'legend.fontsize': 8,
        'font.size': 9,
        'lines.linewidth': 1.4,
        'lines.solid_capstyle': 'round',
        'figure.dpi': 110,
        'savefig.dpi': 200,
    })


def despine(ax, keep=('left', 'bottom')):
    for side in ('top', 'right', 'left', 'bottom'):
        ax.spines[side].set_visible(side in keep)


def note(ax, text, loc='upper left'):
    '''A short annotation in muted ink, never in a series colour.'''
    xy = {'upper left': (0.012, 0.975), 'upper right': (0.988, 0.975),
          'lower left': (0.012, 0.03), 'lower right': (0.988, 0.03)}[loc]
    ax.text(*xy, text, transform=ax.transAxes, fontsize=7.5, color=INK_MUTED,
            ha='right' if 'right' in loc else 'left',
            va='bottom' if 'lower' in loc else 'top', linespacing=1.5)


def save(fig, out_dir, stem, dpi=None):
    '''Write the figure as a PNG, at print resolution in thesis mode.

    Every figure here is sized in millimetres for the page it will sit on, so
    a 600 dpi raster at that size is indistinguishable from the vector it was
    drawn as, and one file per figure is one thing to keep track of instead
    of two.
    '''
    out_dir.mkdir(parents=True, exist_ok=True)
    png = out_dir / f'{stem}.png'
    fig.savefig(png, dpi=dpi or (600 if THESIS else None))
    plt.close(fig)
    print(f'  wrote {png.name}')
    return png


# --------------------------------------------------------------------------
# thesis mode
# --------------------------------------------------------------------------
#
# The campaign's figures are read in two places: on a screen while the work
# is going on, where a headline on the figure is the fastest way to see what
# it says, and in the thesis, where a headline on the figure is wrong - the
# finding belongs in the caption's first sentence and the figure carries only
# panel letters. One switch covers both, so no figure has to be rewritten to
# move between them.

THESIS = True

# Print widths, in inches. A thesis text block is 160 mm; a half-width figure
# that sits two to a row is 80 mm.
FULL = 160 / 25.4
HALF = 80 / 25.4
TWO_THIRDS = 110 / 25.4


def use_thesis():
    """Draw at print size, in the body font, with nothing decorative."""
    plt.style.use(str(Path(__file__).resolve().parent / 'thesis.mplstyle'))


def headline(fig, text, **kwargs):
    """The finding, shown on screen and omitted from the thesis.

    In thesis mode this does nothing: the same sentence is the first line of
    the caption, and having it twice is how a figure comes to look like a
    slide.
    """
    if THESIS:
        return None
    return fig.suptitle(text, fontsize=11, fontweight='bold', x=0.008,
                        ha='left', **kwargs)


def panel(ax, letter, dx=-0.02, dy=1.0):
    """The (a), (b) that the text refers to, outside the axes at top left."""
    return ax.annotate(f'({letter})', (dx, dy), xycoords='axes fraction',
                       ha='right', va='bottom', fontweight='bold')


def corner(ax, text, loc='upper left', **kwargs):
    """A small label inside the panel, for what would otherwise be a title."""
    xy = {'upper left': (0.025, 0.96), 'upper right': (0.975, 0.96),
          'lower left': (0.025, 0.04), 'lower right': (0.975, 0.04)}[loc]
    return ax.annotate(text, xy, xycoords='axes fraction',
                       ha='right' if 'right' in loc else 'left',
                       va='bottom' if 'lower' in loc else 'top', **kwargs)
