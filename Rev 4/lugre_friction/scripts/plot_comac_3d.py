#!/usr/bin/env python3
'''Create 3D counterparts of the Global MAC and COMAC result figures.

The figures are rendered directly from the saved NPZ result arrays. This
keeps every bar at the same numerical value as the existing 2D figures.
Slots rejected by the Split Method are retained as fully opaque black bars.
'''

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.cm import ScalarMappable
from matplotlib.colors import LinearSegmentedColormap, Normalize


LUGRE_DIR = Path(__file__).resolve().parent.parent
COMAC_DIR = LUGRE_DIR / 'rendered_assets' / 'temp' / 'Co-MAC'
SPLIT_DIR = COMAC_DIR / 'Split Method'
KINEMATIC_DIR = COMAC_DIR / 'Kinematically Scaled'

PALETTE = ['#454040', '#605B51', '#D8D365', '#E6F082']
BEST_MATCH_COLOR = '#48A111'
THRESHOLD_COLOR = '#605B51'
FLAGGED_ROT = {(1, 0)}
FLAGGED_TRANS = {
    (1, 0),
    (2, 4), (2, 5),
    (3, 4), (3, 5),
    (4, 4), (4, 5),
    (5, 4), (5, 5),
}

CMAP = LinearSegmentedColormap.from_list('comac_palette', PALETTE)
NORM = Normalize(vmin=0.0, vmax=1.0)


def _style_3d_axis(ax) -> None:
    '''Use a light plotting box that leaves the colored bars dominant.'''
    ax.set_facecolor('white')
    for axis in (ax.xaxis, ax.yaxis, ax.zaxis):
        axis.pane.set_facecolor((1.0, 1.0, 1.0, 1.0))
        axis.pane.set_edgecolor('#b8b8b8')
        axis._axinfo['grid']['color'] = (0.78, 0.78, 0.78, 0.65)
        axis._axinfo['grid']['linewidth'] = 0.7
    ax.tick_params(axis='both', which='major', labelsize=7, pad=0)
    ax.tick_params(axis='z', which='major', labelsize=7, pad=1)
    ax.set_zlim(0.0, 1.05)
    ax.set_zticks(np.arange(0.0, 1.01, 0.2))
    ax.view_init(elev=28, azim=-56)


def _bar3d_group(ax, x, y, values, colors, edges, width=0.72, depth=0.72) -> None:
    '''Draw bars separately so every cuboid can have its own edge color.'''
    for xi, yi, value, color, edge in zip(x, y, values, colors, edges):
        ax.bar3d(
            xi - width / 2,
            yi - depth / 2,
            0.0,
            width,
            depth,
            float(value),
            color=color,
            edgecolor=edge,
            linewidth=1.35 if edge == BEST_MATCH_COLOR else 0.55,
            alpha=1.0,
            shade=True,
        )


def draw_mac_matrix_3d(ax, matrix: np.ndarray, title: str, flagged=frozenset()) -> None:
    '''Render one Set-A x Set-B MAC matrix as height-coded 3D cuboids.'''
    n_rows, n_cols = matrix.shape
    rows, cols = np.indices(matrix.shape)
    x = cols.ravel()
    y = rows.ravel()
    values = matrix.ravel()
    best = {(i, int(np.argmax(matrix[i]))) for i in range(n_rows)}

    colors = []
    edges = []
    for i, j, value in zip(y, x, values):
        ij = (int(i), int(j))
        colors.append('#000000' if ij in flagged else CMAP(NORM(value)))
        edges.append(BEST_MATCH_COLOR if ij in best else '#f5f5f5')

    _bar3d_group(ax, x, y, values, colors, edges)

    # Labels sit above the true bar height; bar geometry remains exact.
    for xi, yi, value in zip(x, y, values):
        is_flagged = (int(yi), int(xi)) in flagged
        ax.text(
            xi,
            yi,
            min(float(value) + 0.025, 1.035),
            f'{value:.2f}',
            ha='center',
            va='bottom',
            fontsize=5.5,
            color=PALETTE[3] if is_flagged else PALETTE[0],
            zorder=20,
        )

    ax.set_xticks(np.arange(n_cols))
    ax.set_xticklabels([f'B-mode{i + 1}' for i in range(n_cols)], rotation=-18, ha='left')
    ax.set_yticks(np.arange(n_rows))
    ax.set_yticklabels([f'A-mode{i + 1}' for i in range(n_rows)], rotation=18, ha='right')
    ax.set_xlabel('Set B', fontsize=8, labelpad=5)
    ax.set_ylabel('Set A', fontsize=8, labelpad=6)
    ax.set_title(title, fontsize=10, pad=10)
    ax.set_xlim(-0.55, n_cols - 0.25)
    ax.set_ylim(n_rows - 0.25, -0.55)
    _style_3d_axis(ax)
    ax.set_zlabel('Global MAC', fontsize=8, labelpad=4)


def draw_comac_bars_3d(ax, labels: list[str], values: np.ndarray, title: str) -> None:
    '''Render a one-row COMAC vector with the same value-mapped palette.'''
    x = np.arange(len(labels))
    y = np.zeros(len(labels))
    colors = [CMAP(NORM(v)) for v in values]
    edges = ['#f5f5f5'] * len(labels)
    _bar3d_group(ax, x, y, values, colors, edges, width=0.62, depth=0.68)

    for xi, value in zip(x, values):
        ax.text(
            xi,
            0.0,
            min(float(value) + 0.025, 1.035),
            f'{value:.2f}',
            ha='center',
            va='bottom',
            fontsize=7,
            color=PALETTE[0],
            zorder=20,
        )

    # Lift the dashed reference into the x-z plane behind the bars.
    ax.plot(
        [-0.45, len(labels) - 0.55],
        [0.44, 0.44],
        [0.9, 0.9],
        color=THRESHOLD_COLOR,
        linestyle='--',
        linewidth=1.4,
    )
    ax.text(len(labels) - 0.5, 0.44, 0.92, '0.9 threshold', fontsize=6.5,
            color=THRESHOLD_COLOR, ha='right')

    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=-18, ha='left')
    ax.set_yticks([])
    ax.set_xlim(-0.5, len(labels) - 0.35)
    ax.set_ylim(-0.5, 0.55)
    ax.set_title(title, fontsize=10, pad=10)
    _style_3d_axis(ax)
    ax.set_zlabel('COMAC', fontsize=8, labelpad=4)


def save_split_mac(data) -> Path:
    fig = plt.figure(figsize=(15.0, 7.4), facecolor='white')
    axes = [fig.add_subplot(1, 2, i + 1, projection='3d') for i in range(2)]
    draw_mac_matrix_3d(axes[0], data['mac_rot'], 'Torsional Global MAC', FLAGGED_ROT)
    draw_mac_matrix_3d(
        axes[1], data['mac_trans'], 'Linear/Translational Global MAC', FLAGGED_TRANS
    )
    fig.suptitle(
        'Global MAC Matrices (Rotational vs. Translational Sub-Vectors) -- 3D',
        fontsize=14,
        y=0.97,
    )
    sm = ScalarMappable(norm=NORM, cmap=CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=axes, fraction=0.025, pad=0.02, shrink=0.72)
    cbar.set_label('Global MAC', fontsize=9)
    fig.text(
        0.5,
        0.025,
        'Solid black cuboids: excluded / unreliable slots; green edges: row-wise best match',
        ha='center',
        fontsize=8.5,
    )
    out = SPLIT_DIR / 'comac_mac_matrices_3d.png'
    fig.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out


def save_kinematic_mac(data) -> Path:
    fig = plt.figure(figsize=(9.0, 7.4), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    draw_mac_matrix_3d(ax, data['mac_scaled'], 'Kinematically Scaled Global MAC (6-DOF)')
    fig.suptitle('Global MAC Matrix -- Kinematically Scaled (6-DOF, unified) -- 3D',
                 fontsize=13, y=0.97)
    sm = ScalarMappable(norm=NORM, cmap=CMAP)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.04, shrink=0.72)
    cbar.set_label('Global MAC', fontsize=9)
    fig.text(0.5, 0.025, 'Green edges: row-wise best match', ha='center', fontsize=8.5)
    out = KINEMATIC_DIR / 'comac_kinematic_mac_matrix_3d.png'
    fig.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out


def save_split_comac(data) -> Path:
    fig = plt.figure(figsize=(12.0, 6.7), facecolor='white')
    grid = fig.add_gridspec(1, 2, width_ratios=[4, 2])
    axes = [fig.add_subplot(grid[0, i], projection='3d') for i in range(2)]
    draw_comac_bars_3d(
        axes[0], [str(v) for v in data['rot_labels']], data['comac_rot'][:, 3],
        'Rotational COMAC',
    )
    draw_comac_bars_3d(
        axes[1], [str(v) for v in data['trans_labels']], data['comac_trans'][:, 3],
        'Translational COMAC',
    )
    fig.suptitle('Coordinate Modal Assurance Criterion (COMAC) per Degree of Freedom -- 3D',
                 fontsize=13, y=0.97)
    out = SPLIT_DIR / 'comac_stage3_3d.png'
    fig.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out


def save_kinematic_comac(data) -> Path:
    fig = plt.figure(figsize=(9.0, 6.7), facecolor='white')
    ax = fig.add_subplot(111, projection='3d')
    draw_comac_bars_3d(
        ax, [str(v) for v in data['state_labels']], data['comac'][:, 3],
        'Unified 6-DOF COMAC (Kinematically Scaled)',
    )
    fig.suptitle('Coordinate Modal Assurance Criterion (COMAC) -- Kinematically Scaled -- 3D',
                 fontsize=13, y=0.97)
    out = KINEMATIC_DIR / 'comac_kinematic_comac_3d.png'
    fig.savefig(out, dpi=180, bbox_inches='tight', facecolor='white')
    plt.close(fig)
    return out


def main() -> None:
    split_npz = SPLIT_DIR / 'npz' / 'comac_mode_extraction_data.npz'
    kinematic_npz = KINEMATIC_DIR / 'npz' / 'comac_kinematic_scaled_data.npz'
    with np.load(split_npz) as split_data, np.load(kinematic_npz) as kinematic_data:
        outputs = [
            save_split_mac(split_data),
            save_split_comac(split_data),
            save_kinematic_mac(kinematic_data),
            save_kinematic_comac(kinematic_data),
        ]
    for output in outputs:
        print(f'Wrote {output}')


if __name__ == '__main__':
    main()
