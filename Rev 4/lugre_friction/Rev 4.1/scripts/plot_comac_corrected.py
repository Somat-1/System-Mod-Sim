#!/usr/bin/env python3
'''Generate corrected Global MAC and COMAC results in a separate temp folder.

Corrections:
1. Use the common physical mass matrix as the MAC metric, avoiding dependence
   on the arbitrary coordinate scaling used by an ordinary Euclidean MAC.
2. Obtain a one-to-one mode map with the Hungarian assignment.
3. Apply one mass-metric polarity sign to every unique paired mode.
4. Evaluate standard signed-sum COMAC on the raw coordinates. Per-row scaling
   is deliberately omitted because it cancels from COMAC for fixed pairs.
'''

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
from scipy.optimize import linear_sum_assignment


LUGRE_DIR = Path(__file__).resolve().parent.parent
COMAC_DIR = LUGRE_DIR / 'rendered_assets' / 'temp' / 'Co-MAC'
SPLIT_DIR = COMAC_DIR / 'Split Method'
KINEMATIC_DIR = COMAC_DIR / 'Kinematically Scaled'
OUT_DIR = KINEMATIC_DIR / 'temp_corrected'
NPZ_DIR = OUT_DIR / 'npz'

PALETTE = ['#454040', '#605B51', '#D8D365', '#E6F082']
BEST_MATCH_COLOR = '#48A111'
THRESHOLD_COLOR = '#605B51'
CMAP = LinearSegmentedColormap.from_list('comac_palette', PALETTE)


def mass_weighted_mac(phi_a: np.ndarray, phi_b: np.ndarray,
                      mass: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    '''Return the MAC matrix and signed cross-inner-product matrix.'''
    cross = phi_a.T @ mass @ phi_b
    norm_a = np.sum(phi_a * (mass @ phi_a), axis=0)
    norm_b = np.sum(phi_b * (mass @ phi_b), axis=0)
    mac = np.abs(cross) ** 2 / np.outer(norm_a, norm_b)
    return mac, cross


def signed_comac(phi_a: np.ndarray, phi_b: np.ndarray,
                 pairs: np.ndarray, signs: np.ndarray) -> np.ndarray:
    '''Return columns cross-sum, A square-sum, B square-sum, and COMAC.'''
    a = phi_a[:, pairs[:, 0]]
    b = phi_b[:, pairs[:, 1]] * signs[None, :]
    cross_sum = np.sum(a * b, axis=1)
    den_a = np.sum(a * a, axis=1)
    den_b = np.sum(b * b, axis=1)
    values = cross_sum ** 2 / (den_a * den_b)
    return np.column_stack((cross_sum, den_a, den_b, values))


def absolute_product_comac(phi_a: np.ndarray, phi_b: np.ndarray,
                           pairs: np.ndarray) -> np.ndarray:
    '''Sensitivity variant with magnitudes inside the modal cross-sum.'''
    a = phi_a[:, pairs[:, 0]]
    b = phi_b[:, pairs[:, 1]]
    return np.sum(np.abs(a * b), axis=1) ** 2 / (
        np.sum(a * a, axis=1) * np.sum(b * b, axis=1)
    )


def draw_mac(mac: np.ndarray, pairs: np.ndarray, output: Path) -> None:
    n = mac.shape[0]
    fig, ax = plt.subplots(figsize=(7.2, 6.1))
    image = ax.imshow(mac, cmap=CMAP, vmin=0.0, vmax=1.0, aspect='auto')

    ax.set_xticks(np.arange(-0.5, n, 1.0), minor=True)
    ax.set_yticks(np.arange(-0.5, n, 1.0), minor=True)
    ax.grid(which='minor', color='white', linewidth=1.6)
    ax.tick_params(which='minor', length=0)
    ax.tick_params(which='major', length=0)
    ax.set_xticks(np.arange(n))
    ax.set_xticklabels([f'B-mode{i + 1}' for i in range(n)], rotation=45, ha='left', fontsize=8)
    ax.xaxis.tick_top()
    ax.set_yticks(np.arange(n))
    ax.set_yticklabels([f'A-mode{i + 1}' for i in range(n)], fontsize=8)

    for i in range(n):
        for j in range(n):
            value = mac[i, j]
            color = PALETTE[0] if value > 0.6 else PALETTE[3]
            ax.text(j, i, f'{value:.2f}', ha='center', va='center', fontsize=7.5, color=color)

    for i, j in pairs:
        ax.add_patch(Rectangle(
            (j - 0.5, i - 0.5), 1.0, 1.0, fill=False,
            edgecolor=BEST_MATCH_COLOR, linewidth=2.5, zorder=5, clip_on=False,
        ))

    ax.text(0.5, -0.08, 'Common-mass metric; green border = optimal bijective assignment',
            transform=ax.transAxes, ha='center', va='top', fontsize=9)
    fig.suptitle('Coordinate-Invariant Mass-Weighted Global MAC', fontsize=12)
    fig.tight_layout(rect=[0.0, 0.07, 0.88, 0.93])
    cbar_ax = fig.add_axes([0.90, 0.16, 0.03, 0.66])
    fig.colorbar(image, cax=cbar_ax, label='Mass-weighted MAC')
    fig.savefig(output, dpi=180, facecolor='white')
    plt.close(fig)


def draw_comac(labels: list[str], values: np.ndarray, output: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.5))
    x = np.arange(len(labels))
    colors = [CMAP(value) for value in values]
    ax.bar(x, values, color=colors, width=0.62, zorder=3)
    ax.axhline(0.9, color=THRESHOLD_COLOR, linewidth=1.3, linestyle='--', zorder=2)

    for xi, value in zip(x, values):
        va = 'top' if value >= 0.9 else 'bottom'
        offset = -0.02 if value >= 0.9 else 0.02
        ax.text(xi, value + offset, f'{value:.3f}', ha='center', va=va,
                fontsize=9, color=PALETTE[0])

    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=9)
    ax.set_ylim(0.0, 1.05)
    ax.set_yticks(np.arange(0.0, 1.01, 0.2))
    ax.set_ylabel('COMAC', fontsize=9.5)
    ax.grid(axis='y', color='#cccccc', linewidth=0.8, zorder=0)
    ax.text(0.5, -0.12, 'Six unique mass-MAC-paired modes; standard signed-sum COMAC',
            transform=ax.transAxes, ha='center', va='top', fontsize=9)
    fig.suptitle('Corrected Coordinate Modal Assurance Criterion', fontsize=12)
    fig.tight_layout(rect=[0.0, 0.07, 1.0, 0.93])
    fig.savefig(output, dpi=180, facecolor='white')
    plt.close(fig)


def markdown_table(headers: list[str], rows: list[list[str]]) -> str:
    head = '| ' + ' | '.join(headers) + ' |\n'
    sep = '|' + '|'.join(['---'] * len(headers)) + '|\n'
    body = ''.join('| ' + ' | '.join(row) + ' |\n' for row in rows)
    return head + sep + body


def write_report(output: Path, labels: list[str], mac: np.ndarray,
                 pairs: np.ndarray, signs: np.ndarray, comac: np.ndarray,
                 comac_abs: np.ndarray, freq_a: np.ndarray, freq_b: np.ndarray,
                 r_value: float, fraction_a: np.ndarray, fraction_b: np.ndarray,
                 original_comac: np.ndarray) -> None:
    mode_rows = []
    for (a, b), sign in zip(pairs, signs):
        mode_rows.append([
            f'A-mode{a + 1}', f'{freq_a[a]:.3f}',
            f'B-mode{b + 1}', f'{freq_b[b]:.3f}',
            f'{mac[a, b]:.4f}', '+' if sign > 0 else '-',
        ])

    mac_rows = []
    for i in range(mac.shape[0]):
        mac_rows.append([f'A-mode{i + 1}'] + [f'{value:.4f}' for value in mac[i]])

    comac_rows = []
    for label, row, abs_value, old_value in zip(labels, comac, comac_abs, original_comac):
        comac_rows.append([
            label, f'{row[0]:.6e}', f'{row[1]:.6e}', f'{row[2]:.6e}',
            f'{row[3]:.6f}', f'{abs_value:.6f}', f'{old_value:.6f}',
        ])

    fraction_rows = []
    for i in range(len(fraction_a)):
        fraction_rows.append([
            f'mode{i + 1}', f'{fraction_a[i]:.6f}', f'{fraction_b[i]:.6f}',
        ])

    text = '''# Corrected Kinematically Scaled Co-MAC Audit

This directory is a separate corrected result set. The original Kinematically
Scaled files are unchanged.

## Corrections implemented

1. **Coordinate-invariant MAC metric.** The two systems have exactly the same
   physical mass matrix. Global mode correlation therefore uses the common-mass
   inner product instead of Euclidean MAC after an arbitrary coordinate scaling:

   $$MAC_M(a,b)=\\frac{|\\phi_{A,a}^T M\\phi_{B,b}|^2}
   {(\\phi_{A,a}^T M\\phi_{A,a})(\\phi_{B,b}^T M\\phi_{B,b})}.$$

   This is invariant to a consistent coordinate transformation and does not
   claim that equal units imply equal DOF weighting.

2. **Bijective pairing.** `linear_sum_assignment(..., maximize=True)` selects
   one unique Set-B mode for every Set-A mode. No mode column is reused.

3. **Consistent polarity.** Each unique pair receives one sign from
   $\\operatorname{sign}(\\phi_A^T M\\phi_B)$. The standard signed-sum COMAC is
   then evaluated. The per-term-absolute result is reported only as a sensitivity
   comparison, not as the primary definition.

4. **COMAC is evaluated on raw coordinates.** For fixed pairs and signs, scaling
   a row in both datasets cancels identically from COMAC. The lead-screw scale is
   relevant to the superseded Euclidean pairing method, not to the row-wise
   COMAC normalization itself.

## Mass-weighted Global MAC

'''
    text += markdown_table(
        ['Set A / Set B'] + [f'B-mode{i + 1}' for i in range(mac.shape[1])],
        mac_rows,
    )
    text += '\nSee `comac_corrected_mac_matrix.png`. Green borders show the optimal bijection.\n\n'
    text += '## Optimal bijective mode pairs\n\n'
    text += markdown_table(
        ['Set A', 'f_A [Hz]', 'Set B', 'f_B [Hz]', 'MAC_M', 'sign'], mode_rows,
    )
    text += f'\nTotal assigned mass-weighted MAC = `{sum(mac[a, b] for a, b in pairs):.6f}`.\n\n'
    text += '## Corrected COMAC\n\n'
    text += markdown_table(
        ['DOF', 'sum(phiA phiB)', 'sum(phiA^2)', 'sum(phiB^2)',
         'signed COMAC', 'absolute-product check', 'original duplicate'],
        comac_rows,
    )
    text += '''
The corrected result does not support interpreting the original low
`theta_sb = 0.267` value as a direct friction signature. With unique paired
modes it becomes approximately `0.782`. The result should still be interpreted
as coordinate correlation conditional on this mode map, rather than as a direct
damage or friction-localization measurement.

## Scaling audit

The superseded kinematic method used
'''
    text += f'`R = {r_value:.12e} m/rad` and `R^2 = {r_value ** 2:.12e}`. '
    text += '''Its fraction of scaled Euclidean modal norm in the two
translational rows was:\n\n'''
    text += markdown_table(['Mode', 'Set A translational fraction',
                            'Set B translational fraction'], fraction_rows)
    text += '''
These fractions show that equalizing units did not equalize modal weighting.
The corrected Global MAC therefore uses the common physical mass metric.

## Split Method sign clarification

The Split Method already applied signs `(+,-,+)` to its three translational
pairs. Its near-zero `x_s` and `x_n` results were caused by the selected
three-pair amplitude correspondence and the small sample size, not by unresolved
eigenvector signs. No claim to the contrary is carried into this result set.
'''
    output.write_text(text, encoding='utf-8')


def main() -> None:
    split_path = SPLIT_DIR / 'npz' / 'comac_mode_extraction_data.npz'
    scaled_path = KINEMATIC_DIR / 'npz' / 'comac_kinematic_scaled_data.npz'

    with np.load(split_path) as split, np.load(scaled_path) as scaled:
        phi_a = split['phi0']
        phi_b = split['phi1']
        mass_a = split['M0']
        mass_b = split['M1']
        if not np.array_equal(mass_a, mass_b):
            raise ValueError('A common mass metric requires identical Set-A and Set-B mass matrices.')

        mac, cross = mass_weighted_mac(phi_a, phi_b, mass_a)
        rows, cols = linear_sum_assignment(mac, maximize=True)
        pairs = np.column_stack((rows, cols))
        signs = np.sign(cross[rows, cols])
        signs[signs == 0.0] = 1.0

        comac = signed_comac(phi_a, phi_b, pairs, signs)
        comac_abs = absolute_product_comac(phi_a, phi_b, pairs)
        labels = [str(value) for value in split['state_labels']]

        phi_a_scaled = scaled['phi0_scaled']
        phi_b_scaled = scaled['phi1_scaled']
        fraction_a = np.sum(phi_a_scaled[4:, :] ** 2, axis=0) / np.sum(phi_a_scaled ** 2, axis=0)
        fraction_b = np.sum(phi_b_scaled[4:, :] ** 2, axis=0) / np.sum(phi_b_scaled ** 2, axis=0)

        OUT_DIR.mkdir(parents=True, exist_ok=True)
        NPZ_DIR.mkdir(parents=True, exist_ok=True)

        mac_figure = OUT_DIR / 'comac_corrected_mac_matrix.png'
        comac_figure = OUT_DIR / 'comac_corrected_comac.png'
        report = OUT_DIR / 'comac_corrected_report.md'

        draw_mac(mac, pairs, mac_figure)
        draw_comac(labels, comac[:, 3], comac_figure)
        write_report(
            report, labels, mac, pairs, signs, comac, comac_abs,
            split['freq0_hz'], split['freq1_hz'], float(scaled['R']),
            fraction_a, fraction_b, scaled['comac'][:, 3],
        )

        np.savez(
            NPZ_DIR / 'comac_corrected_data.npz',
            state_labels=np.array(labels), mass=mass_a,
            phi0=phi_a, phi1=phi_b,
            freq0_hz=split['freq0_hz'], freq1_hz=split['freq1_hz'],
            mac_mass=mac, cross_mass=cross,
            pairs=pairs, signs=signs,
            comac=comac, comac_absolute_product=comac_abs,
            original_comac=scaled['comac'][:, 3],
            scaled_trans_fraction_a=fraction_a,
            scaled_trans_fraction_b=fraction_b,
            R=scaled['R'],
        )

    for path in (mac_figure, comac_figure, report,
                 NPZ_DIR / 'comac_corrected_data.npz'):
        print(f'Wrote {path}')


if __name__ == '__main__':
    main()
