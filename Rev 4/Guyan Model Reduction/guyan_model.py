#!/usr/bin/env python3
'''Classical two-master Guyan reduction supplied by the user.

Masters: [theta_m, x_n].
Slaves:  [theta_c, theta_s, theta_sb, x_s].
No mass tuning or IRS correction is applied.
'''

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
from scipy.linalg import eigh


ROOT = Path(__file__).resolve().parent
REV4_ROOT = ROOT.parent
PARAMETER_FILE = REV4_ROOT / 'model_parameters.json'
MASTER_INDICES = np.array([0, 5])
SLAVE_INDICES = np.array([1, 2, 3, 4])
DOF_LABELS = ('theta_m', 'theta_c', 'theta_s', 'theta_sb', 'x_s', 'x_n')


def load_parameters() -> dict[str, float]:
    return json.loads(PARAMETER_FILE.read_text(encoding='utf-8'))['parameters']


def build_full_matrices(p: dict[str, float]):
    '''Build the frictionless 6-DOF M/C/K and physical b=[k_EM,0,...].'''
    ell = p['L'] / (2.0 * np.pi)
    k_em = p['N_r'] * p['T_hold']
    k_d = 4.0 * p['N_r'] * p['T_d']
    mass = np.diag([
        p['I_m'], p['I_c'], p['I_s'], p['I_sb'],
        p['M_screw'], p['M_s'],
    ])
    k_c, k_s1, k_s2 = p['k_c'], p['k_s1'], p['k_s2']
    k_nut, k_brg = p['k_nut'], p['k_brg']
    stiffness = np.array([
        [k_c + k_em + k_d, -k_c, 0.0, 0.0, 0.0, 0.0],
        [-k_c, k_c + k_s1, -k_s1, 0.0, 0.0, 0.0],
        [0.0, -k_s1, k_s1 + k_s2 + ell**2 * k_nut,
         -k_s2, ell * k_nut, -ell * k_nut],
        [0.0, 0.0, -k_s2, k_s2, 0.0, 0.0],
        [0.0, 0.0, ell * k_nut, 0.0, k_brg + k_nut, -k_nut],
        [0.0, 0.0, -ell * k_nut, 0.0, -k_nut, k_nut],
    ])
    c_c, c_s1, c_s2 = p['c_c'], p['c_s1'], p['c_s2']
    c_nut, c_brg, c_em = p['c_nut'], p['c_brg'], p['c_EM']
    damping = np.array([
        [c_c + c_em, -c_c, 0.0, 0.0, 0.0, 0.0],
        [-c_c, c_c + c_s1, -c_s1, 0.0, 0.0, 0.0],
        [0.0, -c_s1, c_s1 + c_s2 + ell**2 * c_nut,
         -c_s2, ell * c_nut, -ell * c_nut],
        [0.0, 0.0, -c_s2, c_s2, 0.0, 0.0],
        [0.0, 0.0, ell * c_nut, 0.0, c_brg + c_nut, -c_nut],
        [0.0, 0.0, -ell * c_nut, 0.0, -c_nut, c_nut],
    ])
    command = np.array([k_em, 0.0, 0.0, 0.0, 0.0, 0.0])
    return mass, damping, stiffness, command


def closed_form_ratios(p: dict[str, float]) -> dict[str, float]:
    ell = p['L'] / (2.0 * np.pi)
    beta = p['k_nut'] / (p['k_nut'] + p['k_brg'])
    kappa = 1.0 / (1.0 / p['k_nut'] + 1.0 / p['k_brg'])
    k_ch = 1.0 / (1.0 / p['k_c'] + 1.0 / p['k_s1'])
    nu = ell**2 * kappa / (k_ch + ell**2 * kappa)
    mu = p['k_s1'] / (p['k_c'] + p['k_s1'])
    k_ax = kappa * (1.0 - nu)
    gamma = k_ax / p['k_nut']
    return {
        'ell': ell,
        'beta': beta,
        'kappa': kappa,
        'k_ch': k_ch,
        'nu': nu,
        'mu': mu,
        'k_ax': k_ax,
        'gamma': gamma,
    }


def closed_form_transformation(p: dict[str, float]) -> np.ndarray:
    r = closed_form_ratios(p)
    ell, beta, nu, mu = r['ell'], r['beta'], r['nu'], r['mu']
    return np.array([
        [1.0, 0.0],
        [1.0 - mu * nu, mu * nu / ell],
        [1.0 - nu, nu / ell],
        [1.0 - nu, nu / ell],
        [-ell * beta * (1.0 - nu), beta * (1.0 - nu)],
        [0.0, 1.0],
    ])


def numerical_transformation(stiffness: np.ndarray) -> np.ndarray:
    k_ss = stiffness[np.ix_(SLAVE_INDICES, SLAVE_INDICES)]
    k_sm = stiffness[np.ix_(SLAVE_INDICES, MASTER_INDICES)]
    slave_map = -np.linalg.solve(k_ss, k_sm)
    transformation = np.zeros((6, 2))
    transformation[MASTER_INDICES] = np.eye(2)
    transformation[SLAVE_INDICES] = slave_map
    return transformation


def closed_form_stiffness(p: dict[str, float]) -> np.ndarray:
    r = closed_form_ratios(p)
    ell, k_ax = r['ell'], r['k_ax']
    k_em = p['N_r'] * p['T_hold']
    k_d = 4.0 * p['N_r'] * p['T_d']
    return np.array([
        [k_em + k_d + ell**2 * k_ax, -ell * k_ax],
        [-ell * k_ax, k_ax],
    ])


def closed_form_mass(p: dict[str, float]) -> np.ndarray:
    r = closed_form_ratios(p)
    ell, beta, nu, mu = r['ell'], r['beta'], r['nu'], r['mu']
    i_sum = p['I_s'] + p['I_sb']
    m11 = (
        p['I_m']
        + p['I_c'] * (1.0 - mu * nu) ** 2
        + i_sum * (1.0 - nu) ** 2
        + p['M_screw'] * ell**2 * beta**2 * (1.0 - nu) ** 2
    )
    m12 = (
        p['I_c'] * (1.0 - mu * nu) * mu * nu / ell
        + i_sum * (1.0 - nu) * nu / ell
        - p['M_screw'] * ell * beta**2 * (1.0 - nu) ** 2
    )
    m22 = (
        p['I_c'] * (mu * nu / ell) ** 2
        + i_sum * (nu / ell) ** 2
        + p['M_screw'] * beta**2 * (1.0 - nu) ** 2
        + p['M_s']
    )
    return np.array([[m11, m12], [m12, m22]])


def reduce_model(p: dict[str, float]) -> dict[str, np.ndarray]:
    mass, damping, stiffness, command = build_full_matrices(p)
    transformation = closed_form_transformation(p)
    reduced = {
        'T': transformation,
        'M': transformation.T @ mass @ transformation,
        'C': transformation.T @ damping @ transformation,
        'K': transformation.T @ stiffness @ transformation,
        'b': transformation.T @ command,
        'M_full': mass,
        'C_full': damping,
        'K_full': stiffness,
        'b_full': command,
    }
    ell = p['L'] / (2.0 * np.pi)
    ports = {
        'way': np.array([0.0, 0.0, 0.0, 0.0, 0.0, 1.0]),
        'nut': np.array([0.0, 0.0, ell, 0.0, 1.0, -1.0]),
        'sb': np.array([0.0, 0.0, 0.0, 1.0, 0.0, 0.0]),
    }
    reduced['J_way'] = ports['way'] @ transformation
    reduced['J_nut'] = ports['nut'] @ transformation
    reduced['J_sb'] = ports['sb'] @ transformation
    return reduced


def fixed_interface_frequencies_hz(
    mass: np.ndarray, stiffness: np.ndarray
) -> np.ndarray:
    m_ss = mass[np.ix_(SLAVE_INDICES, SLAVE_INDICES)]
    k_ss = stiffness[np.ix_(SLAVE_INDICES, SLAVE_INDICES)]
    eigenvalues = eigh(k_ss, m_ss, eigvals_only=True)
    return np.sqrt(np.maximum(eigenvalues, 0.0)) / (2.0 * np.pi)


def frequency_response(
    mass: np.ndarray,
    damping: np.ndarray,
    stiffness: np.ndarray,
    command: np.ndarray,
    output: np.ndarray,
    frequencies_hz: np.ndarray,
) -> np.ndarray:
    response = np.empty(len(frequencies_hz), dtype=np.complex128)
    for index, frequency in enumerate(frequencies_hz):
        omega = 2.0 * np.pi * frequency
        dynamic_stiffness = stiffness + 1j * omega * damping - omega**2 * mass
        response[index] = output @ np.linalg.solve(dynamic_stiffness, command)
    return response


def modal_frequencies_hz(mass: np.ndarray, stiffness: np.ndarray) -> np.ndarray:
    eigenvalues = eigh(stiffness, mass, eigvals_only=True)
    return np.sqrt(np.maximum(eigenvalues, 0.0)) / (2.0 * np.pi)
