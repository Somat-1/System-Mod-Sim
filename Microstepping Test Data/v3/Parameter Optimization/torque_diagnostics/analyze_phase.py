#!/usr/bin/env python3
"""Quantify the actual phase relationship between motor drive torque and
detent torque in the cruise-window simulation, rather than eyeballing
Figure B/C. Fits A*sin(2*pi*f*t + phi) + C to each of motor_force and the
RAW (unsigned, constitutive-law) detent torque at the known detent
frequency, and reports the phase difference -- exactly 180 deg would mean
perfectly out of phase; anything else is a real, physically meaningful
deviation from that idealization.
"""
from __future__ import annotations

import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from plot_cruise_zoom import RUN_INDEX, RATES, cruise_zoom_data  # noqa: E402


def fit_sine(t, y, f_hz):
    omega = 2.0 * np.pi * f_hz
    amp0 = 0.5 * (np.max(y) - np.min(y))
    off0 = np.mean(y)

    def model(t, amp, phi, off):
        return amp * np.sin(omega * t + phi) + off

    popt, _ = curve_fit(model, t, y, p0=[amp0, 0.0, off0], maxfev=20000)
    amp, phi, off = popt
    if amp < 0:
        amp, phi = -amp, phi + np.pi
    return amp, np.degrees(phi) % 360.0, off


def main():
    started = time.perf_counter()
    results = {}
    with ProcessPoolExecutor(max_workers=2) as pool:
        futures = {pool.submit(cruise_zoom_data, RUN_INDEX, rate): rate for rate in RATES}
        for fut in as_completed(futures):
            rate = futures[fut]
            results[rate] = fut.result()
            print(f'done D_{rate} ({time.perf_counter() - started:.1f}s)', flush=True)

    print(f'\n{"rate":>6s}  {"motor amp (N)":>14s}  {"motor phase":>12s}  '
          f'{"detent amp (N)":>15s}  {"detent phase":>13s}  {"phase diff":>11s}  {"dev from 180":>13s}')
    for rate in RATES:
        res = results[rate]
        t = res['t']
        f_hz = res['detent_hz']
        motor_amp, motor_phase, _ = fit_sine(t, res['motor_force'], f_hz)
        # detent_force here is the SIGNED (-T_detent/lead) contribution;
        # fit the raw constitutive law's own sinusoid by flipping sign back.
        detent_amp, detent_phase, _ = fit_sine(t, -res['detent_force'], f_hz)
        diff = (motor_phase - detent_phase) % 360.0
        dev_from_180 = min(abs(diff - 180.0), abs(diff - 180.0 - 360.0), abs(diff - 180.0 + 360.0))
        print(f'D_{rate:>4s}  {motor_amp:>14.4f}  {motor_phase:>11.2f}°  '
              f'{detent_amp:>15.4f}  {detent_phase:>12.2f}°  {diff:>10.2f}°  {dev_from_180:>12.2f}°')


if __name__ == '__main__':
    main()
