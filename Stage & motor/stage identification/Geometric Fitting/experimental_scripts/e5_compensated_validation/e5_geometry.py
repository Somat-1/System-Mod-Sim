"""Standalone evaluator for the frozen E2 geometric model.

Reads model/geom_model_E2_frozen.json and evaluates G(x, theta) without the
E2 pipeline: the same cubic B-spline bases and the same projected slow-trend
transform that render_e2_findings.prepare_design() built when the model was
fitted, stored as numbers. e5_freeze_model.py checks this evaluator against
the E2 code to machine precision before the file is written.

Convention (E2): deviation is measured minus commanded, so G > 0 means the
stage travels further than commanded; theta = 2 pi (usteps mod 3200) / 3200.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import BSpline

HERE = Path(__file__).resolve().parent
MODEL_FILE = HERE / "model" / "geom_model_E2_frozen.json"
U_PER_MM = 1600
SCREW_REV_U = 3200


def theta_of(usteps: np.ndarray) -> np.ndarray:
    return 2 * np.pi * np.mod(np.asarray(usteps, dtype=np.int64), SCREW_REV_U) / SCREW_REV_U


def bspline_basis(x: np.ndarray, interior: int, domain: tuple[float, float]) -> np.ndarray:
    degree = 3
    lo, hi = domain
    inner = np.linspace(lo, hi, interior + 2)[1:-1] if interior else np.empty(0)
    knots = np.r_[[lo] * (degree + 1), inner, [hi] * (degree + 1)]
    return BSpline.design_matrix(np.clip(x, lo, hi), knots, degree,
                                 extrapolate=True).toarray()


class FrozenModel:
    def __init__(self, path: Path = MODEL_FILE):
        raw = path.read_bytes()
        self.sha256 = hashlib.sha256(raw).hexdigest()
        self.data = json.loads(raw)
        s = self.data["state"]
        self.xbar = float(s["xbar"])
        self.domain = (float(s["domain"][0]), float(s["domain"][1]))
        self.H, self.Ks, self.Kc = int(s["H"]), int(s["Ks"]), s["Kc"]
        self.slow_projection = np.asarray(s["slow_projection"])
        self.slow_transform = np.asarray(s["slow_transform"])
        self.beta = np.asarray(self.data["beta"])
        self.a0 = float(self.beta[0])
        self.beta0_um = float(self.data["reversal"]["beta0_um"])

    def design(self, x: np.ndarray, theta: np.ndarray) -> np.ndarray:
        x = np.asarray(x, dtype=float); theta = np.asarray(theta, dtype=float)
        q = np.column_stack((np.ones(x.size), x - self.xbar))
        columns = [q]
        if self.Ks:
            basis = bspline_basis(x, self.Ks, self.domain)
            columns.append((basis - q @ self.slow_projection) @ self.slow_transform)
        if self.Kc is None:
            for k in range(1, self.H + 1):
                columns.append(np.column_stack((np.cos(k * theta), np.sin(k * theta))))
        else:
            cb = bspline_basis(x, int(self.Kc), self.domain)
            for k in (1, 2):
                columns.append(cb * np.cos(k * theta)[:, None])
                columns.append(cb * np.sin(k * theta)[:, None])
            for k in range(3, self.H + 1):
                columns.append(np.column_stack((np.cos(k * theta), np.sin(k * theta))))
        return np.column_stack(columns)

    def G(self, x_mm: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """Direction-averaged geometric error in um, including the offset a0."""
        return self.design(x_mm, theta) @ self.beta

    def G_rel(self, x_mm: np.ndarray, theta: np.ndarray) -> np.ndarray:
        """G - a0: the part a command correction acts on (datum dropped)."""
        return self.G(x_mm, theta) - self.a0

    def in_domain(self, x_mm: np.ndarray, margin_mm: float = 0.0) -> np.ndarray:
        x = np.asarray(x_mm)
        return (x >= self.domain[0] + margin_mm) & (x <= self.domain[1] - margin_mm)
