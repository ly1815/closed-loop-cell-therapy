"""Batch-perfusion process in the large chamber of the FACER, after Egan et al. (2025).

Sedimentation period (t_I <= t <= t_sed) followed by the perfusion period (eqs 2-4, 8; S1-S24):

    d[S]/dt = H(t - t_sed) D(t) ([S_I] - [S(t)]) - mu(t)[X(t)] / Y_xs           (eq. 3)
    d[P]/dt = Y_px mu(t)[X(t)] - H(t - t_sed) D(t) [P(t)]                        (eq. 4)
    dD/dt   = H(t - t_sed) ( -K_p d[S]/dt + K_i ([S_sp] - [S(t)]) ),  D >= D_min (eq. 2)
    dC'/dt  = (H(t - t_sed) D(t) + kLa)(1 - C'(t)) - mu(t)[X(t)] / Y'_xc          (eq. 8)

where C' = [C]/[C*] is the normalised DO. [X(t)] and mu(t) are given analytically (growth.py).
The ODEs are approximated by forward finite difference equations and solved via the forward
Euler method with a step size of 10 s, with the 'max' function enforcing the physical
constraints at each time step (S22-S24). This is the model fitted by the digital shadow.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .growth import cell_concentration, specific_growth_rate
from .parameters import FacerRun

DT_10S = 10.0 / 86400.0  # days


@dataclass
class Trajectory:
    t: np.ndarray
    X: np.ndarray
    mu: np.ndarray
    S: np.ndarray
    P: np.ndarray
    D: np.ndarray
    C: np.ndarray  # normalised DO, C'(t) = [C(t)]/[C*]

    def at(self, t: float) -> dict:
        i = int(np.clip(np.searchsorted(self.t, t), 0, len(self.t) - 1))
        return {k: float(getattr(self, k)[i]) for k in ("t", "X", "mu", "S", "P", "D", "C")}


def simulate_model(
    run: FacerRun,
    *,
    Y_xs: float | None = None,
    Y_px: float | None = None,
    t_lag: float | None = None,
    mu_max: float | None = None,
    nu: float | None = None,
    kLa: float | None = None,
    Y_xc: float | None = None,
    t_end: float | None = None,
    dt: float = DT_10S,
) -> Trajectory:
    """Noise-free model trajectory from t_I to t_end (default t_F). Keyword arguments override
    the run's parameters, which is how the shadow evaluates candidate parameter sets."""
    Y_xs = run.Y_xs if Y_xs is None else Y_xs
    Y_px = run.Y_px if Y_px is None else Y_px
    t_lag = run.t_lag if t_lag is None else t_lag
    nu = run.nu if nu is None else nu
    mu_max = run.mu_max if mu_max is None else mu_max
    kLa = run.kLa if kLa is None else kLa
    Y_xc = run.Y_xc if Y_xc is None else Y_xc
    t_end = run.t_F if t_end is None else t_end
    g = run.gains

    t = np.arange(run.t_I, t_end + 0.5 * dt, dt)
    X = cell_concentration(t, run.X_I, t_lag, mu_max, nu)
    mu = specific_growth_rate(t, t_lag, mu_max, nu)
    growth = mu * X  # 1e6 cells/mL/day

    n = len(t)
    S = np.empty(n)
    P = np.empty(n)
    D = np.empty(n)
    C = np.empty(n)
    S[0], P[0], D[0], C[0] = run.S_I, 0.0, 0.0, 1.0
    s, p, d, c = run.S_I, 0.0, 0.0, 1.0
    started = False
    for k in range(n - 1):
        perf = t[k] >= run.t_sed
        if perf and not started:
            d, started = g.D_min, True
        dk = d if perf else 0.0
        s_next = max(0.0, s + dt * (dk * (run.S_I - s) - growth[k] / Y_xs))
        if perf:
            d = max(g.D_min, d - g.kp * (s_next - s) + dt * g.ki * (run.S_sp - s))
        p = max(0.0, p + dt * (Y_px * growth[k] - dk * p))
        c = min(1.0, max(0.0, c + dt * ((dk + kLa) * (1.0 - c) - growth[k] / Y_xc)))
        s = s_next
        S[k + 1], P[k + 1], D[k + 1], C[k + 1] = s, p, (d if t[k + 1] >= run.t_sed else 0.0), c
    return Trajectory(t, X, mu, S, P, D, C)
