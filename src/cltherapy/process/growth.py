"""Cell growth model: the special case of the Richards growth model with a lag period.

Egan et al., Biotechnology Progress 2025 (doi:10.1002/btpr.70045), eqs 10-13 and S60-S72.
In the early stages of cell growth, when [X(t)] is much lower than the carrying capacity [K],
the specific growth rate decreases over time as the cells compete for space:

    [X(t)] = [X_I] / (1 - H(t - t_lag) mu_max nu (t - t_lag))^(1/nu)      (eq. 10, S60)
    mu(t)  = H(t - t_lag) mu_max / (1 - H(t - t_lag) mu_max nu (t - t_lag)) (eq. 11, S61)

for t >= t_I and nu < 0, where mu_max is the maximum specific growth rate at t_lag and nu
determines the decrease of the specific growth rate over time. nu = -1 gives linear cell growth,
nu = -1/2 quadratic growth, and nu -> 0 the limiting case of exponential growth. Time t in days,
[X] in 10^6 cells/mL. Function arguments X_I, X_F, X_T denote [X_I], [X_F], [X_T].
"""

from __future__ import annotations

import numpy as np

NU_EXP = 1e-6  # |nu| below this is treated as exponential growth (the nu -> 0 limit)


def _elapsed(t, t_lag):
    return np.maximum(np.asarray(t, float) - t_lag, 0.0)


def cell_concentration(t, X_I: float, t_lag: float, mu_max: float, nu: float):
    """Cell concentration [X(t)], eq. 10 (S60); equal to [X_I] during the lag period."""
    tau = _elapsed(t, t_lag)
    if abs(nu) < NU_EXP:
        return X_I * np.exp(mu_max * tau)
    return X_I / (1.0 - mu_max * nu * tau) ** (1.0 / nu)


def specific_growth_rate(t, t_lag: float, mu_max: float, nu: float):
    """Specific growth rate mu(t) (day^-1), eq. 11 (S61); zero during the lag period."""
    t = np.asarray(t, float)
    tau = _elapsed(t, t_lag)
    on = (t >= t_lag).astype(float)
    if abs(nu) < NU_EXP:
        return on * mu_max * np.ones_like(tau)
    return on * mu_max / (1.0 - mu_max * nu * tau)


def mu_max_from_final(X_I: float, X_F: float, t_F: float, t_lag: float, nu: float) -> float:
    """Maximum specific growth rate from [X_I] and [X_F], eq. 12 (S66)."""
    if abs(nu) < NU_EXP:
        return np.log(X_F / X_I) / (t_F - t_lag)
    return (1.0 - (X_I / X_F) ** nu) / (nu * (t_F - t_lag))


def time_to_threshold(X_I: float, X_T: float, t_lag: float, mu_max: float, nu: float) -> float:
    """Time t_T at which [X] reaches the threshold cell concentration [X_T], eq. 13."""
    if X_T <= X_I:
        return t_lag
    if abs(nu) < NU_EXP:
        return t_lag + np.log(X_T / X_I) / mu_max
    return t_lag + (1.0 - (X_I / X_T) ** nu) / (mu_max * nu)


def lactate_steady_state(Y_ps: float, S_I: float, S_sp: float) -> float:
    """Lactate concentration steady state, eq. 18: [P_ss] = Y_ps ([S_I] - [S_sp])."""
    return Y_ps * (S_I - S_sp)


def sensitivities(t, X_I: float, t_lag: float, mu_max: float, nu: float):
    """Partial derivatives of [X] and mu with respect to mu_max and nu, eqs S67-S70.

    Returns two arrays of shape (len(t), 2) used for the confidence bands (S71, S72).
    """
    t = np.asarray(t, float)
    tau = _elapsed(t, t_lag)
    X = cell_concentration(t, X_I, t_lag, mu_max, nu)
    mu = specific_growth_rate(t, t_lag, mu_max, nu)
    nu_eff = nu if abs(nu) >= NU_EXP else -NU_EXP
    a = 1.0 - mu_max * nu_eff * tau
    dX_dmu = X * tau / a
    dX_dnu = X / nu_eff**2 * (np.log(a) + mu_max * nu_eff * tau / a)
    dmu_dmu = np.where(mu_max > 0, mu**2 / mu_max, 0.0)
    dmu_dnu = mu**2 * tau
    return np.column_stack([dX_dmu, dX_dnu]), np.column_stack([dmu_dmu, dmu_dnu])
