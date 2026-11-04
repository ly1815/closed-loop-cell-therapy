"""Digital shadow of cell growth (Egan et al., 2025): model fitting to online data, then prediction.

The digital shadow uses the online data of the perfusion rate D(t) and the glucose and lactate
concentrations [S(t)] and [P(t)], plus the off-line initial cell concentration [X_I]. It simulates the whole controlled process, PI controller included,
and fits it to all three online signals at once, with equal weights and no scaling (section 2.3).

Two modes, as in the paper:

- digital shadow training (section 3.2): all off-line and online data of a run, so that [X_F] is
  known; fitted parameters Y_xs, Y_px, t_lag and nu, with mu_max given by eq. 12.
- digital shadow testing (section 3.3): only [X_I] and the online data available at the time of
  model fitting; t_lag and Y_xs fixed at their mean training values (Table S1: the real-time
  model is structurally identifiable only if mu_max or t_lag is known a priori), and Y_px,
  mu_max and nu fitted.

Least squares by scipy (trust-region reflective, which handles the bounds; the paper used
Levenberg-Marquardt via lmfit). The model is evaluated on a 60 s grid rather than 10 s to keep
fits to a fraction of a second; the data are compared at the same instants. Confidence bands for
X(t) and mu(t) use eqs S67-S72: first-order propagation of the (mu_max, nu) covariance, scaled by
the reduced chi-squared and a Student t factor.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import optimize, stats

from ..process.facer import simulate_model
from ..process.growth import (
    cell_concentration,
    mu_max_from_final,
    sensitivities,
    specific_growth_rate,
    time_to_threshold,
)
from ..process.parameters import MEAN_LAG, MEAN_Y_XS, FacerRun

DT_FIT = 60.0 / 86400.0


@dataclass
class OnlineData:
    t: np.ndarray
    D: np.ndarray
    S: np.ndarray
    P: np.ndarray

    @classmethod
    def from_dict(cls, d: dict) -> "OnlineData":
        return cls(np.asarray(d["t"], float), np.asarray(d["D"], float),
                   np.asarray(d["S"], float), np.asarray(d["P"], float))

    def until(self, t_now: float) -> "OnlineData":
        k = self.t <= t_now + 1e-12
        return OnlineData(self.t[k], self.D[k], self.S[k], self.P[k])


@dataclass
class ShadowFit:
    mode: str
    t_now: float
    params: dict[str, float]
    cov: np.ndarray  # covariance of the fitted parameters, in fit order
    names: tuple[str, ...]
    chi2_red: float
    n_points: int
    run: FacerRun
    X_I: float

    def stderr(self) -> dict[str, float]:
        return {n: float(np.sqrt(max(self.cov[i, i], 0.0))) for i, n in enumerate(self.names)}

    def predict_X(self, t) -> np.ndarray:
        p = self.params
        return cell_concentration(t, self.X_I, p["t_lag"], p["mu_max"], p["nu"])

    def predict_mu(self, t) -> np.ndarray:
        p = self.params
        return specific_growth_rate(t, p["t_lag"], p["mu_max"], p["nu"])

    def threshold_time(self, X_T: float) -> float:
        """Predicted time t_T at which [X] reaches the threshold [X_T] (eq. 13)."""
        p = self.params
        return float(time_to_threshold(self.X_I, X_T, p["t_lag"], p["mu_max"], p["nu"]))

    def band(self, t, alpha: float = 0.05):
        """95% confidence bands (S71, S72) on X and mu. Needs mu_max and nu in the covariance."""
        t = np.asarray(t, float)
        p = self.params
        idx = [self.names.index(n) for n in ("mu_max", "nu") if n in self.names]
        if len(idx) < 2:
            raise ValueError("bands need mu_max and nu as fitted parameters (real-time mode)")
        cov = self.cov[np.ix_(idx, idx)]
        gX, gmu = sensitivities(t, self.X_I, p["t_lag"], p["mu_max"], p["nu"])
        dof = max(self.n_points - len(self.names), 1)
        tau = stats.t.ppf(1 - alpha / 2, dof)
        X = self.predict_X(t)
        mu = self.predict_mu(t)
        # S72 multiplies by the reduced chi-squared; cov here already carries that scaling
        hX = tau * np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", gX, cov, gX), 0.0))
        hmu = tau * np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", gmu, cov, gmu), 0.0))
        return {"X": X, "X_lo": X - hX, "X_hi": X + hX, "mu": mu, "mu_lo": mu - hmu, "mu_hi": mu + hmu}


def _grid_values(traj, t_data: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    return (np.interp(t_data, traj.t, traj.D), np.interp(t_data, traj.t, traj.S),
            np.interp(t_data, traj.t, traj.P))


def _prepare(data: OnlineData, run: FacerRun, dt: float) -> OnlineData:
    """Keep samples from seeding onwards, drop gaps, thin to the model grid."""
    ok = (data.t >= run.t_I) & np.isfinite(data.S) & np.isfinite(data.P) & np.isfinite(data.D)
    d = OnlineData(data.t[ok], data.D[ok], data.S[ok], data.P[ok])
    step = max(1, int(round(dt / np.median(np.diff(d.t))))) if len(d.t) > 1 else 1
    return OnlineData(d.t[::step], d.D[::step], d.S[::step], d.P[::step])


def _fit(run, X_I, data, t_end, names, x0, lower, upper, build, mode, dt):
    data = _prepare(data, run, dt)
    obs = np.concatenate([data.D, data.S, data.P])
    shadow_run = run.__class__(**{**run.__dict__, "X_I": X_I})

    def residuals(theta):
        kw = build(theta)
        traj = simulate_model(shadow_run, t_end=t_end, dt=dt, **kw)
        return np.concatenate(_grid_values(traj, data.t)) - obs

    res = optimize.least_squares(residuals, x0, bounds=(lower, upper), method="trf",
                                 x_scale="jac", max_nfev=400)
    n, k = len(obs), len(names)
    chi2_red = float(np.sum(res.fun**2) / max(n - k, 1))
    J = res.jac
    try:
        cov = np.linalg.inv(J.T @ J) * chi2_red
    except np.linalg.LinAlgError:
        cov = np.full((k, k), np.nan)
    params = build(res.x)
    params.setdefault("t_lag", run.t_lag)
    return ShadowFit(mode, float(data.t[-1]), {**params, "Y_xs": params.get("Y_xs", run.Y_xs)},
                     cov, tuple(names), chi2_red, n, shadow_run, X_I)


def fit_training(run: FacerRun, data: OnlineData, X_I: float, X_F: float, dt: float = DT_FIT) -> ShadowFit:
    """Fit a complete run (both counts known): Y_xs, Y_px, t_lag - t_I and nu."""

    def build(theta):
        Y_xs, Y_px, lag, nu = theta
        t_lag = run.t_I + lag
        return {"Y_xs": Y_xs, "Y_px": Y_px, "t_lag": t_lag, "nu": nu,
                "mu_max": mu_max_from_final(X_I, X_F, run.t_F, t_lag, nu)}

    names = ("Y_xs", "Y_px", "lag", "nu")
    x0 = np.array([0.25, 7.5, 0.1, -0.5])
    lower = np.array([0.01, 0.1, 0.0, -3.0])
    upper = np.array([5.0, 50.0, run.t_F - run.t_I - 0.5, -1e-6])
    return _fit(run, X_I, data, run.t_F, names, x0, lower, upper, build, "training", dt)


def fit_realtime(run: FacerRun, data: OnlineData, X_I: float, t_now: float,
                 lag: float = MEAN_LAG, Y_xs: float = MEAN_Y_XS, dt: float = DT_FIT) -> ShadowFit:
    """Fit the data available at t_now: Y_px, mu_max and nu, with t_lag and Y_xs fixed."""
    t_lag = run.t_I + lag

    def build(theta):
        Y_px, mu_max, nu = theta
        return {"Y_xs": Y_xs, "Y_px": Y_px, "t_lag": t_lag, "mu_max": mu_max, "nu": nu}

    names = ("Y_px", "mu_max", "nu")
    x0 = np.array([7.5, 1.5, -0.5])
    lower = np.array([0.1, 0.01, -3.0])
    upper = np.array([50.0, 10.0, -1e-6])
    return _fit(run, X_I, data.until(t_now), t_now, names, x0, lower, upper, build, "realtime", dt)
