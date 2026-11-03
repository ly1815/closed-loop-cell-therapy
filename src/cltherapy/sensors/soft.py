"""Soft measurements derived from online signals: oxygen uptake and acid production.

These are what the EnKF will use between at-line samples (day 6).

OUR from the dissolved-oxygen balance over a window:
    OUR = kLa (C* - C) - dC/dt
with kLa from the rock rate and angle via a correlation that may be biased
for this particular bag (kla_bias), and C* from the O2 fraction in the gas.

Acid production from the CO2 the pH loop needs: invert Henderson-Hasselbalch
for bicarbonate, then
    acid = -(dB/dt - D (B_in - B)) / acid_to_bicarb
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from ..process.model import kla
from ..process.parameters import DEFAULT_PROCESS, ProcessParams


@dataclass
class SoftSensors:
    process: ProcessParams = DEFAULT_PROCESS
    our_window_h: float = 0.5  # DO is precise enough for a short window
    acid_window_h: float = 6.0  # pH has 0.01 resolution, so the bicarbonate slope needs a long window
    kla_bias: float = 1.0  # multiplicative error of the kLa correlation for this bag

    def _bicarbonate(self, pH: np.ndarray, y_CO2: np.ndarray) -> np.ndarray:
        p = self.process
        return 10.0 ** (pH - p.pKa_bicarb) * p.CO2_solubility * np.maximum(y_CO2, 1e-4) * p.pCO2_per_frac

    @staticmethod
    def _slope(t: np.ndarray, y: np.ndarray, window_h: float) -> np.ndarray:
        """Least-squares slope over a trailing window (causal, as it would run live)."""
        out = np.full_like(y, np.nan, dtype=float)
        j = 0
        for i in range(len(t)):
            while t[i] - t[j] > window_h:
                j += 1
            if t[i] - t[j] >= 0.8 * window_h and i - j >= 2:
                tt, yy = t[j : i + 1], y[j : i + 1]
                out[i] = np.polyfit(tt - tt[-1], yy, 1)[0]
        return out

    @staticmethod
    def _mean(t: np.ndarray, y: np.ndarray, window_h: float) -> np.ndarray:
        out = np.full_like(y, np.nan, dtype=float)
        j = 0
        for i in range(len(t)):
            while t[i] - t[j] > window_h:
                j += 1
            out[i] = np.mean(y[j : i + 1])
        return out

    def estimate(self, t_h, DO_pct, pH, y_O2, y_CO2, rock_rpm, rock_angle_deg, perfusion_vvd) -> dict:
        """Arrays in, arrays out (mM/h). Inputs are what the platform logs."""
        p = self.process
        t = np.asarray(t_h, float)
        C = np.asarray(DO_pct, float) / 100.0 * p.O2_sat_air
        C_star = p.O2_sat_air * np.asarray(y_O2, float) / 0.21
        k = np.array([kla(p, r, a) for r, a in zip(rock_rpm, rock_angle_deg)]) * self.kla_bias
        w = self.our_window_h
        our = self._mean(t, k * (C_star - C), w) - self._slope(t, C, w)

        B = self._bicarbonate(np.asarray(pH, float), np.asarray(y_CO2, float))
        D = np.asarray(perfusion_vvd, float) / 24.0
        w = self.acid_window_h
        acid = -(self._slope(t, B, w) - self._mean(t, D * (p.hco3_in - B), w)) / p.acid_to_bicarb
        return {"OUR_hat": our, "acid_hat": acid}
