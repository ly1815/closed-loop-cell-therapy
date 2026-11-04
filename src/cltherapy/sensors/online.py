"""Online sensors of the FACER cartridge, as described by Egan et al. (2025), section 2.2.

- Glucose and lactate: single-use enzymatic sensors, a reading every 10 s.
- Bubbles in the chamber pass the sensors and cause artificial drops in both readings, followed
  by a recovery once the bubble moves away (seconds to minutes).
- The raw glucose and lactate signals go through a 2nd-order low-pass Butterworth filter
  (cut-off 0.001 Hz) before the PID sees them. Not every bubble artefact is filtered out, which
  is why the published data show downward glucose spikes and upward perfusion spikes.
- Dissolved oxygen is measured online too (normalised to saturation here).

Noise levels and bubble statistics are placeholders chosen to look like the published traces;
the filter design is the published one. The filter is causal and runs sample by sample, as on
the instrument. The paper quotes a 1 Hz sampling frequency for the filter; this simulator
filters the 10 s readings, keeping the same 0.001 Hz cut-off (a time constant of about 160 s).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy import signal


class CausalButterworth:
    """2nd-order low-pass Butterworth, applied one sample at a time."""

    def __init__(self, cutoff_hz: float = 0.001, fs_hz: float = 0.1, order: int = 2) -> None:
        self.sos = signal.butter(order, cutoff_hz, btype="low", fs=fs_hz, output="sos")
        self._zi_unit = signal.sosfilt_zi(self.sos)
        self.zi: list[list[float]] | None = None

    def __call__(self, x: float) -> float:
        if self.zi is None:  # start at steady state on the first sample
            self.zi = (self._zi_unit * x).tolist()
        y = x
        for (b0, b1, b2, _a0, a1, a2), z in zip(self.sos, self.zi):  # transposed direct form II
            out = b0 * y + z[0]
            z[0] = b1 * y - a1 * out + z[1]
            z[1] = b2 * y - a2 * out
            y = out
        return float(y)


@dataclass
class BubbleModel:
    rate_per_day: float = 3.0  # bubbles passing the sensors
    depth: tuple[float, float] = (0.3, 0.8)  # reading drops to this fraction of the truth
    duration_s: tuple[float, float] = (20.0, 300.0)


@dataclass
class OnlineSensors:
    seed: int = 0
    glucose_sd: float = 0.08  # mmol/L
    lactate_sd: float = 0.08  # mmol/L
    DO_sd: float = 0.01  # normalised
    bubbles: BubbleModel = field(default_factory=BubbleModel)
    dt_s: float = 10.0

    def __post_init__(self) -> None:
        self.rng = np.random.default_rng(self.seed)
        self.filt_S = CausalButterworth(fs_hz=1.0 / self.dt_s)
        self.filt_P = CausalButterworth(fs_hz=1.0 / self.dt_s)
        self._bubble_left_s = 0.0
        self._bubble_depth = 1.0
        self.offline = False  # True while data acquisition is down (a mechanical stop)

    def _bubble_factor(self) -> float:
        if self._bubble_left_s <= 0.0:
            p = self.bubbles.rate_per_day * self.dt_s / 86400.0
            if self.rng.random() < p:
                self._bubble_left_s = self.rng.uniform(*self.bubbles.duration_s)
                self._bubble_depth = self.rng.uniform(*self.bubbles.depth)
        if self._bubble_left_s > 0.0:
            self._bubble_left_s -= self.dt_s
            return self._bubble_depth
        return 1.0

    def read(self, S: float, P: float, C: float) -> dict:
        """One 10 s reading: raw and filtered glucose and lactate, and DO. None while offline."""
        if self.offline:
            return {"S_raw": None, "P_raw": None, "S_filt": None, "P_filt": None, "C_meas": None}
        b = self._bubble_factor()
        S_raw = max(0.0, S * b + self.rng.normal(0.0, self.glucose_sd))
        P_raw = max(0.0, P * b + self.rng.normal(0.0, self.lactate_sd))
        C_meas = float(np.clip(C + self.rng.normal(0.0, self.DO_sd), 0.0, 1.2))
        return {"S_raw": S_raw, "P_raw": P_raw, "S_filt": self.filt_S(S_raw),
                "P_filt": self.filt_P(P_raw), "C_meas": C_meas}


@dataclass
class OfflineCounts:
    """Cell counts need the cells in suspension, so they are taken only at seeding and after
    perfusion ends; resuspending mid-run would disturb growth (Egan et al. 2025, section 2.2)."""

    seed: int = 1
    cv: float = 0.05

    def __post_init__(self) -> None:
        self.rng = np.random.default_rng(self.seed)

    def count(self, X: float) -> float:
        return float(X * (1.0 + self.rng.normal(0.0, self.cv)))
