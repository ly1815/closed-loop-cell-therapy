"""Online sensors on the rocking bag: what the controller and the platform actually see.

Each channel has white noise, a slow drift that grows with time since
calibration, its own update period (a reading is held between updates) and a
display resolution. Values are placeholders typical of single-use optical
patches, a load cell and gas mass-flow controllers; see docs/sensors.md.

Faults (used from day 9): add a bias or a drift step, freeze a channel
(stale), or drop it (None, as a lost connection would).
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


@dataclass
class ChannelSpec:
    key: str  # name in the truth dict from model.observables
    noise_sd: float
    drift_per_day: float  # signed drift rate, applied as offset = rate * days
    period_h: float  # how often a new value arrives
    resolution: float  # display rounding
    unit: str = ""


DEFAULT_CHANNELS: tuple[ChannelSpec, ...] = (
    ChannelSpec("pH", noise_sd=0.01, drift_per_day=-0.005, period_h=10 / 3600, resolution=0.01, unit=""),
    ChannelSpec("DO_pct", noise_sd=1.0, drift_per_day=0.5, period_h=10 / 3600, resolution=0.1, unit="% air sat."),
    ChannelSpec("weight_kg", noise_sd=0.002, drift_per_day=0.0, period_h=1 / 3600, resolution=0.001, unit="kg"),
)


@dataclass
class Fault:
    kind: str  # "bias", "drift", "stale", "dropout"
    value: float = 0.0  # bias size, or extra drift per day
    start_h: float = 0.0


@dataclass
class OnlineSensors:
    channels: tuple[ChannelSpec, ...] = DEFAULT_CHANNELS
    seed: int = 0
    faults: dict[str, Fault] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.rng = np.random.default_rng(self.seed)
        self._held: dict[str, tuple[float, float | None]] = {}  # key -> (time of last update, value)

    def inject(self, key: str, fault: Fault) -> None:
        self.faults[key] = fault

    def clear(self, key: str) -> None:
        self.faults.pop(key, None)

    def _raw(self, spec: ChannelSpec, t_h: float, truth: float) -> float:
        value = truth + spec.drift_per_day * t_h / 24.0 + self.rng.normal(0.0, spec.noise_sd)
        fault = self.faults.get(spec.key)
        if fault is not None and t_h >= fault.start_h:
            if fault.kind == "bias":
                value += fault.value
            elif fault.kind == "drift":
                value += fault.value * (t_h - fault.start_h) / 24.0
        return round(value / spec.resolution) * spec.resolution

    def read(self, t_h: float, truth: dict) -> dict:
        """Readings at time t_h, keyed like the truth dict, with a timestamp."""
        truth = dict(truth)
        truth.setdefault("weight_kg", truth.get("volume_L", 1.0) * 1.0)
        out: dict = {"t_h": t_h}
        for spec in self.channels:
            fault = self.faults.get(spec.key)
            active = fault is not None and t_h >= fault.start_h
            if active and fault.kind == "dropout":
                out[spec.key] = None
                continue
            last_t, last_v = self._held.get(spec.key, (-np.inf, None))
            stale = active and fault.kind == "stale" and last_v is not None
            if not stale and t_h - last_t >= spec.period_h - 1e-12:
                last_t, last_v = t_h, self._raw(spec, t_h, truth[spec.key])
                self._held[spec.key] = (last_t, last_v)
            out[spec.key] = last_v
            out[f"{spec.key}_age_h"] = t_h - last_t
        return out
