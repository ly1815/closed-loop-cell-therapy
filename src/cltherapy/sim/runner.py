"""The virtual FACER: one run of the large chamber, stepped every 10 s.

The process follows the Egan et al. (2025) balances (process/facer.py), but here the PID acts on
filtered, noisy, bubble-affected glucose readings, as on the real instrument. Everything is
recorded so a run can be replayed, plotted, streamed to virtual devices or fed to the shadow.

Faults:
    pump_stop(start, duration)  the pump stops and data acquisition is lost, as in run 3
    growth_noise                 day-to-day variation of the growth rate (off by default)
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from ..control.perfusion import GlucosePerfusionPID
from ..process.facer import DT_10S
from ..process.growth import specific_growth_rate
from ..process.parameters import FacerRun
from ..sensors.online import OfflineCounts, OnlineSensors

FIELDS = ("t", "X", "mu", "S", "P", "D", "C", "S_raw", "P_raw", "S_filt", "P_filt", "C_meas",
          "D_cmd", "pump_on", "phase")


@dataclass
class PumpStop:
    start: float  # day
    duration: float  # days


@dataclass
class VirtualFacer:
    run: FacerRun
    sensors: OnlineSensors = field(default_factory=OnlineSensors)
    counts: OfflineCounts = field(default_factory=OfflineCounts)
    dt: float = DT_10S
    growth_noise: float = 0.0  # sd of a slow random walk on log(mu), per sqrt(day)
    pump_stops: list[PumpStop] = field(default_factory=list)
    seed: int = 0

    def __post_init__(self) -> None:
        r = self.run
        self.rng = np.random.default_rng(self.seed)
        self.pid = GlucosePerfusionPID(S_sp=r.S_sp, gains=r.gains)
        self.t = r.t_I
        self.X, self.S, self.P, self.D, self.C = r.X_I, r.S_I, 0.0, 0.0, 1.0
        self._log_mu_dev = 0.0
        self.history: dict[str, list] = {k: [] for k in FIELDS}
        self.X_I_count = self.counts.count(r.X_I)  # suspension count at seeding
        self.X_F_count: float | None = None
        self._reading = self.sensors.read(self.S, self.P, self.C)
        self._record(D_cmd=0.0, pump_on=False)

    # ------------------------------------------------------------------ helpers
    def phase(self) -> str:
        if self.t < self.run.t_sed:
            return "sedimentation"
        if self.t < self.run.t_F:
            return "perfusion"
        return "done"

    def _pump_stopped(self) -> bool:
        return any(p.start <= self.t < p.start + p.duration for p in self.pump_stops)

    def _mu(self) -> float:
        r = self.run
        tau = self.t - r.t_lag
        if tau < 0:
            return 0.0
        mu = float(specific_growth_rate(self.t, r.t_lag, r.mu_max, r.nu))
        return mu * np.exp(self._log_mu_dev)

    def _record(self, D_cmd: float, pump_on: bool) -> None:
        h = self.history
        for k, v in (("t", self.t), ("X", self.X), ("mu", self._mu()), ("S", self.S), ("P", self.P),
                     ("D", self.D), ("C", self.C), ("D_cmd", D_cmd), ("pump_on", pump_on),
                     ("phase", self.phase())):
            h[k].append(v)
        for k in ("S_raw", "P_raw", "S_filt", "P_filt", "C_meas"):
            v = self._reading[k]
            h[k].append(np.nan if v is None else v)

    # ------------------------------------------------------------------ stepping
    def step(self) -> None:
        r, dt = self.run, self.dt
        if self.phase() == "perfusion" and not self.pid.running:
            self.pid.start()
        stopped = self._pump_stopped()
        self.sensors.offline = stopped
        D_cmd = self.pid.update(self._reading["S_filt"], dt) if not stopped else self.pid.D
        self.D = D_cmd if (self.phase() == "perfusion" and not stopped) else 0.0

        if self.growth_noise > 0:
            self._log_mu_dev += self.growth_noise * np.sqrt(dt) * self.rng.normal()
        mu = self._mu()
        growth = mu * self.X
        self.S = max(0.0, self.S + dt * (self.D * (r.S_I - self.S) - growth / r.Y_xs))
        self.P = max(0.0, self.P + dt * (r.Y_px * growth - self.D * self.P))
        self.C = min(1.0, max(0.0, self.C + dt * ((self.D + r.kLa) * (1.0 - self.C) - growth / r.Y_xc)))
        self.X = self.X + dt * growth
        self.t += dt
        self._reading = self.sensors.read(self.S, self.P, self.C)
        self._record(D_cmd=D_cmd, pump_on=self.D > 0)

    def run_until(self, t_end: float) -> "VirtualFacer":
        while self.t < t_end - 0.5 * self.dt:
            self.step()
        return self

    def run_to_end(self) -> "VirtualFacer":
        self.run_until(self.run.t_F)
        self.X_F_count = self.counts.count(self.X)  # suspension count once perfusion ends
        return self

    # ------------------------------------------------------------------ access
    def series(self, key: str) -> np.ndarray:
        return np.asarray(self.history[key], dtype=object if key == "phase" else float)

    def online_data(self) -> dict[str, np.ndarray]:
        """What the platform logs: time, perfusion command, filtered glucose and lactate, DO."""
        return {"t": self.series("t"), "D": self.series("D"), "S": self.series("S_filt"),
                "P": self.series("P_filt"), "C": self.series("C_meas")}

    def to_csv(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(FIELDS)
            w.writerows(zip(*(self.history[k] for k in FIELDS)))
        return path


def run_published(run: FacerRun, *, seed: int = 0, realistic: bool = True, **kwargs) -> VirtualFacer:
    """Simulate a published run. realistic=False gives noise-free sensors without bubbles."""
    if not realistic:
        kwargs.setdefault("sensors", OnlineSensors(seed=seed, glucose_sd=0.0, lactate_sd=0.0,
                                                   DO_sd=0.0, bubbles=_no_bubbles()))
    else:
        kwargs.setdefault("sensors", OnlineSensors(seed=seed))
    if run.run == 3:
        kwargs.setdefault("pump_stops", [PumpStop(run.t_sed + 2.5, 0.5)])
    return VirtualFacer(run=run, seed=seed, **kwargs).run_to_end()


def _no_bubbles():
    from ..sensors.online import BubbleModel

    return BubbleModel(rate_per_day=0.0)
