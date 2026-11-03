"""Run one batch: process model plus the rocker's inner loop.

The controller acts every `control_dt_h`; between actions the inputs are held
and the ODEs are integrated with scipy. Every interval is recorded, so a run
can be replayed, plotted or streamed to devices later.
"""

from __future__ import annotations

import csv
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

from ..control.rocker import RockerController
from ..process.model import Inputs, initial_state, observables, rhs, STATE_NAMES
from ..process.parameters import DEFAULT_PROCESS, PatientParams, ProcessParams


@dataclass
class Batch:
    patient: PatientParams
    process: ProcessParams = DEFAULT_PROCESS
    controller: RockerController = field(default_factory=RockerController)
    control_dt_h: float = 5.0 / 60.0
    t_h: float = 0.0
    history: list[dict] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.x = initial_state(self.process, self.patient)
        self.u = Inputs(rock_angle_deg=self.controller.rock_angle_deg)
        self._record()

    def _record(self) -> None:
        obs = observables(self.t_h, self.x, self.u, self.process, self.patient)
        row = {"t_h": self.t_h, "t_day": self.t_h / 24.0}
        row.update({f"x_{n}": float(v) for n, v in zip(STATE_NAMES, self.x)})
        row.update({k: float(v) for k, v in obs.items()})
        row.update({
            "perfusion_vvd": self.u.perfusion_vvd,
            "rock_rpm": self.u.rock_rpm,
            "y_O2": self.u.y_O2,
            "y_CO2": self.u.y_CO2,
            "pH_setpoint": self.controller.pH_setpoint,
            "DO_setpoint": self.controller.DO_setpoint,
        })
        self.history.append(row)

    def step(self) -> dict:
        """Advance one control interval and return the newest record."""
        obs = observables(self.t_h, self.x, self.u, self.process, self.patient)
        self.u = self.controller.step(self.t_h, obs["pH"], obs["DO_pct"], self.control_dt_h)
        sol = solve_ivp(
            rhs, (self.t_h, self.t_h + self.control_dt_h), self.x,
            args=(self.u, self.process, self.patient), method="LSODA", rtol=1e-6, atol=1e-9,
        )
        self.x = np.maximum(sol.y[:, -1], 0.0)
        self.t_h += self.control_dt_h
        self._record()
        return self.history[-1]

    def run(self, days: float) -> "Batch":
        n = int(round(days * 24.0 / self.control_dt_h))
        for _ in range(n):
            self.step()
        return self

    def series(self, key: str) -> np.ndarray:
        return np.array([row[key] for row in self.history])

    def at_day(self, day: float) -> dict:
        i = int(np.argmin(np.abs(self.series("t_day") - day)))
        return self.history[i]

    def to_csv(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        with path.open("w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=list(self.history[0].keys()))
            writer.writeheader()
            writer.writerows(self.history)
        return path


def run_batch(patient: PatientParams, days: float = 10.0, **kwargs) -> Batch:
    return Batch(patient=patient, **kwargs).run(days)
