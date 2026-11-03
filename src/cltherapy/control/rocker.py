"""The rocker's embedded controller: the inner loop, running every few minutes.

Loops:
    pH    PI on the CO2 fraction in the headspace (more CO2 lowers pH)
    DO    PI with split range: rock rate first (gentle for T cells), then O2
    perfusion  a schedule: off, then a ramp in vessel volumes per day
Temperature is assumed perfectly held at 37 C by the tray heater.

The supervisory layer (later: the platform and the digital twin) may change
setpoints and the perfusion schedule, but only within hard limits enforced here.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..process.model import Inputs
from .pid import PI


@dataclass
class Limits:
    rock_rpm: tuple[float, float] = (6.0, 25.0)
    y_O2: tuple[float, float] = (0.21, 0.60)
    y_CO2: tuple[float, float] = (0.01, 0.10)
    perfusion_vvd: tuple[float, float] = (0.0, 2.0)
    pH_setpoint: tuple[float, float] = (6.9, 7.4)
    DO_setpoint: tuple[float, float] = (20.0, 80.0)


@dataclass
class PerfusionSchedule:
    start_day: float = 3.0
    ramp_days: float = 2.0
    max_vvd: float = 1.5

    def rate(self, t_day: float) -> float:
        if t_day < self.start_day:
            return 0.0
        frac = min(1.0, (t_day - self.start_day) / max(self.ramp_days, 1e-6))
        return frac * self.max_vvd


@dataclass
class RockerController:
    pH_setpoint: float = 7.20
    DO_setpoint: float = 50.0  # % air saturation
    rock_angle_deg: float = 6.0
    perfusion: PerfusionSchedule = field(default_factory=PerfusionSchedule)
    limits: Limits = field(default_factory=Limits)

    def __post_init__(self) -> None:
        # pH: higher CO2 lowers pH, so the loop is reverse acting
        self.pH_loop = PI(kp=0.04, ti_h=0.25, setpoint=self.pH_setpoint,
                          out_min=self.limits.y_CO2[0], out_max=self.limits.y_CO2[1],
                          bias=0.05, reverse=True)
        # DO: output 0-100%; 0-50 drives rock rate, 50-100 drives O2 fraction
        self.DO_loop = PI(kp=1.5, ti_h=1.0, setpoint=self.DO_setpoint,
                          out_min=0.0, out_max=100.0, bias=20.0)

    def set_setpoints(self, pH: float | None = None, DO: float | None = None) -> None:
        """Supervisory change, clipped to the hard limits."""
        if pH is not None:
            lo, hi = self.limits.pH_setpoint
            self.pH_setpoint = self.pH_loop.setpoint = min(max(pH, lo), hi)
        if DO is not None:
            lo, hi = self.limits.DO_setpoint
            self.DO_setpoint = self.DO_loop.setpoint = min(max(DO, lo), hi)

    def step(self, t_h: float, pH: float, DO_pct: float, dt_h: float) -> Inputs:
        y_CO2 = self.pH_loop.update(pH, dt_h)
        u = self.DO_loop.update(DO_pct, dt_h)
        rpm_lo, rpm_hi = self.limits.rock_rpm
        o2_lo, o2_hi = self.limits.y_O2
        if u <= 50.0:
            rock = rpm_lo + (rpm_hi - rpm_lo) * u / 50.0
            y_O2 = o2_lo
        else:
            rock = rpm_hi
            y_O2 = o2_lo + (o2_hi - o2_lo) * (u - 50.0) / 50.0
        vvd = self.perfusion.rate(t_h / 24.0)
        vvd = min(max(vvd, self.limits.perfusion_vvd[0]), self.limits.perfusion_vvd[1])
        return Inputs(perfusion_vvd=vvd, rock_rpm=rock, rock_angle_deg=self.rock_angle_deg,
                      y_O2=y_O2, y_CO2=y_CO2)
