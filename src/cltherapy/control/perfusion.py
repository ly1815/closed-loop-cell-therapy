"""PID control of the perfusion rate, with the glucose concentration as the process variable.

Velocity form of eqs S12 and S13 of Egan et al. (2025), discretised as in eq. S23, acting on
the filtered online glucose data rather than the true [S(t)]:

    e_k   = [S_sp] - [S_meas](t_k)
    D_k+1 = max(D_min, D_k + K_p (e_k+1 - e_k) + K_i dt e_k + K_d (e_k+1 - 2 e_k + e_k-1) / dt)

With K_d = 0 this is exactly S23. The derivative gain is 3-4 orders of magnitude below K_p and
K_i (Table 1) so that bubble artefacts do not kick the pump. During sedimentation the pump is
off; perfusion starts at D_min when sedimentation ends.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from ..process.parameters import ControllerGains


@dataclass
class GlucosePerfusionPID:
    S_sp: float
    gains: ControllerGains = field(default_factory=ControllerGains)
    D_max: float = 20.0  # 1/day, pump limit
    use_derivative: bool = True

    def __post_init__(self) -> None:
        self.D = 0.0
        self.running = False
        self._e1: float | None = None  # previous error
        self._e2: float | None = None  # error before that

    def start(self) -> None:
        """End of sedimentation: start perfusing at the minimum rate."""
        self.running = True
        self.D = self.gains.D_min

    def update(self, S_meas: float | None, dt: float) -> float:
        """New perfusion rate (1/day) from the latest filtered glucose reading (mmol/L)."""
        if not self.running:
            return 0.0
        if S_meas is None:  # no reading: hold the last output
            return self.D
        e = self.S_sp - S_meas
        if self._e1 is not None:
            g = self.gains
            dD = g.kp * (e - self._e1) + g.ki * dt * self._e1
            if self.use_derivative and self._e2 is not None:
                dD += g.kd * (e - 2.0 * self._e1 + self._e2) / dt
            self.D = min(self.D_max, max(g.D_min, self.D + dD))
        self._e2, self._e1 = self._e1, e
        return self.D

    def set_setpoint(self, S_sp: float, limits: tuple[float, float] = (10.0, 25.0)) -> float:
        """Supervisory change of the glucose set point, clipped to hard limits."""
        self.S_sp = min(max(S_sp, limits[0]), limits[1])
        return self.S_sp
