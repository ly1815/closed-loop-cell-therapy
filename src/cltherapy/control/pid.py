"""A discrete PI controller with output limits and anti-windup."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class PI:
    kp: float
    ti_h: float  # integral time in hours; inf disables integral action
    setpoint: float
    out_min: float
    out_max: float
    bias: float = 0.0
    reverse: bool = False  # True when a higher output lowers the measurement
    _integral: float = field(default=0.0, init=False)

    def update(self, measurement: float, dt_h: float) -> float:
        error = self.setpoint - measurement
        if self.reverse:
            error = -error
        p = self.kp * error
        self._integral += self.kp * error * dt_h / self.ti_h
        unsaturated = self.bias + p + self._integral
        out = min(max(unsaturated, self.out_min), self.out_max)
        # anti-windup by clamping: keep the integral where the output just reaches its limit
        self._integral -= unsaturated - out
        return out

    def reset(self) -> None:
        self._integral = 0.0
