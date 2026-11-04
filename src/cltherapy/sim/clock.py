"""An accelerated clock: simulated hours per wall-clock second.

At speed 1.0, one simulated hour passes per second, so a 6-day FACER run
takes 2.4 minutes. Used later by the live devices and the dashboard.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field


@dataclass
class SimClock:
    speed_h_per_s: float = 1.0
    _start_wall: float = field(default_factory=time.monotonic, init=False)
    _offset_h: float = field(default=0.0, init=False)
    _paused_at: float | None = field(default=None, init=False)

    def now_h(self) -> float:
        wall = self._paused_at if self._paused_at is not None else time.monotonic()
        return self._offset_h + (wall - self._start_wall) * self.speed_h_per_s

    def pause(self) -> None:
        if self._paused_at is None:
            self._paused_at = time.monotonic()

    def resume(self) -> None:
        if self._paused_at is not None:
            self._start_wall += time.monotonic() - self._paused_at
            self._paused_at = None

    def set_speed(self, speed_h_per_s: float) -> None:
        current = self.now_h()
        self._start_wall = time.monotonic()
        if self._paused_at is not None:
            self._paused_at = self._start_wall
        self._offset_h = current
        self.speed_h_per_s = speed_h_per_s

    def jump_to(self, t_h: float) -> None:
        self._offset_h = t_h
        self._start_wall = time.monotonic()
        if self._paused_at is not None:
            self._paused_at = self._start_wall
