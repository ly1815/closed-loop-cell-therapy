"""At-line analysers: a cell counter and a metabolite analyser.

A sample is drawn through the closed port (Batch.take_sample), queued, and
its result becomes available only after the instrument's turnaround time.
Noise is multiplicative for counts and concentrations. Values are
placeholders; see docs/sensors.md.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import count

import numpy as np


@dataclass
class Instrument:
    name: str
    turnaround_h: float
    cv: dict[str, float]  # measured key -> coefficient of variation
    abs_sd: dict[str, float] = field(default_factory=dict)  # measured key -> absolute sd


CELL_COUNTER = Instrument("cell_counter", turnaround_h=0.25,
                          cv={"viable_density": 0.08}, abs_sd={"viability": 0.02})
METABOLITE_ANALYSER = Instrument("metabolite_analyser", turnaround_h=0.5,
                                 cv={"glucose": 0.03, "lactate": 0.04})


@dataclass
class Result:
    sample_id: int
    instrument: str
    drawn_h: float
    ready_h: float
    values: dict[str, float]


@dataclass
class AtLineLab:
    instruments: tuple[Instrument, ...] = (CELL_COUNTER, METABOLITE_ANALYSER)
    seed: int = 1
    sample_volume_mL: float = 2.0

    def __post_init__(self) -> None:
        self.rng = np.random.default_rng(self.seed)
        self._ids = count(1)
        self._pending: list[Result] = []
        self.delivered: list[Result] = []

    def submit(self, drawn_h: float, drawn_truth: dict) -> int:
        """Queue one drawn sample on every instrument. Returns the sample id."""
        sid = next(self._ids)
        for inst in self.instruments:
            values = {}
            for key, cv in inst.cv.items():
                values[key] = max(0.0, drawn_truth[key] * (1.0 + self.rng.normal(0.0, cv)))
            for key, sd in inst.abs_sd.items():
                values[key] = float(np.clip(drawn_truth[key] + self.rng.normal(0.0, sd), 0.0, 1.0))
            self._pending.append(Result(sid, inst.name, drawn_h, drawn_h + inst.turnaround_h, values))
        return sid

    def poll(self, now_h: float) -> list[Result]:
        """Results that have become ready since the last poll."""
        ready = [r for r in self._pending if r.ready_h <= now_h]
        self._pending = [r for r in self._pending if r.ready_h > now_h]
        self.delivered.extend(ready)
        return ready


@dataclass
class FlowCytometry:
    """Offline phenotype at harvest: the delayed, expensive label."""

    seed: int = 2
    abs_sd: float = 0.02
    turnaround_h: float = 12.0

    def __post_init__(self) -> None:
        self.rng = np.random.default_rng(self.seed)

    def measure(self, drawn_h: float, drawn_truth: dict) -> Result:
        mf = float(np.clip(drawn_truth["memory_fraction"] + self.rng.normal(0.0, self.abs_sd), 0.0, 1.0))
        return Result(0, "flow_cytometry", drawn_h, drawn_h + self.turnaround_h, {"memory_fraction": mf})
