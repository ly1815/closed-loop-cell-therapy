"""Step-response check of the inner loop: change the pH and DO setpoints mid-run.

    python scripts/step_response.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cltherapy.process.parameters import PATIENT_A
from cltherapy.sim.runner import Batch

OUT = Path(__file__).resolve().parents[1] / "outputs"


def main() -> None:
    OUT.mkdir(exist_ok=True)
    b = Batch(patient=PATIENT_A)
    b.run(7.0)
    b.controller.set_setpoints(pH=7.10, DO=40.0)  # a supervisory change on day 7
    b.run(1.0)

    t = b.series("t_day")
    fig, axes = plt.subplots(2, 2, figsize=(11, 6), sharex=True)
    axes[0, 0].plot(t, b.series("pH"), color="#0F6E56")
    axes[0, 0].plot(t, b.series("pH_setpoint"), "k--", lw=1)
    axes[0, 0].set_ylabel("pH")
    axes[1, 0].plot(t, b.series("y_CO2"), color="#0F6E56")
    axes[1, 0].set_ylabel("CO2 in gas (fraction)")
    axes[0, 1].plot(t, b.series("DO_pct"), color="#378ADD")
    axes[0, 1].plot(t, b.series("DO_setpoint"), "k--", lw=1)
    axes[0, 1].set_ylabel("DO (% air sat.)")
    axes[1, 1].plot(t, b.series("rock_rpm"), color="#378ADD")
    axes[1, 1].set_ylabel("Rock rate (rpm)")
    for ax in axes.flat:
        ax.set_xlim(6.5, 8.0)
        ax.grid(alpha=0.25)
    for ax in axes[-1]:
        ax.set_xlabel("Culture day")
    fig.suptitle("Inner loop: setpoint change on day 7 (pH 7.20 to 7.10, DO 50% to 40%)")
    fig.tight_layout()
    fig.savefig(OUT / "step_response.png", dpi=150)
    print(f"Saved {OUT / 'step_response.png'}")


if __name__ == "__main__":
    main()
