"""Simulate patients A and B for 10 days and save plots and CSVs to outputs/.

    python scripts/run_two_patients.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from cltherapy.process.parameters import MEMORY_ANCHORS, PATIENT_A, PATIENT_B
from cltherapy.sim.runner import run_batch

OUT = Path(__file__).resolve().parents[1] / "outputs"
COLORS = {"A": "#0F6E56", "B": "#C4532E"}


def main() -> None:
    OUT.mkdir(exist_ok=True)
    runs = {"A": run_batch(PATIENT_A, days=10.0), "B": run_batch(PATIENT_B, days=10.0)}
    for key, b in runs.items():
        b.to_csv(OUT / f"patient_{key}.csv")

    panels = [
        ("viable_density", "Viable cells (1e6/mL)"),
        ("memory_fraction", "Memory-like share (hidden)"),
        ("viability", "Viability"),
        ("glucose", "Glucose (mM)"),
        ("lactate", "Lactate (mM)"),
        ("OUR", "Oxygen uptake (mM/h)"),
        ("pH", "pH"),
        ("DO_pct", "DO (% air sat.)"),
        ("y_CO2", "CO2 in gas (fraction)"),
        ("rock_rpm", "Rock rate (rpm)"),
        ("y_O2", "O2 in gas (fraction)"),
        ("perfusion_vvd", "Perfusion (VVD)"),
    ]
    fig, axes = plt.subplots(4, 3, figsize=(13, 12), sharex=True)
    for ax, (key, label) in zip(axes.flat, panels):
        for name, b in runs.items():
            ax.plot(b.series("t_day"), b.series(key), color=COLORS[name], lw=1.6, label=f"Patient {name}")
        ax.set_ylabel(label)
        ax.grid(alpha=0.25)
        if key == "memory_fraction":
            ax.scatter(list(MEMORY_ANCHORS), list(MEMORY_ANCHORS.values()), color="black", zorder=3,
                       label="Cadinanos-Garai 2025")
    axes[0, 1].legend(frameon=False, fontsize=9)
    for ax in axes[-1]:
        ax.set_xlabel("Culture day")
    fig.suptitle("Virtual CAR-T expansion in a perfused rocking bag (synthetic)", y=0.995)
    fig.tight_layout()
    fig.savefig(OUT / "two_patients.png", dpi=150)

    for name, b in runs.items():
        r9 = b.at_day(9)
        print(f"Patient {name}, day 9: {r9['total_viable_cells'] / 1e3:.1f}e9 viable cells, "
              f"memory share {r9['memory_fraction']:.1%}, viability {r9['viability']:.1%}")
    print(f"Saved plots and CSVs to {OUT}")


if __name__ == "__main__":
    main()
