"""Simulate one FACER run and plot it in the layout of Egan et al. (2025), Figure 3.

    python scripts/simulate_run.py          # run 1
    python scripts/simulate_run.py 3        # run 3, including its mechanical issue (approx. 0.5 days)
    python scripts/simulate_run.py all
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cltherapy.process.facer import simulate_model
from cltherapy.process.parameters import RUNS
from cltherapy.sim.runner import run_published

OUT = Path(__file__).resolve().parents[1] / "outputs"
BLUE, RED, SED = "#4472C4", "#E0533C", "#DDEFE6"


def plot_run(n: int) -> Path:
    r = RUNS[n]
    v = run_published(r, seed=n)
    m = simulate_model(r)
    t = v.series("t")
    panels = [
        ("D", "D", "(a) Perfusion rate, D(t)", "day$^{-1}$"),
        ("S_filt", "S", "(b) Glucose concentration, [S(t)]", "mmol/L"),
        ("P_filt", "P", "(c) Lactate concentration, [P(t)]", "mmol/L"),
        ("C_meas", "C", "(d) Normalised dissolved oxygen, [C(t)]/[C$^*$]", "unitless"),
        (None, "X", "(e) Cell concentration, [X(t)]", "10$^6$ cells/mL"),
        (None, "mu", "(f) Specific growth rate, $\\mu(t)$", "day$^{-1}$"),
    ]
    fig, axes = plt.subplots(len(panels), 1, figsize=(8, 11), sharex=True)
    for ax, (online, model, title, unit) in zip(axes, panels):
        ax.axvspan(r.t_I, r.t_sed, color=SED, zorder=0, label="Sedimentation period")
        if online is not None:
            ax.plot(t, v.series(online), color=BLUE, lw=0.8, label="On-line data (simulated)")
        ax.plot(m.t, getattr(m, model), color=RED, lw=1.4, label="Model best-fit (published parameters)")
        if model == "X":
            ax.plot([r.t_I, r.t_F], [v.X_I_count, v.X_F_count], "o", color=BLUE, mec="k",
                    label="Off-line data")
        ax.set_title(title, fontsize=10)
        ax.set_ylabel(unit)
        ax.grid(alpha=0.3)
    axes[0].legend(frameon=False, fontsize=8, ncol=3, loc="upper left")
    axes[4].legend(frameon=False, fontsize=8, loc="upper left")
    axes[-1].set_xlabel("Time from start of process (days)")
    note = " (derived parameters, mechanical issue at day %.2f)" % (r.t_sed + 2.5) if r.derived else ""
    fig.suptitle(f"Virtual FACER, run {n}: {r.cells}, media {r.media}{note}", fontsize=11)
    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    path = OUT / f"run{n}.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    v.to_csv(OUT / f"run{n}.csv")
    return path


def main() -> None:
    arg = sys.argv[1] if len(sys.argv) > 1 else "1"
    runs = list(RUNS) if arg == "all" else [int(arg)]
    for n in runs:
        print("Saved", plot_run(n))


if __name__ == "__main__":
    main()
