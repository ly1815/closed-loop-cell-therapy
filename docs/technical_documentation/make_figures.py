"""Figures for the technical documentation. Run from the repository root:

    python docs/technical_documentation/make_figures.py
"""

from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "scripts"))

from cltherapy.process.growth import cell_concentration, specific_growth_rate  # noqa: E402
from cltherapy.process.parameters import RUNS  # noqa: E402
from cltherapy.sim.runner import run_published  # noqa: E402

FIG = Path(__file__).resolve().parent / "figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.3, "legend.frameon": False})
BLUE, RED, GREEN, GREY = "#4472C4", "#E0533C", "#0F6E56", "#777777"


def growth_family() -> None:
    t = np.linspace(0, 5, 400)
    fig, ax = plt.subplots(1, 2, figsize=(8, 3))
    for nu, c in ((-1.0, RED), (-0.5, BLUE), (-0.25, "#7A4FA0"), (0.0, GREEN)):
        label = "exponential" if nu == 0 else f"$\\nu$ = {nu:g}"
        ax[0].plot(t, cell_concentration(t, 1.0, 0.2, 1.5, nu), color=c, label=label)
        ax[1].plot(t, specific_growth_rate(t, 0.2, 1.5, nu), color=c, label=label)
    ax[0].set_ylim(0, 25)
    ax[0].set_ylabel("X / X$_I$")
    ax[1].set_ylabel("$\\mu$ (day$^{-1}$)")
    for a in ax:
        a.set_xlabel("Time after seeding (days)")
    ax[0].legend()
    fig.tight_layout()
    fig.savefig(FIG / "growth_family.png", dpi=200)
    plt.close(fig)


def filter_zoom() -> None:
    r = RUNS[1]
    v = run_published(r, seed=1)
    t = v.series("t")
    raw, filt, D = v.series("S_raw"), v.series("S_filt"), v.series("D")
    bubble = np.argmax((t > r.t_sed + 1.0) & (raw < 0.8 * v.series("S")))
    t0 = t[bubble]
    k = (t > t0 - 0.04) & (t < t0 + 0.1)
    fig, ax = plt.subplots(2, 1, figsize=(8, 4), sharex=True)
    ax[0].plot((t[k] - t0) * 24 * 60, raw[k], color="#B5C7E8", lw=0.8, label="Raw enzymatic reading")
    ax[0].plot((t[k] - t0) * 24 * 60, filt[k], color=BLUE, lw=1.5, label="After Butterworth filter")
    ax[0].plot((t[k] - t0) * 24 * 60, v.series("S")[k], color=RED, lw=1.0, ls="--", label="True glucose")
    ax[0].set_ylabel("Glucose (mmol/L)")
    ax[0].legend(fontsize=8)
    ax[1].plot((t[k] - t0) * 24 * 60, D[k], color=GREEN)
    ax[1].set_ylabel("Perfusion (day$^{-1}$)")
    ax[1].set_xlabel("Minutes from the bubble")
    fig.tight_layout()
    fig.savefig(FIG / "bubble_filter.png", dpi=200)
    plt.close(fig)


def copy_script_figures() -> None:
    import realtime_test
    import simulate_run

    out = ROOT / "outputs"
    for n in (1, 3):
        simulate_run.plot_run(n)
        (FIG / f"run{n}.png").write_bytes((out / f"run{n}.png").read_bytes())
    realtime_test.main()
    (FIG / "realtime_test.png").write_bytes((out / "realtime_test.png").read_bytes())
    (FIG / "realtime_errors.csv").write_bytes((out / "realtime_errors.csv").read_bytes())


if __name__ == "__main__":
    growth_family()
    filter_zoom()
    copy_script_figures()
    print("Figures in", FIG)
