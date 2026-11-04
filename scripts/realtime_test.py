"""Pseudo-real-time test of the digital shadow, as in Egan et al. (2025), Figures 5 and 6.

For every run, fit the real-time shadow with 1.5 to 4 days of perfusion data (run 3: up to 2.5
days, before its mechanical issue) and record the relative error of the predicted final cell
concentration. Then plot the 2.5-day forecasts with their 95% bands and the predicted time to
reach 1e7 cells/mL.

    python scripts/realtime_test.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cltherapy.process.parameters import X_T, RUNS
from cltherapy.shadow import OnlineData, fit_realtime
from cltherapy.sim.runner import run_published

OUT = Path(__file__).resolve().parents[1] / "outputs"
COLORS = {1: "#0F6E56", 2: "#2F6DB5", 3: "#111111", 4: "#C4532E", 5: "#7A4FA0", 6: "#B8860B"}
WINDOWS = (1.5, 2.0, 2.5, 3.0, 3.5, 4.0)


def main() -> None:
    OUT.mkdir(exist_ok=True)
    errors = {w: {} for w in WINDOWS}
    fits_25 = {}
    sims = {}
    for n, r in RUNS.items():
        v = run_published(r, seed=n)
        sims[n] = v
        data = OnlineData.from_dict(v.online_data())
        for w in WINDOWS:
            if n == 3 and w > 2.5:
                continue
            f = fit_realtime(r, data, v.X_I_count, r.t_sed + w)
            errors[w][n] = 100 * (float(f.predict_X(r.t_F)) / v.X_F_count - 1)
            if w == 2.5:
                fits_25[n] = f
    lines = ["days_of_data," + ",".join(f"run{n}" for n in RUNS) + ",mean_abs"]
    for w in WINDOWS:
        e = errors[w]
        lines.append(f"{w}," + ",".join(f"{e[n]:.1f}" if n in e else "" for n in RUNS)
                     + f",{np.mean(np.abs(list(e.values()))):.1f}")
    (OUT / "realtime_errors.csv").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))

    fig, ax = plt.subplots(1, 2, figsize=(13, 5))
    ax[0].boxplot([np.abs(list(errors[w].values())) for w in WINDOWS], tick_labels=[str(w) for w in WINDOWS])
    for i, w in enumerate(WINDOWS, start=1):
        for n, e in errors[w].items():
            ax[0].plot(i + 0.15, abs(e), "o", ms=4, color=COLORS[n], alpha=0.8)
    ax[0].axhline(30, color="grey", ls=":", lw=1)
    ax[0].set_xlabel("Online data (days)")
    ax[0].set_ylabel("Relative error of predicted [X$_F$] (%)")
    ax[0].set_title("Digital shadow testing (cf. Egan et al., Fig. 5)")

    for n, f in fits_25.items():
        r, v = RUNS[n], sims[n]
        t = np.linspace(r.t_I, r.t_F, 300)
        b = f.band(t)
        past = t <= f.t_now
        ax[1].plot(t[past], b["X"][past], color=COLORS[n], lw=1.6, label=f"Run {n}")
        ax[1].plot(t[~past], b["X"][~past], color=COLORS[n], lw=1.6, ls="-.")
        ax[1].fill_between(t, b["X_lo"], b["X_hi"], color=COLORS[n], alpha=0.12, lw=0)
        ax[1].plot(r.t_F, v.X_F_count, "D" if "patient" in r.cells else "o", color=COLORS[n], mec="k")
        ax[1].axvline(f.t_now, color=COLORS[n], ls=":", lw=0.8)
    ax[1].axhline(X_T, color="grey", lw=1, ls="--")
    ax[1].text(0.05, X_T + 0.5, "[X$_T$] = 10$^7$ cells/mL", color="grey", fontsize=8)
    ax[1].set_xlabel("Time from start of process (days)")
    ax[1].set_ylabel("[X(t)] (10$^6$ cells/mL)")
    ax[1].set_title("Predictions from 2.5 days of online data, 95% bands (cf. Fig. 6)")
    ax[1].legend(frameon=False, fontsize=8)
    for a in ax:
        a.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(OUT / "realtime_test.png", dpi=150)
    print("Saved", OUT / "realtime_test.png")

    print("\nPredicted time to 1e7 cells/mL from 2.5 days of data (days after seeding):")
    for n, f in fits_25.items():
        r = RUNS[n]
        print(f"  run {n}: {f.threshold_time(X_T) - r.t_I:.2f}")


if __name__ == "__main__":
    main()
