"""What the platform actually sees: noisy online signals, daily at-line samples,
soft measurements, and flow cytometry only at harvest.

    python scripts/sensor_demo.py
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from cltherapy.process.parameters import PATIENT_A
from cltherapy.sensors import AtLineLab, FlowCytometry, OnlineSensors, SoftSensors
from cltherapy.sim.runner import Batch

OUT = Path(__file__).resolve().parents[1] / "outputs"
SAMPLE_HOURS = [24 * d + 9 for d in range(1, 9)]  # one closed sample a day at 09:00
HARVEST_H = 9 * 24


def main() -> None:
    OUT.mkdir(exist_ok=True)
    batch = Batch(patient=PATIENT_A, sensors=OnlineSensors(seed=0))
    lab = AtLineLab(seed=1)
    results = []
    while batch.t_h < HARVEST_H - 1e-9:
        batch.step()
        if any(abs(batch.t_h - h) < batch.control_dt_h / 2 for h in SAMPLE_HOURS):
            lab.submit(batch.t_h, batch.take_sample(lab.sample_volume_mL))
        results.extend(lab.poll(batch.t_h))
    flow = FlowCytometry(seed=2).measure(batch.t_h, batch.truth())

    h = batch.history[1:]
    t = np.array([r["t_h"] for r in h])
    soft = SoftSensors().estimate(
        t, [r["meas_DO_pct"] for r in h], [r["meas_pH"] for r in h], [r["y_O2"] for r in h],
        [r["y_CO2"] for r in h], [r["rock_rpm"] for r in h], [r["rock_angle_deg"] for r in h],
        [r["perfusion_vvd"] for r in h],
    )
    day = t / 24.0
    cc = [r for r in results if r.instrument == "cell_counter"]
    ma = [r for r in results if r.instrument == "metabolite_analyser"]

    fig, ax = plt.subplots(3, 2, figsize=(12, 10), sharex=True)
    ax[0, 0].plot(day, [r["meas_pH"] for r in h], color="#9FC9BB", lw=0.6, label="Optical patch")
    ax[0, 0].plot(day, [r["pH"] for r in h], color="#0F6E56", lw=1.5, label="True")
    ax[0, 0].set_ylabel("pH")
    ax[0, 1].plot(day, [r["meas_DO_pct"] for r in h], color="#B5D4F4", lw=0.6, label="Optical patch")
    ax[0, 1].plot(day, [r["DO_pct"] for r in h], color="#185FA5", lw=1.5, label="True")
    ax[0, 1].set_ylabel("DO (% air sat.)")
    ax[1, 0].plot(day, [r["viable_density"] for r in h], color="#0F6E56", label="True")
    ax[1, 0].errorbar([r.ready_h / 24 for r in cc], [r.values["viable_density"] for r in cc],
                      fmt="o", color="#378ADD", label="Cell counter (ready time)")
    ax[1, 0].set_ylabel("Viable cells (1e6/mL)")
    ax[1, 1].plot(day, [r["glucose"] for r in h], color="#0F6E56", label="Glucose, true")
    ax[1, 1].plot(day, [r["lactate"] for r in h], color="#C4532E", label="Lactate, true")
    ax[1, 1].plot([r.ready_h / 24 for r in ma], [r.values["glucose"] for r in ma], "o", color="#0F6E56")
    ax[1, 1].plot([r.ready_h / 24 for r in ma], [r.values["lactate"] for r in ma], "o", color="#C4532E")
    ax[1, 1].set_ylabel("mM")
    ax[2, 0].plot(day, soft["OUR_hat"], color="#B5D4F4", lw=0.8, label="Soft sensor")
    ax[2, 0].plot(day, [r["OUR"] for r in h], color="#185FA5", lw=1.5, label="True")
    ax[2, 0].set_ylabel("Oxygen uptake (mM/h)")
    ax[2, 1].plot(day, soft["acid_hat"], color="#F5C4B3", lw=0.8, label="Soft sensor")
    ax[2, 1].plot(day, [r["acid_production"] for r in h], color="#C4532E", lw=1.5, label="True")
    ax[2, 1].set_ylabel("Acid production (mM/h)")
    for a in ax.flat:
        a.grid(alpha=0.25)
        a.legend(frameon=False, fontsize=8)
    for a in ax[-1]:
        a.set_xlabel("Culture day")
    fig.suptitle(
        f"What the platform sees (patient A). Flow cytometry at harvest: memory share "
        f"{flow.values['memory_fraction']:.1%}, ready {flow.ready_h / 24:.1f} d; "
        f"{batch.n_samples} samples took {batch.cells_sampled:.0f}e6 cells"
    )
    fig.tight_layout()
    fig.savefig(OUT / "sensors.png", dpi=150)
    print(f"Saved {OUT / 'sensors.png'}")


if __name__ == "__main__":
    main()
