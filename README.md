# Closed-loop cell therapy

A virtual CAR-T expansion in a glucose-controlled perfusion bioreactor, built entirely on

> Egan JR, Marí-Buyé N, Benítez-Cano EV, et al. *A digital shadow of CAR T cell expansion in a
> perfusion bioreactor: informing optimal harvest times for autologous cell therapy.*
> Biotechnology Progress 2025;41(6):e70045. [doi:10.1002/btpr.70045](https://doi.org/10.1002/btpr.70045)
> (open access, CC BY 4.0)

The repository re-implements the digital shadow of cell growth in a PID controlled perfusion
bioreactor, the Aglaris FACER, in which cell retention is via sedimentation and the PID controller
maintains the glucose concentration [S(t)] at its set-point [S_sp] by adjusting the perfusion rate
D(t). The digital shadow is fitted to online data of D(t), [S(t)] and [P(t)] and the off-line initial
cell concentration [X_I], and provides nowcasts and short-term forecasts of the cell concentration
[X(t)], and hence of the time t_T at which a threshold cell concentration [X_T] is reached.

Since the original datasets are proprietary, a virtual bioreactor generates online and off-line data
for the six bioreactor runs from the published observed and best-fit parameters (Tables 1, 2 and S4).

## What is in the repository

| Layer | Module | Source in Egan et al. |
| --- | --- | --- |
| Cell growth model: special case of the Richards growth model with a lag period | `process/growth.py` | eqs 10-13, S60-S70 |
| Observed and best-fit parameters of the six bioreactor runs | `process/parameters.py` | Tables 1, 2, S4 |
| Batch-perfusion process model, forward Euler method, 10 s step | `process/facer.py` | eqs 2-8, S19-S24 |
| PID control of the perfusion rate | `control/perfusion.py` | eqs 2, S12, S13, S23 |
| Enzymatic sensors, bubble artefacts, Butterworth filter, off-line cell counts | `sensors/online.py` | section 2.2 |
| Virtual bioreactor run, including the mechanical issue of run 3 | `sim/runner.py` | sections 2.2, 3.1 |
| Digital shadow training and testing, confidence bands, t_T | `shadow/fit.py` | sections 2.3, 3.2, 3.3; S71, S72 |
| Accelerated clock for live runs | `sim/clock.py` | |
| Generic PI controller with anti-windup (control refresher notes) | `control/pid.py` | |

## Quick start

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest
python scripts/simulate_run.py 1     # outputs/run1.png, laid out like the paper's Figure 3
python scripts/simulate_run.py all
python scripts/realtime_test.py      # outputs/realtime_test.png, cf. Figures 5 and 6
```

## Results

- The final cell concentrations [X_F] of all six runs, and the maximum specific growth rates mu_max
  and required expansion periods t_T - t_I for [X_T] = 10^7 cells/mL of Table 2, are reproduced
  (`tests/test_growth.py`, `tests/test_facer.py`).
- The glucose concentration remains at its set-point and the lactate concentration approximately
  stabilises at the steady state [P_ss] = Y_ps ([S_I] - [S_sp]) (eq. 18).
- Bubble artefacts not removed by the filter give downward spikes of the glucose and lactate
  concentrations followed by upward spikes in the perfusion rate.
- The mechanical issue of run 3 (no perfusion and no online data for approximately 0.5 days,
  2.5 days after the start of perfusion) is reproduced.
- Digital shadow testing: with 2.5 days of online data, [X_F] is predicted within a 30% relative
  error for all runs, with a mean relative error of 8% (Egan et al.: 13%) (`scripts/realtime_test.py`).
- As in the original study, the confidence bands are too narrow, since the uncertainty in the fixed
  parameters Y_xs and t_lag is not captured.

## Documentation

- [docs/model.md](docs/model.md): mathematical model, parameters and their sources
- [docs/sensors.md](docs/sensors.md): online and off-line data
- [docs/architecture.md](docs/architecture.md): system architecture
- [docs/technical_documentation/](docs/technical_documentation/): technical documentation (LaTeX source and PDF)
- [docs/control_refresher/](docs/control_refresher/): control theory refresher notes

The paper and its supporting information go in `Egan/` locally (ignored by git).

## Plan

| Day | Step | State |
| --- | --- | --- |
| 1-3 | Batch-perfusion process model, PID control, online and off-line data, virtual bioreactor runs, digital shadow | done |
| 4 | Virtual bioreactor as an OPC UA server (D(t), [S(t)], [P(t)], C'(t), [S_sp]) | next |
| 5 | Off-line cell counter as a REST API, data collection | |
| 6 | Ensemble Kalman filter digital shadow with uncertain Y_xs and t_lag | |
| 7 | Run state machine: seeding, sedimentation period, perfusion period, harvest | |
| 8 | Decision support for the harvest time, dashboard | |
| 9 | Process disturbances, data validation, audit trail | |
| 10 | PID tuning via the digital shadow; Bayesian optimisation across bioreactor runs | |
| 11 | Validation page, packaging, backup video | |
| 12 | Rehearsal | |

## Licence

MIT for the code. Parameter values are taken from Egan et al. (2025), CC BY 4.0.
