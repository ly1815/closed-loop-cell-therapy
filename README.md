# Closed-loop cell therapy

A virtual CAR-T expansion in a glucose-controlled perfusion bioreactor, built entirely on

> Egan JR, Marí-Buyé N, Benítez-Cano EV, et al. *A digital shadow of CAR T cell expansion in a
> perfusion bioreactor: informing optimal harvest times for autologous cell therapy.*
> Biotechnology Progress 2025;41(6):e70045. [doi:10.1002/btpr.70045](https://doi.org/10.1002/btpr.70045)
> (open access, CC BY 4.0)

The process is the Aglaris FACER: T cells sediment onto a gas-permeable membrane in a
134 mL chamber, then fresh medium is perfused at a rate set by a PID controller that holds glucose
at a set point. Cell concentration cannot be measured during the run without resuspending the
cells, so the only online signals are the perfusion rate, glucose, lactate and dissolved oxygen.
The digital shadow fits a mechanistic model to those signals and forecasts the cell concentration,
and with it the time the batch reaches the required dose.

Everything runs offline on synthetic data generated from the published parameters of the six
FACER runs (Tables 1, 2 and S4 of the paper). The paper's raw data are proprietary.

## What is in the repository

| Layer | Module | Source in Egan et al. |
| --- | --- | --- |
| Growth model: Richards special case with a lag | `process/growth.py` | eqs 10-13, S60-S70 |
| The six published runs, with sources | `process/parameters.py` | Tables 1, 2, S4 |
| Mass balances, forward Euler at 10 s | `process/facer.py` | eqs 2-8, S19-S24 |
| Glucose PID setting the perfusion rate | `control/perfusion.py` | eqs 2, S12, S13, S23 |
| Enzymatic sensors, bubbles, Butterworth filter, offline counts | `sensors/online.py` | section 2.2 |
| Virtual FACER run, with the run 3 pump stop | `sim/runner.py` | sections 2.2, 3.1 |
| Digital shadow: training and real-time fits, bands, harvest time | `shadow/fit.py` | sections 2.3, 3.2, 3.3; S71, S72 |
| Accelerated clock for live runs | `sim/clock.py` | |
| Generic PI with anti-windup (control notes example) | `control/pid.py` | |

## Quick start

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest
python scripts/simulate_run.py 1     # outputs/run1.png, laid out like the paper's Figure 3
python scripts/simulate_run.py all
python scripts/realtime_test.py      # outputs/realtime_test.png, cf. Figures 5 and 6
```

## What the virtual runs reproduce

- Final cell concentrations of all six runs, and Table 2's maximum growth rates and expansion
  periods to 1e7 cells/mL (tests in `tests/test_growth.py`, `tests/test_facer.py`).
- Glucose held at its set point and lactate settling at Y_ps (S_I - S_sp), eq. 18.
- Bubble artefacts: downward glucose and lactate spikes and upward perfusion spikes after filtering.
- Run 3's half-day pump stop 2.5 days into perfusion, with the online data lost.
- The real-time shadow with 2.5 days of data predicting the final concentration within 30% for
  every run, mean absolute error about 8% here against 13% in the paper (`scripts/realtime_test.py`).
- Confidence bands that are too narrow, as the paper found: they ignore the uncertainty in the
  fixed Y_xs and t_lag.

## Documentation

- [docs/model.md](docs/model.md): equations, parameters and where each number comes from
- [docs/sensors.md](docs/sensors.md): what the platform measures, how often, and how well
- [docs/architecture.md](docs/architecture.md): devices, protocols and the loops
- [docs/technical_documentation/](docs/technical_documentation/): full LaTeX write-up (PDF included)
- [docs/control_refresher/](docs/control_refresher/): control theory refresher notes

The paper and its supporting information go in `Egan/` locally (ignored by git).

## Plan

| Day | Step | State |
| --- | --- | --- |
| 1-3 | FACER model, glucose PID, sensors and filtering, virtual runs, digital shadow | done |
| 4 | Virtual FACER as an OPC UA server (D, S, P, DO, set point) | next |
| 5 | Offline counter as a REST API, data collector | |
| 6 | EnKF version of the shadow, sampling Y_xs and t_lag (the paper's suggested improvement) | |
| 7 | Run state machine: seed, sediment, perfuse, harvest | |
| 8 | Harvest-time recommendations, first dashboard pages | |
| 9 | Faults (pump stop, sensor bubbles, drift), watchdog, audit log | |
| 10 | PID tuning with the shadow; Bayesian optimisation across runs | |
| 11 | Validation page, packaging, backup video | |
| 12 | Rehearsal | |

## Licence

MIT for the code. Parameter values are taken from Egan et al. (2025), CC BY 4.0.
