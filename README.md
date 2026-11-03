# Closed-loop cell therapy

A virtual closed-loop CAR-T expansion. One autologous batch grows in a
perfused rocking bag; the rocker's controller holds pH and DO; later steps
add industrial device interfaces, an EnKF soft sensor for the memory-cell
share and the dose, harvest and sampling recommendations an operator
approves, and Bayesian optimisation across process-development runs.

Everything is synthetic and runs offline. Parameters are calibrated to
published CAR-T data where possible; see [docs/calibration.md](docs/calibration.md).

## The question it answers

Longer culture gives more cells but fewer memory cells, and the memory
phenotype is usually measured only at harvest. Can we estimate it during the
run, with uncertainty, and use that to decide when to sample and when to
harvest?

## Quick start

```bash
uv venv && source .venv/bin/activate
uv pip install -e ".[dev]"
pytest
python scripts/run_two_patients.py   # outputs/two_patients.png and CSVs
python scripts/step_response.py      # outputs/step_response.png
python scripts/sensor_demo.py        # outputs/sensors.png: what the platform sees
```

Plain `pip` works too: `python -m venv .venv`, then `pip install -e ".[dev]"`.

## Layout

```
src/cltherapy/
  process/parameters.py  parameters with units, status and sources; patients A and B
  process/model.py       two-population model, rates, observables
  control/pid.py         PI controller with anti-windup
  control/rocker.py      inner loop: pH by CO2, DO by rock rate then O2, perfusion
  sensors/online.py      optical pH and DO patches, load cell: noise, drift, update rate, faults
  sensors/atline.py      cell counter and metabolite analyser with turnaround; flow at harvest
  sensors/soft.py        soft measurements: oxygen uptake and acid production
  sim/runner.py          batch simulation with recorded history and CSV export
  sim/clock.py           accelerated clock for live runs
scripts/                 example runs and plots
tests/                   checks on physics, anchors and control
docs/                    architecture, calibration plan, sensor table
```

## Status

| Day | Step | State |
| --- | --- | --- |
| 1 | Repo, architecture, two-population model | done |
| 2 | Patients A and B, accelerated clock, rocking-bag controller | done (literature fit pending data) |
| 3 | Sensor layer: noise, drift, sampling, delays | done; see [docs/sensors.md](docs/sensors.md) |
| 4 | Rocker as an OPC UA server | next |
| 5 | Analyser as a REST API, data collector | |
| 6 | EnKF soft sensor, sampling by information gain | |
| 7 | Run state machine and scheduler | |
| 8 | Harvest recommendations, first dashboard pages | |
| 9 | Faults, watchdog, audit log | |
| 10 | Patient B page, Bayesian optimisation | |
| 11 | Validation page, packaging, backup video | |
| 12 | Rehearsal | |

## Licence

MIT
