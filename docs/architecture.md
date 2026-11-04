# Architecture

The demo follows one autologous CAR-T batch through a FACER perfusion expansion. Each layer is
its own component, placed on the ISA-95 levels, so the software around the model looks like a
real deployment.

```mermaid
flowchart BT
    subgraph L0["Level 0-1: process and devices"]
        P["FACER large chamber<br/>sedimented T cells, glucose,<br/>lactate, DO (Egan et al. 2025)"]
        R["Virtual FACER controller<br/>glucose PID on the pump;<br/>OPC UA server (day 4)"]
        A["Offline cell counter<br/>REST API (day 5)"]
        P --- R
        P --- A
    end
    subgraph L2["Level 2: supervision"]
        C["Collector and orchestrator<br/>run state machine (days 5, 7)"]
        H[("Historian and audit log<br/>SQLite (day 9)")]
    end
    subgraph L3["Level 3: operations"]
        T["Digital shadow<br/>least squares now, EnKF (day 6)"]
        D["Harvest-time recommendation<br/>(day 8)"]
        B["PID tuning and Bayesian<br/>optimisation across runs (day 10)"]
    end
    U["Operator and QA<br/>Streamlit dashboard, thin client"]
    R -->|"OPC UA: D, S, P, DO,<br/>set point, pump state"| C
    A -->|"REST: counts at seeding and end"| C
    C --> H
    C --> T
    T --> D
    D -->|"advisory, needs approval"| U
    U -->|"approved set-point changes, within limits"| R
    H --> B
    B -->|"next run conditions, PID gains"| U
```

## The three loops

| Loop | Timescale | Where it lives | Decides |
| --- | --- | --- | --- |
| Inner | 10 s | `control/perfusion.py`: PID on filtered glucose, setting the perfusion rate | the controller, within D_min and the pump limit |
| Run | hours, one patient | the digital shadow and its harvest forecast | operator and QA |
| Campaign | between runs | PID tuning and Bayesian optimisation (day 10) | the process development team |

A digital shadow only receives data (physical to virtual); a digital twin also acts back on the
process. Egan et al. argue the shadow's advisory role keeps regulatory risk low: it predicts, a
person decides.

## What exists now

- `process/`: growth model, published runs, mass balances.
- `control/perfusion.py`: the glucose PID; `control/pid.py`: a generic PI kept for the control notes.
- `sensors/`: enzymatic glucose and lactate with bubbles and filtering, DO, offline counts.
- `sim/`: the virtual FACER run and the accelerated clock.
- `shadow/`: training and real-time fits, confidence bands, harvest time.
