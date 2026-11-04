# System architecture

The virtual bioreactor, the digital shadow and the supervisory components are arranged on the
levels of the ISA-95 (IEC 62264) model. A digital shadow receives data from the physical component
only (one-way transfer); a digital twin also transfers data back and can adaptively control the
process (Egan et al., 2025). The digital shadow is an advisory tool for a human operator, who
decides the harvest time.

```mermaid
flowchart BT
    subgraph L0["Levels 0-1: process and devices"]
        P["Large chamber of the FACER<br/>batch-perfusion process<br/>(Egan et al. 2025)"]
        R["Virtual bioreactor controller<br/>PID control of D(t) on [S(t)];<br/>OPC UA server (day 4)"]
        A["Off-line cell counter<br/>REST API (day 5)"]
        P --- R
        P --- A
    end
    subgraph L2["Level 2: supervisory control"]
        C["Data collection and orchestration<br/>run state machine (days 5, 7)"]
        H[("Historian and audit trail<br/>(day 9)")]
    end
    subgraph L3["Level 3: manufacturing operations"]
        T["Digital shadow<br/>least squares; ensemble Kalman filter (day 6)"]
        D["Predicted time to [X_T]<br/>(day 8)"]
        B["PID tuning; Bayesian optimisation<br/>across bioreactor runs (day 10)"]
    end
    U["Operator and quality assurance<br/>dashboard"]
    R -->|"OPC UA: D(t), [S(t)], [P(t)], C'(t),<br/>[S_sp], pump state"| C
    A -->|"REST: [X_I], [X_F]"| C
    C --> H
    C --> T
    T --> D
    D -->|"advisory"| U
    U -->|"approved set-point changes"| R
    H --> B
    B -->|"process conditions, PID gains"| U
```

## Control and decision layers

| Layer | Time scale | Implementation | Decision |
| --- | --- | --- | --- |
| PID control of the perfusion rate | 10 s | `control/perfusion.py` | controller, within D_min and the pump limit |
| Digital shadow: nowcasts and short-term forecasts of [X(t)] | hours, per bioreactor run | `shadow/fit.py` | operator and quality assurance (harvest time) |
| Process development across runs | between bioreactor runs | day 10 | process development |

## Implemented components

- `process/`: cell growth model, observed and best-fit parameters, batch-perfusion process model.
- `control/perfusion.py`: PID controller of the perfusion rate; `control/pid.py`: generic PI
  controller used as the example in the control refresher notes.
- `sensors/`: enzymatic glucose and lactate sensors with bubble artefacts and the Butterworth
  filter, DO sensor, off-line cell counts.
- `sim/`: virtual bioreactor run and accelerated clock.
- `shadow/`: digital shadow training and testing, confidence bands, time to a threshold cell
  concentration.
