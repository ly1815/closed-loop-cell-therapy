# Sensors

What the platform can see, how often, and how well. Noise, drift and timing
are placeholders typical of single-use rocking-bag systems; replace them with
the instrument specifications of a real partner setup.

| Signal | Where it comes from | Update | Noise | Drift | Delay | In code |
| --- | --- | --- | --- | --- | --- | --- |
| pH | optical patch in the bag floor, read through the tray | 10 s | 0.01 | -0.005 per day | none | `sensors/online.py` |
| DO (% air saturation) | optical patch | 10 s | 1% | +0.5% per day | none | `sensors/online.py` |
| Bag weight | load cell in the tray | 1 s | 2 g | none | none | `sensors/online.py` |
| Rock rate and angle, gas fractions, pump rates | controller outputs | each control step | exact | none | none | `control/rocker.py` |
| Viable density, viability | at-line cell counter, from a 2 mL closed sample | on request | 8% CV; 2 points | none | 15 min | `sensors/atline.py` |
| Glucose, lactate | at-line metabolite analyser, same sample | on request | 3% and 4% CV | none | 30 min | `sensors/atline.py` |
| Memory-cell share | flow cytometry at harvest | once | 2 points | none | 12 h | `sensors/atline.py` |
| Oxygen uptake rate | soft sensor: oxygen balance from DO, gas and kLa | each control step | depends on kLa error | inherits DO drift | 30 min window | `sensors/soft.py` |
| Acid production | soft sensor: bicarbonate from pH and CO2 demand | each control step | large, smoothed | inherits pH drift | 6 h window | `sensors/soft.py` |

## Measured, derived or inferred

| State or quantity | How the EnKF will treat it (day 6) |
| --- | --- |
| Total viable cells | measured at-line, sparse; also informed by oxygen uptake |
| Memory-like and effector-like cells | inferred: never measured until flow at harvest |
| Glucose, lactate | measured at-line, sparse; lactate also informed by acid production |
| Bicarbonate | derived from pH and CO2 |
| Dissolved oxygen | measured online |
| Patient parameters (growth, differentiation, activation) | inferred, augmented into the state |

## Things the demo shows on purpose

- **A drifting sensor moves the real process.** The controller holds the
  *reading* at setpoint, so a pH patch drifting down by 0.005 a day lets the
  true pH rise by about the same amount (see `test_controller_follows_the_drifting_sensor`).
- **Soft measurements inherit sensor errors.** Oxygen uptake reads slightly
  low when the DO patch drifts high; acid production needs a long window
  because pH has 0.01 resolution.
- **Every sample costs cells.** `Batch.take_sample` dilutes the bag, as the
  scale-controlled feed refills it with fresh medium, and counts the cells removed.
