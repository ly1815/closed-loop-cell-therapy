# Sensors

What the FACER platform measures, following Egan et al. (2025), section 2.2. Noise levels and
bubble statistics are placeholders chosen to resemble the published traces; the rest is as
published.

| Signal | How | Rate | In the simulator |
| --- | --- | --- | --- |
| Glucose | single-use enzymatic sensor in the cartridge | every 10 s | noise 0.08 mmol/L, bubbles, then the Butterworth filter |
| Lactate | single-use enzymatic sensor | every 10 s | same as glucose |
| Dissolved oxygen | online sensor, reported as % of saturation | online | normalised C/C*, noise 0.01 |
| pH | online sensor | online | not modelled (not part of the paper's model) |
| Perfusion rate | the PID's own output | every 10 s | exact |
| Cell concentration | offline counter, cells resuspended | at seeding and after perfusion only | 5% CV |
| Phenotype | flow cytometry; can use sedimented cells, but gives no count | offline | not modelled |

## Bubbles and filtering

Bubbles passing the sensors cause artificial drops in both glucose and lactate, followed by a
recovery when the bubble moves on. To stop them driving the pump, the raw signals go through a
2nd-order low-pass Butterworth filter (cut-off 0.001 Hz) before the PID. Some artefacts survive
the filter, which is why the published data show downward glucose spikes followed by upward
perfusion spikes. The derivative gain is kept 3-4 orders of magnitude below K_p and K_i for the
same reason.

The simulator draws bubbles as a Poisson process (3 a day, depth 30-80% of the true value,
20-300 s long), filters the 10 s readings causally, and feeds the filtered glucose to the PID. The
paper quotes a 1 Hz filter sampling frequency; filtering the 10 s readings with the same
0.001 Hz cut-off gives the same time constant (about 160 s).

## Why there is no online cell count

Resuspending sedimented cells to count them disturbs their growth and risks the sterile barrier.
Bio-capacitance probes need suspended cells and space for a probe. So the cell concentration is
inferred, which is what the digital shadow does.

## Faults the simulator can reproduce

| Fault | Source | How |
| --- | --- | --- |
| Pump stop with data loss | run 3: 0.5 day, 2.5 days into perfusion | `PumpStop(start, duration)`; applied to run 3 by default |
| Bubble artefacts | section 2.2 | `BubbleModel` |
| Oscillating perfusion | run 6 (reduced chi-squared 3.2, attributed to PID tuning) | change `ControllerGains` |
