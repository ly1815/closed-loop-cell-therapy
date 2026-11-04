# Online and off-line data

Measurements of the virtual bioreactor, following section 2.2 of Egan et al. (2025). The noise
levels and bubble statistics are assumed values chosen to resemble the published online data; the
measurement frequencies and the filter specification are those of the original study.

| Variable | Measurement | Frequency | Virtual bioreactor |
| --- | --- | --- | --- |
| Glucose concentration [S(t)] | single-use enzymatic sensor (online) | 10 s | Gaussian noise, sd 0.08 mmol/L; bubble artefacts; Butterworth filter |
| Lactate concentration [P(t)] | single-use enzymatic sensor (online) | 10 s | as for glucose |
| Normalised DO, [C(t)]/[C*] | single-use sensor (online) | online | Gaussian noise, sd 0.01 |
| pH | single-use sensor (online) | online | not modelled (not part of the model) |
| Perfusion rate D(t) | PID controller output (online) | 10 s | exact |
| Cell concentration [X] | cell counter, suspended cells (off-line) | t_I and t_F only | coefficient of variation 5% |
| Phenotype | flow cytometry, suspended or sedimented cells (off-line) | | not modelled |

## Bubble artefacts and filtering

Bubbles generated during the expansion process interfere with the online sensors, leading to
artificial drops in the glucose and lactate concentration readings, followed by rises after the
bubbles have moved away from the sensors. To prevent these fluctuations impacting the perfusion
rate, the raw glucose and lactate data are filtered before being passed to the PID controller by a
lowpass Butterworth filter of order 2 with a cut-off frequency of 0.001 Hz. Not all bubble
artefacts are filtered out, which results in downward spikes of the glucose (and lactate)
concentration followed by upward spikes in the perfusion rate. The derivative gain K_d is set
3-4 orders of magnitude lower than K_p and K_i to mitigate the impact of the bubbles on the
perfusion rate.

In the virtual bioreactor, bubbles occur as a Poisson process (3 per day); during a bubble the
readings fall to a fraction of 0.3-0.8 of the true concentration for 20-300 s. The filter is applied
causally to the 10 s readings with the published cut-off frequency. (The original filter is
specified with a sampling frequency of 1 Hz; the cut-off frequency, and hence the time constant of
approximately 160 s, is the same.)

## Off-line cell concentration

Accurately measuring the cell concentration requires the cells to be in suspension, and
resuspending the cells after they have sedimented can disturb their growth. Suspension samples are
therefore taken only prior to sedimentation and when perfusion has ended, giving [X_I] and [X_F].
Hard sensors such as dielectric spectroscopy are less suited to bioreactors in which the cells grow
whilst sedimented on a membrane, which motivates the soft sensor of cell concentration.

## Process disturbances

| Disturbance | Origin in the original study | Virtual bioreactor |
| --- | --- | --- |
| Mechanical issue: perfusion stopped, no online data | run 3, 2.5 days after the start of perfusion, approximately 0.5 days | `PumpStop(start, duration)`, applied to run 3 by default |
| Bubble artefacts | all runs (section 2.2) | `BubbleModel` |
| Oscillations in the online data | run 6 (reduced chi-squared 3.2), possibly induced by the choice of PID parameters | modify `ControllerGains` |
