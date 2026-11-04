# Mathematical model

The model follows Egan et al., *Biotechnology Progress* 2025;41(6):e70045
([doi:10.1002/btpr.70045](https://doi.org/10.1002/btpr.70045)). Equation numbers refer to the main
text (1-18) and the supporting information (S1-S80). The full derivation and results are in
[technical_documentation/technical_documentation.pdf](technical_documentation/technical_documentation.pdf).

## Batch-perfusion process

Cells with concentration [X(t)] are seeded at time t_I and retained in the bioreactor via
sedimentation. During the sedimentation period (t_I <= t <= t_sed) there is no perfusion. During the
perfusion period (t_sed < t <= t_F), fresh media with glucose concentration [S_I] is perfused into the
bioreactor at the perfusion rate D(t) = Q(t)/V, and spent media with glucose and lactate
concentrations [S(t)] and [P(t)] is perfused out at the same rate, so that the volume V remains
constant. A PID controller adjusts D(t) to maintain [S(t)] at the set-point [S_sp]. The contents of
the bioreactor are modelled as homogeneous (an unstructured and unsegregated model).

## Equations

Cell growth: a lag period followed by a special case of the Richards growth model (eqs 10, 11):

```
[X(t)] = [X_I] / (1 - H(t - t_lag) mu_max nu (t - t_lag))^(1/nu)       for t >= t_I and nu < 0
mu(t)  = H(t - t_lag) mu_max / (1 - H(t - t_lag) mu_max nu (t - t_lag)) for t >= t_I and nu < 0
```

Glucose, lactate and the PI controlled perfusion rate (eqs 2-4), and the normalised dissolved
oxygen C'(t) = [C(t)]/[C*] of the extended model (eq. 8):

```
d[S]/dt = H(t - t_sed) D(t) ([S_I] - [S(t)]) - mu(t)[X(t)] / Y_xs
d[P]/dt = Y_px mu(t)[X(t)] - H(t - t_sed) D(t) [P(t)]
dD/dt   = H(t - t_sed) ( K_p d([S_sp] - [S(t)])/dt + K_i ([S_sp] - [S(t)]) )
dC'/dt  = (H(t - t_sed) D(t) + kLa)(1 - C'(t)) - mu(t)[X(t)] / Y'_xc
```

The ODEs are approximated by forward finite difference equations and solved via the forward Euler
method with a step size of 10 s, with [S], [P] >= 0 and D >= D_min enforced at each time step
(eqs S22-S24).

Derived quantities:

| Quantity | Expression | Equation |
| --- | --- | --- |
| Yield of lactate from glucose | Y_ps = Y_px Y_xs | 5 |
| Maximum specific growth rate from two cell concentrations | mu_max = (1 - ([X_I]/[X_F])^nu) / (nu (t_F - t_lag)) | 12 |
| Time to reach a threshold cell concentration [X_T] | t_T = t_lag + (1 - ([X_I]/[X_T])^nu) / (mu_max nu) | 13 |
| Lactate concentration steady state | [P_ss] = Y_ps ([S_I] - [S_sp]) | 18 |

## Observed and best-fit parameters

| Run | Cells, media | t_I | t_sed | t_F | [S_I] | [S_sp] | [X_I] | [X_F] | Y_xs | Y_px | t_lag - t_I | nu | mu_max | kLa | Y'_xc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | donor 1, A | 0.11 | 0.34 | 5.05 | 25 | 20 | 1.49 | 19.0 | 0.231 | 7.93 | 0.160 | -0.922 | 2.14 | 3.00 | 3.68 |
| 2 | donor 2, A | 0.14 | 1.02 | 6.03 | 25 | 20 | 0.37 | 19.0 | 0.275 | 6.25 | 0.172 | -0.477 | 2.02 | 741 | 0.0220 |
| 3 | donor 3, A | 0.11 | 0.86 | 6.01 | 25 | 20 | 1.49 | 32.5 | *0.259* | *7.62* | *0.092* | *-0.537* | *1.36* | *3.00* | *3.68* |
| 4 | patient 1, A | 2.09 | 2.84 | 7.00 | 25 | 20 | 0.62 | 10.9 | 0.272 | 9.09 | 0.128 | -0.524 | 1.39 | 1.73 | 6.27 |
| 5 | donor 4, A | 1.92 | 2.92 | 6.95 | 25 | 17 | 1.13 | 8.73 | 0.207 | 6.83 | 0.000 | -0.000 | 0.405 | 42.9 | 0.189 |
| 6 | donor 5, B | 0.08 | 0.83 | 4.85 | 26 | 21 | 1.49 | 18.2 | 0.311 | 7.99 | 0.000 | -0.761 | 1.58 | 20.2 | 0.893 |

Units: times in days from the start of the process; [S] in mmol/L; [X] in 10^6 cells/mL; Y_xs in
cells/pmol; Y_px in pmol/cell; mu_max and kLa in day^-1; Y'_xc in 10^6 cells/mL.

- t_I to [X_F]: Table 1. Y_xs to nu: Table 2. mu_max: eq. 12, as in Table 2. kLa and Y'_xc: Table
  S4; the 95% confidence intervals are very wide for runs 2, 5 and 6, indicating that these
  parameters are not practically identifiable.
- All runs: K_p = 0.215 L/(mmol day), K_i = 5.16 L/(mmol day^2), K_d = 7.46 x 10^-5 L/mmol,
  D_min = 0.067 day^-1, V = 134 mL.
- Media A: TexMACS with IL-7 and IL-15. Media B: Advanced TCM with IL-2. Run 4: patient CAR T cells
  (lentiviral transduction); all other runs: donor T cells.
- *Italics*: run 3 is the independent test dataset of the original study and has no best-fit
  parameters. It is assigned the mean training values, with mu_max from eq. 12 so that
  [X(t_F)] = [X_F], and the kLa and Y'_xc of run 1.

## Digital shadow training and testing

| | Known | Fitted | Fixed |
| --- | --- | --- | --- |
| Training | [X_I], [X_F], online D(t), [S(t)], [P(t)] | Y_xs, Y_px, t_lag, nu | mu_max via eq. 12 |
| Testing (real time) | [X_I], online data at the time of model fitting | Y_px, mu_max, nu | t_lag - t_I = 0.092 days, Y_xs = 0.259 cells/pmol (mean training values) |

In a real-time scenario either t_lag or mu_max must be known a priori for the model to be
structurally identifiable (Table S1); fixing t_lag alone leaves the model practically
unidentifiable, so Y_xs is also fixed. The extended model with the DO equation is not fitted, since
kLa and Y'_xc are not practically identifiable from these data.

## Assumptions

- The bioreactor contents are homogeneous (well mixed), although the cells are sedimented on the
  membrane and the influent and effluent are situated higher up.
- The glucose and lactate concentrations have a negligible inhibitory effect on cell growth; the
  dominant inhibitory effect is due to the cells competing for space.
- Cell death and the rates of glucose and lactate consumption and production associated with cell
  maintenance are negligible; the growth yield factors are constant.
- All cells are homogeneous (no CD4+/CD8+ segregation).
- Only a negligible number of cells are removed in the effluent.
