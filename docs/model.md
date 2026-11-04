# The FACER model

Everything here follows Egan et al., Biotechnology Progress 2025 (doi:10.1002/btpr.70045),
main text and supporting information. Units: time in days, cells in 1e6 cells/mL, glucose S and
lactate P in mmol/L, Y_xs in cells/pmol, Y_px in pmol/cell; then mu X / Y_xs is in mmol/L/day.

## Process

1. Cells are seeded at t_I into the 134 mL chamber of the FACER cartridge (runs 4 and 5 first
   spend 2 days in the 13.4 mL chamber, not modelled).
2. **Sedimentation**, t_I to t_sed: no perfusion; the cells settle onto the gas-permeable membrane.
3. **Perfusion**, t_sed to t_F: fresh medium in and spent medium out at the same rate D, set by a
   PID that holds glucose at S_sp. Cells stay sedimented, so none leave in the effluent.
4. Cells are resuspended and counted only at seeding and at the end.

## Equations

Growth, the special case of the Richards model with a lag (eqs 10, 11; S60, S61), with nu < 0:

```
X(t)  = X_I / (1 - H(t - t_lag) mu_max nu (t - t_lag))^(1/nu)
mu(t) = H(t - t_lag) mu_max / (1 - H(t - t_lag) mu_max nu (t - t_lag))
```

nu = -1 gives linear growth, nu = -1/2 quadratic, nu -> 0 exponential. Balances (eqs 2-4, 8):

```
dS/dt  = H(t - t_sed) D (S_I - S) - mu X / Y_xs
dP/dt  = Y_px mu X - H(t - t_sed) D P
dD/dt  = H(t - t_sed) [ -K_p dS/dt + K_i (S_sp - S) ],   D >= D_min
dC'/dt = (H(t - t_sed) D + kLa)(1 - C') - mu X / Y'_xc      (normalised DO, extended model)
```

Solved by forward Euler with a 10 s step, clipping every state at its bound (S22-S24).
Useful closed forms: mu_max from the two counts (eq. 12), time to a threshold X_T (eq. 13), and
the lactate steady state P_ss = Y_ps (S_I - S_sp) with Y_ps = Y_px Y_xs (eqs 5, 18).

## Parameters

| Run | Cells, media | t_I | t_sed | t_F | S_I | S_sp | X_I | X_F | Y_xs | Y_px | t_lag - t_I | nu | mu_max | kLa | Y'_xc |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | donor 1, A | 0.11 | 0.34 | 5.05 | 25 | 20 | 1.49 | 19.0 | 0.231 | 7.93 | 0.160 | -0.922 | 2.14 | 3.00 | 3.68 |
| 2 | donor 2, A | 0.14 | 1.02 | 6.03 | 25 | 20 | 0.37 | 19.0 | 0.275 | 6.25 | 0.172 | -0.477 | 2.02 | 741 | 0.022 |
| 3 | donor 3, A | 0.11 | 0.86 | 6.01 | 25 | 20 | 1.49 | 32.5 | *0.259* | *7.62* | *0.092* | *-0.537* | *1.36* | *3.00* | *3.68* |
| 4 | patient 1, A | 2.09 | 2.84 | 7.00 | 25 | 20 | 0.62 | 10.9 | 0.272 | 9.09 | 0.128 | -0.524 | 1.39 | 1.73 | 6.27 |
| 5 | donor 4, A | 1.92 | 2.92 | 6.95 | 25 | 17 | 1.13 | 8.73 | 0.207 | 6.83 | 0.000 | -0.000 | 0.405 | 42.9 | 0.189 |
| 6 | donor 5, B | 0.08 | 0.83 | 4.85 | 26 | 21 | 1.49 | 18.2 | 0.311 | 7.99 | 0.000 | -0.761 | 1.58 | 20.2 | 0.893 |

- Columns t_I to X_F: Table 1. Y_xs to nu: Table 2. mu_max: eq. 12 (as in Table 2). kLa, Y'_xc:
  Table S4, whose intervals are very wide for runs 2, 5 and 6 (not practically identifiable).
- All runs: K_p = 0.215 L/(mmol day), K_i = 5.16 L/(mmol day^2), K_d = 7.46e-5, D_min = 0.067/day.
- Media A: TexMACS with IL-7 and IL-15; media B: Advanced TCM with IL-2. Run 4 is patient CAR-T
  cells (lentiviral transduction); the others are donor T cells.
- *Italics*: run 3 was the paper's independent test run and was never fitted. Its growth and yield
  parameters here are the training means, with mu_max from eq. 12 so it ends at the published
  X_F; kLa and Y'_xc are borrowed from run 1. Treat it as derived.

## What the digital shadow fits

| Mode | Known | Fitted | Fixed |
| --- | --- | --- | --- |
| Training (complete run) | X_I, X_F, online D, S, P | Y_xs, Y_px, t_lag, nu | mu_max from eq. 12 |
| Real time (data so far) | X_I, online D, S, P | Y_px, mu_max, nu | t_lag - t_I = 0.092, Y_xs = 0.259 (training means) |

The real-time model is structurally identifiable only with mu_max or t_lag known (Table S1), and
fixing t_lag alone was not practically identifiable, so Y_xs is fixed too. The DO equation is left
out of the shadow because kLa and Y'_xc are not practically identifiable from these data.

## Assumptions, as listed by the paper

- The chamber is treated as well mixed, although cells are sedimented and medium enters and leaves
  higher up. This is probably the strongest assumption.
- Glucose stays high and lactate low enough not to inhibit growth; growth slows only because cells
  compete for space.
- No cell death and no maintenance consumption; constant yields.
- All cells are alike (no CD4/CD8 or phenotype segregation).
