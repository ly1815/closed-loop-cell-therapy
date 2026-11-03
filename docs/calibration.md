# Calibration plan

The data are synthetic, but every parameter should trace to published CAR-T
manufacturing data. `process/parameters.py` marks each value as `anchor`,
`placeholder` or `assumed`. The current placeholders are plausible values
chosen so the model runs and meets the anchors; replace them as the
digitised data come in.

## Anchors the model already meets

| Quantity | Target | Source | Test |
| --- | --- | --- | --- |
| Stem-cell memory share of CAR+ cells, day 5 | about 25% | Cadinanos-Garai et al., Molecular Therapy 2025 | `test_memory_share_passes_literature_anchors` |
| Same, day 10 | about 11.5% | same | same |
| Slow donor's cell count versus a fit donor | 30 to 50% lower | Jagannathan et al. (SMART), Cell Reports Methods 2026, at day 12 | `test_slow_donor_has_30_to_50_percent_fewer_cells` (day 9 here) |

`k_diff0`, `diff_time_scale_h`, `mu_max_M` and `mu_max_E` were fitted to the
two memory anchors and a day-9 density of 20e6 cells/mL for patient A, so
the memory anchors are met by construction. Treat them as a calibration
check, not a validation.

## To digitise next

| Parameter | Paper and figure to use | Notes |
| --- | --- | --- |
| `mu_max_M`, `mu_max_E`, `lag_tau_h` | Costariol 2019 growth curves (ambr 250); Egan 2025 cell counts | fit the total viable density first |
| `q_glc_*`, `Y_lac_glc` | Costariol 2019 specific rates; SMART 2026 glucose and lactate rates | Costariol reports lactate yield from glucose |
| `q_O2_*` | SMART 2026 oxygen consumption estimates | memory cells rely more on oxidative metabolism |
| `activation`, `early_lac_boost` | SMART 2026 early-lactate correlations | keep the sign: early lactate goes with more central memory |
| `seed_memory_fraction` | Cadinanos-Garai 2025, day 0 phenotype | check the share at day 0 |
| `acid_to_bicarb`, `hco3_in` | Egan 2025 pH and perfusion traces | lumped buffer capacity |
| `kLa_ref` | rocking-bag kLa literature, or measure per bag | depends on rock rate and angle |

Tools: WebPlotDigitizer for figure panels; save each digitised series as a
CSV under `data/literature/<paper>/<figure>.csv` with the figure reference
in the file name.
