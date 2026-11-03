"""Process and patient parameters for the two-population CAR-T expansion model.

Every value carries a unit, a status and a source in PARAM_INFO:

- "anchor":      a published number the model is calibrated to reproduce
- "placeholder": a plausible value chosen for the demo, to be replaced by a
                 fit to digitised literature data (see docs/calibration.md)
- "assumed":     a modelling choice, not a measured quantity

Time is in hours, cell densities in 1e6 cells/mL, concentrations in mM,
and specific rates in pmol/cell/h (so q * X gives mM/h directly).
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace


@dataclass(frozen=True)
class PatientParams:
    """Parameters that differ from patient to patient (the EnKF learns these)."""

    name: str = "patient"
    mu_max_M: float = 0.0205  # 1/h, memory-like cells
    mu_max_E: float = 0.0306  # 1/h, effector-like cells
    k_diff0: float = 0.0018  # 1/h, memory -> effector at day 0
    activation: float = 0.5  # 0..1, strength of early activation
    seed_density: float = 0.5  # 1e6 cells/mL at day 0
    seed_memory_fraction: float = 0.45  # share of memory-like cells at day 0
    seed_viability: float = 0.95


@dataclass(frozen=True)
class ProcessParams:
    """Parameters of the cells, the medium and the rocking bag."""

    # growth and death
    lag_tau_h: float = 24.0  # activation lag before full growth
    diff_time_scale_h: float = 210.0  # k_diff doubles every this many hours
    activation_diff_reduction: float = 0.5  # strong activation slows k_diff
    k_death_M: float = 0.0010  # 1/h
    k_death_E: float = 0.0030  # 1/h
    K_glc: float = 0.5  # mM, Monod constant for glucose
    K_i_lac: float = 40.0  # mM, lactate inhibition constant
    pH_opt_low: float = 6.95  # growth falls below this pH

    # metabolism
    q_glc_M: float = 0.035  # pmol/cell/h
    q_glc_E: float = 0.060  # pmol/cell/h, effectors are more glycolytic
    Y_lac_glc: float = 1.6  # mol lactate per mol glucose
    early_lac_boost: float = 1.0  # extra lactate early, scaled by activation
    early_lac_tau_h: float = 36.0  # decay time of the early glycolytic burst
    q_O2_M: float = 0.030  # pmol/cell/h, memory cells rely more on OXPHOS
    q_O2_E: float = 0.020  # pmol/cell/h

    # medium
    glc_in: float = 20.0  # mM in fresh medium
    hco3_in: float = 24.0  # mM bicarbonate in fresh medium
    acid_to_bicarb: float = 0.4  # share of lactic acid taken up by bicarbonate; the rest by HEPES and other buffers
    lac_in: float = 0.0

    # rocking bag
    volume_L: float = 1.0
    O2_sat_air: float = 0.19  # mM dissolved O2 at 37 C under air
    kLa_ref: float = 2.0  # 1/h at the reference rock rate and angle
    rock_ref_rpm: float = 15.0
    angle_ref_deg: float = 6.0
    kLa_rpm_exp: float = 0.8
    kLa_angle_exp: float = 0.5
    pCO2_per_frac: float = 713.0  # mmHg of CO2 per unit gas fraction at 37 C
    pKa_bicarb: float = 6.1
    CO2_solubility: float = 0.0307  # mM per mmHg


PARAM_INFO: dict[str, tuple[str, str, str]] = {
    # name: (unit, status, source or note)
    "mu_max_M": ("1/h", "placeholder", "fit to Costariol 2019 and Egan 2025 growth curves"),
    "mu_max_E": ("1/h", "placeholder", "fit to Costariol 2019 and Egan 2025 growth curves"),
    "k_diff0": ("1/h", "placeholder", "fitted with diff_time_scale_h so patient A matches the day 5 and day 10 anchors"),
    "diff_time_scale_h": ("h", "placeholder", "fitted with k_diff0 to the memory anchors"),
    "activation": ("-", "assumed", "per-patient; links early lactate to later memory (SMART 2026)"),
    "seed_memory_fraction": ("-", "placeholder", "starting phenotype, check against Cadinanos-Garai 2025 day 0"),
    "q_glc_M": ("pmol/cell/h", "placeholder", "fit to SMART 2026 and Costariol 2019 rates"),
    "q_glc_E": ("pmol/cell/h", "placeholder", "fit to SMART 2026 and Costariol 2019 rates"),
    "Y_lac_glc": ("mol/mol", "placeholder", "Costariol 2019 lactate yield from glucose"),
    "q_O2_M": ("pmol/cell/h", "placeholder", "SMART 2026 oxygen consumption estimates"),
    "q_O2_E": ("pmol/cell/h", "placeholder", "SMART 2026 oxygen consumption estimates"),
    "kLa_ref": ("1/h", "placeholder", "rocking-bag kLa, to measure per bag and setting"),
    "glc_in": ("mM", "placeholder", "medium composition, check the medium used in the calibration papers"),
    "hco3_in": ("mM", "assumed", "typical bicarbonate-buffered medium"),
    "acid_to_bicarb": ("-", "assumed", "lumped buffer capacity; check against the pH and CO2 traces in Egan 2025"),
}

# Published phenotype anchors the model is calibrated to (share of CAR+ cells).
# Cadinanos-Garai et al., Molecular Therapy 2025, 6 healthy donors.
MEMORY_ANCHORS: dict[float, float] = {5.0: 0.25, 10.0: 0.115}

DEFAULT_PROCESS = ProcessParams()

PATIENT_A = PatientParams(name="Patient A (fit donor)")

# Slow donor: growth scaled so day-9 density is 40% lower than patient A,
# in the 30-50% range SMART (2026) saw for patient cells at day 12.
PATIENT_B = replace(
    PATIENT_A,
    name="Patient B (slow donor)",
    mu_max_M=0.0172,
    mu_max_E=0.0256,
    activation=0.25,
    seed_viability=0.85,
)

PATIENTS = {"A": PATIENT_A, "B": PATIENT_B}
