"""Two-population T-cell expansion in a perfused rocking bag.

States (index, name, unit):
    0 M    memory-like viable cells     1e6 cells/mL
    1 E    effector-like viable cells   1e6 cells/mL
    2 Xd   dead cells                   1e6 cells/mL
    3 G    glucose                      mM
    4 L    lactate                      mM
    5 B    bicarbonate                  mM
    6 O    dissolved oxygen             mM

Inputs (held constant over one control interval):
    perfusion_vvd   vessel volumes per day of fresh medium (cells retained)
    rock_rpm        rock rate
    rock_angle_deg  rock angle
    y_O2            O2 fraction in the headspace gas
    y_CO2           CO2 fraction in the headspace gas

Memory-like cells differentiate into effector-like cells at a rate that rises
with time in culture. A per-patient activation parameter raises early lactate
production and slows later differentiation, following the correlation that
early lactate goes with more central memory at harvest (Jagannathan et al.,
Cell Reports Methods 2026). The structure is a hypothesis to test.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from .parameters import PatientParams, ProcessParams

STATE_NAMES = ("M", "E", "Xd", "G", "L", "B", "O")
N_STATES = len(STATE_NAMES)


@dataclass
class Inputs:
    perfusion_vvd: float = 0.0
    rock_rpm: float = 15.0
    rock_angle_deg: float = 6.0
    y_O2: float = 0.21
    y_CO2: float = 0.05


def kla(proc: ProcessParams, rock_rpm: float, angle_deg: float) -> float:
    """Oxygen transfer coefficient (1/h) as a power law in rock rate and angle."""
    rpm = max(rock_rpm, 0.0)
    ang = max(angle_deg, 0.0)
    return (
        proc.kLa_ref
        * (rpm / proc.rock_ref_rpm) ** proc.kLa_rpm_exp
        * (ang / proc.angle_ref_deg) ** proc.kLa_angle_exp
    )


def ph_from(proc: ProcessParams, bicarbonate: float, y_CO2: float) -> float:
    """Henderson-Hasselbalch for the bicarbonate buffer at 37 C."""
    pco2 = max(y_CO2, 1e-4) * proc.pCO2_per_frac
    return proc.pKa_bicarb + np.log10(max(bicarbonate, 1e-3) / (proc.CO2_solubility * pco2))


def initial_state(proc: ProcessParams, pat: PatientParams, y_O2: float = 0.21) -> np.ndarray:
    viable = pat.seed_density
    dead = viable * (1.0 - pat.seed_viability) / pat.seed_viability
    return np.array(
        [
            viable * pat.seed_memory_fraction,
            viable * (1.0 - pat.seed_memory_fraction),
            dead,
            proc.glc_in,
            proc.lac_in,
            proc.hco3_in,
            proc.O2_sat_air * y_O2 / 0.21,
        ]
    )


def rates(t_h: float, x: np.ndarray, u: Inputs, proc: ProcessParams, pat: PatientParams) -> dict:
    """Specific and volumetric rates at one instant. Shared by the ODE and the soft sensors."""
    M, E, Xd, G, L, B, O = np.maximum(x, 0.0)

    lag = 1.0 - np.exp(-t_h / proc.lag_tau_h)
    f_glc = G / (proc.K_glc + G)
    f_lac = proc.K_i_lac / (proc.K_i_lac + L)
    pH = ph_from(proc, B, u.y_CO2)
    f_pH = 1.0 if pH >= proc.pH_opt_low else max(0.0, 1.0 - 4.0 * (proc.pH_opt_low - pH))
    f_O2 = O / (0.005 + O)
    f = lag * f_glc * f_lac * f_pH * f_O2

    mu_M = pat.mu_max_M * f
    mu_E = pat.mu_max_E * f
    k_diff = (
        pat.k_diff0
        * 2.0 ** (t_h / proc.diff_time_scale_h)
        * (1.0 - proc.activation_diff_reduction * pat.activation)
    )

    burst = 1.0 + proc.early_lac_boost * pat.activation * np.exp(-t_h / proc.early_lac_tau_h)
    glc_uptake = (proc.q_glc_M * M + proc.q_glc_E * E) * f_glc  # mM/h
    lac_prod = proc.Y_lac_glc * glc_uptake * burst  # mM/h
    our = (proc.q_O2_M * M + proc.q_O2_E * E) * f_O2  # mM/h

    return {
        "mu_M": mu_M,
        "mu_E": mu_E,
        "k_diff": k_diff,
        "glc_uptake": glc_uptake,
        "lac_prod": lac_prod,
        "our": our,
        "pH": pH,
        "kLa": kla(proc, u.rock_rpm, u.rock_angle_deg),
        "O2_star": proc.O2_sat_air * u.y_O2 / 0.21,
    }


def rhs(t_h: float, x: np.ndarray, u: Inputs, proc: ProcessParams, pat: PatientParams) -> np.ndarray:
    M, E, Xd, G, L, B, O = np.maximum(x, 0.0)
    r = rates(t_h, x, u, proc, pat)
    D = u.perfusion_vvd / 24.0  # 1/h; cells are retained by the perfusion filter

    dM = (r["mu_M"] - r["k_diff"] - proc.k_death_M) * M
    dE = (r["mu_E"] - proc.k_death_E) * E + r["k_diff"] * M
    dXd = proc.k_death_M * M + proc.k_death_E * E
    dG = D * (proc.glc_in - G) - r["glc_uptake"]
    dL = D * (proc.lac_in - L) + r["lac_prod"]
    dB = D * (proc.hco3_in - B) - proc.acid_to_bicarb * r["lac_prod"]  # lactic acid consumes bicarbonate
    dO = r["kLa"] * (r["O2_star"] - O) - r["our"]
    return np.array([dM, dE, dXd, dG, dL, dB, dO])


def observables(t_h: float, x: np.ndarray, u: Inputs, proc: ProcessParams, pat: PatientParams) -> dict:
    """Quantities an operator or a sensor would see, plus the hidden truth."""
    M, E, Xd, G, L, B, O = np.maximum(x, 0.0)
    r = rates(t_h, x, u, proc, pat)
    viable = M + E
    return {
        "viable_density": viable,
        "viability": viable / max(viable + Xd, 1e-9),
        "memory_fraction": M / max(viable, 1e-9),  # hidden: only flow cytometry sees it
        "total_viable_cells": viable * proc.volume_L * 1e3,  # 1e6 cells
        "glucose": G,
        "lactate": L,
        "pH": r["pH"],
        "DO_pct": 100.0 * O / proc.O2_sat_air,  # % of air saturation
        "OUR": r["our"],
        "acid_production": r["lac_prod"],
    }
