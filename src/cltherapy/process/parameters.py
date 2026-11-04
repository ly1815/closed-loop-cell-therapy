"""The six FACER runs of Egan et al. (2025), as published.

Sources (Egan et al., Biotechnology Progress 2025, doi:10.1002/btpr.70045):
- Table 1: observed inputs (times, controller gains, media, set points, cell concentrations)
- Table 2: best-fit growth and yield parameters for the five training runs
- Table S4: kLa and the scaled oxygen yield from the extended model fitted to DO data

Units: time in days from the start of data collection, cell concentration in 1e6 cells/mL,
glucose and lactate in mmol/L, Y_xs in cells/pmol, Y_px in pmol/cell, kLa in 1/day,
Y'_xc in 1e6 cells/mL. With these units mu*X/Y_xs and Y_px*mu*X come out in mmol/L/day.

Run 3 was never fitted in the paper (a mechanical stop 2.5 days into perfusion). Its growth
parameters here are the training means, with mu_max set by eq. 12 so the run ends at the
published final concentration: a derived run, flagged as such.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace

from .growth import mu_max_from_final


@dataclass(frozen=True)
class FacerChamber:
    name: str
    volume_mL: float
    area_cm2: float


SMALL_CHAMBER = FacerChamber("small", volume_mL=13.4, area_cm2=10.0)
LARGE_CHAMBER = FacerChamber("large", volume_mL=134.0, area_cm2=100.0)


@dataclass(frozen=True)
class ControllerGains:
    """PID on glucose, acting on the perfusion rate (Table 1, identical for all runs)."""

    kp: float = 0.215  # L/(mmol day)
    ki: float = 5.16  # L/(mmol day^2)
    kd: float = 7.46e-5  # L/mmol
    D_min: float = 0.067  # 1/day


@dataclass(frozen=True)
class FacerRun:
    run: int
    cells: str  # donor or patient
    media: str  # "A": TexMACS + IL-7 + IL-15; "B": Advanced TCM + IL-2
    # Table 1
    t_I: float  # seeding time in the large chamber
    t_sed: float  # end of sedimentation, start of perfusion
    t_F: float  # end of perfusion
    S_I: float  # glucose in fresh medium
    S_sp: float  # glucose set point
    X_I: float  # initial cell concentration
    X_F: float  # final cell concentration (offline count)
    # Table 2 (training fits)
    Y_xs: float
    Y_px: float
    lag: float  # t_lag - t_I
    nu: float
    # Table S4 (extended fit to DO data); wide intervals for runs 2, 5 and 6
    kLa: float
    Y_xc: float
    gains: ControllerGains = field(default_factory=ControllerGains)
    chamber: FacerChamber = LARGE_CHAMBER
    derived: bool = False  # True for run 3, whose parameters are not from Table 2
    note: str = ""

    @property
    def t_lag(self) -> float:
        return self.t_I + self.lag

    @property
    def mu_max(self) -> float:
        """Not fitted directly in training: eq. 12 from X_I, X_F and the fitted t_lag, nu."""
        return mu_max_from_final(self.X_I, self.X_F, self.t_F, self.t_lag, self.nu)

    @property
    def Y_ps(self) -> float:
        """Yield of lactate from glucose, eq. 5."""
        return self.Y_px * self.Y_xs


RUNS: dict[int, FacerRun] = {
    1: FacerRun(1, "donor 1", "A", 0.11, 0.34, 5.05, 25, 20, 1.49, 19.0,
                Y_xs=0.231, Y_px=7.93, lag=0.160, nu=-0.922, kLa=3.00, Y_xc=3.68),
    2: FacerRun(2, "donor 2", "A", 0.14, 1.02, 6.03, 25, 20, 0.37, 19.0,
                Y_xs=0.275, Y_px=6.25, lag=0.172, nu=-0.477, kLa=741.0, Y_xc=0.0220),
    4: FacerRun(4, "patient 1", "A", 2.09, 2.84, 7.00, 25, 20, 0.62, 10.9,
                Y_xs=0.272, Y_px=9.09, lag=0.128, nu=-0.524, kLa=1.73, Y_xc=6.27,
                note="patient CAR-T cells, lentiviral transduction; 2 days in the small chamber first"),
    5: FacerRun(5, "donor 4", "A", 1.92, 2.92, 6.95, 25, 17, 1.13, 8.73,
                Y_xs=0.207, Y_px=6.83, lag=0.000, nu=-0.000, kLa=42.9, Y_xc=0.189,
                note="2 days in the small chamber first; nu fitted as zero, so growth is exponential"),
    6: FacerRun(6, "donor 5", "B", 0.08, 0.83, 4.85, 26, 21, 1.49, 18.2,
                Y_xs=0.311, Y_px=7.99, lag=0.000, nu=-0.761, kLa=20.2, Y_xc=0.893,
                note="media B; oscillating perfusion, reduced chi-squared 3.2"),
}

TRAINING_RUNS = (1, 2, 4, 5, 6)

# Means over the training runs (Table 2), used to fix t_lag and Y_xs in the real-time shadow.
MEAN_LAG = 0.092
MEAN_Y_XS = 0.259
MEAN_Y_PX = 7.62
MEAN_NU = -0.537

RUNS[3] = FacerRun(3, "donor 3", "A", 0.11, 0.86, 6.01, 25, 20, 1.49, 32.5,
                   Y_xs=MEAN_Y_XS, Y_px=MEAN_Y_PX, lag=MEAN_LAG, nu=MEAN_NU,
                   kLa=3.00, Y_xc=3.68, derived=True,
                   note="independent test run; perfusion stopped for 0.5 day, 2.5 days into perfusion")
RUNS = dict(sorted(RUNS.items()))

# Published values the simulator should reproduce (Table 2), for checks.
PUBLISHED_MU_MAX = {1: 2.14, 2: 2.02, 4: 1.39, 5: 0.405, 6: 1.58}
PUBLISHED_EXPANSION_TO_1E7 = {1: 2.59, 2: 4.14, 4: 4.65, 5: 5.38, 6: 2.71}  # t_T - t_I, days

DOSE_THRESHOLD = 10.0  # 1e6 cells/mL, the threshold used in Table 2


def with_changes(run: FacerRun, **kwargs) -> FacerRun:
    return replace(run, **kwargs)
