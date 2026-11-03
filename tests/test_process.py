import numpy as np
import pytest

from cltherapy.process.parameters import MEMORY_ANCHORS, PATIENT_A, PATIENT_B
from cltherapy.sim.runner import run_batch


@pytest.fixture(scope="module")
def runs():
    return {name: run_batch(p, days=10.0) for name, p in (("A", PATIENT_A), ("B", PATIENT_B))}


def test_states_stay_physical(runs):
    for b in runs.values():
        for key in ("x_M", "x_E", "x_G", "x_L", "x_O"):
            assert np.all(b.series(key) >= 0.0)
        assert np.all(b.series("viability") <= 1.0)


def test_memory_share_passes_literature_anchors(runs):
    """Patient A reproduces the day 5 and day 10 stem-cell memory shares (Cadinanos-Garai 2025)."""
    b = runs["A"]
    for day, target in MEMORY_ANCHORS.items():
        assert b.at_day(day)["memory_fraction"] == pytest.approx(target, abs=0.03)


def test_dose_rises_while_memory_falls(runs):
    b = runs["A"]
    days = (3, 5, 7, 9)
    dose = [b.at_day(d)["total_viable_cells"] for d in days]
    memory = [b.at_day(d)["memory_fraction"] for d in days]
    assert all(a < c for a, c in zip(dose, dose[1:]))
    assert all(a > c for a, c in zip(memory, memory[1:]))


def test_slow_donor_has_30_to_50_percent_fewer_cells(runs):
    """SMART (2026) reported 30-50% lower day-12 counts for patient cells; day 9 is used here."""
    ratio = runs["B"].at_day(9)["viable_density"] / runs["A"].at_day(9)["viable_density"]
    assert 0.5 <= ratio <= 0.7


def test_inner_loop_holds_pH_and_DO(runs):
    b = runs["A"]
    t = b.series("t_day")
    do = b.series("DO_pct")
    # early on, few cells use little oxygen: DO sits above setpoint at the minimum rock rate
    reached = np.argmax((t > 2.0) & (do <= 51.0))
    late = (np.arange(len(t)) >= reached) & (t < 9.5)
    assert late.sum() > 100
    assert np.max(np.abs(b.series("pH")[late] - 7.2)) < 0.05
    assert np.max(np.abs(do[late] - 50.0)) < 5.0
