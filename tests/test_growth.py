import numpy as np
import pytest

from cltherapy.process.growth import (
    cell_concentration,
    lactate_steady_state,
    mu_max_from_final,
    sensitivities,
    specific_growth_rate,
    time_to_threshold,
)
from cltherapy.process.parameters import (
    X_T,
    PUBLISHED_EXPANSION_TO_1E7,
    PUBLISHED_MU_MAX,
    RUNS,
    TRAINING_RUNS,
)


def test_nu_minus_one_is_linear_and_minus_half_is_quadratic():
    """Special cases named in Egan et al. (S46-S50)."""
    t = np.linspace(0, 4, 9)
    X = cell_concentration(t, 2.0, 0.0, 1.5, -1.0)
    assert np.allclose(X, 2.0 * (1 + 1.5 * t))
    X = cell_concentration(t, 2.0, 0.0, 1.5, -0.5)
    assert np.allclose(X, 2.0 * (1 + 0.75 * t) ** 2)


def test_nu_to_zero_is_exponential():
    t = np.linspace(0, 3, 7)
    assert np.allclose(cell_concentration(t, 1.0, 0.0, 0.4, 0.0), np.exp(0.4 * t))
    assert np.allclose(cell_concentration(t, 1.0, 0.0, 0.4, -1e-5), np.exp(0.4 * t), rtol=1e-4)


def test_lag_holds_cells_constant():
    t = np.array([0.0, 0.1, 0.2])
    assert np.allclose(cell_concentration(t, 1.5, 0.2, 2.0, -0.9), 1.5)
    assert np.allclose(specific_growth_rate(t[:2], 0.2, 2.0, -0.9), 0.0)


def test_growth_rate_is_derivative_of_log_concentration():
    t = np.linspace(0.5, 4.0, 2000)
    X = cell_concentration(t, 1.2, 0.3, 1.8, -0.6)
    mu = specific_growth_rate(t, 0.3, 1.8, -0.6)
    assert np.allclose(np.gradient(np.log(X), t)[5:-5], mu[5:-5], rtol=1e-3)


@pytest.mark.parametrize("n", TRAINING_RUNS)
def test_mu_max_and_expansion_period_match_table_2(n):
    r = RUNS[n]
    assert r.mu_max == pytest.approx(PUBLISHED_MU_MAX[n], abs=0.015)
    tT = time_to_threshold(r.X_I, X_T, r.t_lag, r.mu_max, r.nu) - r.t_I
    assert tT == pytest.approx(PUBLISHED_EXPANSION_TO_1E7[n], abs=0.04)


def test_eq12_and_eq13_are_inverse():
    mu = mu_max_from_final(1.0, 20.0, 5.0, 0.2, -0.7)
    assert time_to_threshold(1.0, 20.0, 0.2, mu, -0.7) == pytest.approx(5.0)


def test_lactate_steady_state():
    assert lactate_steady_state(2.0, 25.0, 20.0) == pytest.approx(10.0)


def test_sensitivities_match_finite_differences():
    t = np.linspace(0.5, 5, 50)
    args = dict(X_I=1.3, t_lag=0.2)
    gX, gmu = sensitivities(t, mu_max=1.6, nu=-0.6, **args)
    h = 1e-6
    dX_dmu = (cell_concentration(t, mu_max=1.6 + h, nu=-0.6, **args)
              - cell_concentration(t, mu_max=1.6 - h, nu=-0.6, **args)) / (2 * h)
    dX_dnu = (cell_concentration(t, mu_max=1.6, nu=-0.6 + h, **args)
              - cell_concentration(t, mu_max=1.6, nu=-0.6 - h, **args)) / (2 * h)
    assert np.allclose(gX[:, 0], dX_dmu, rtol=1e-4)
    assert np.allclose(gX[:, 1], dX_dnu, rtol=1e-4)
    dmu_dnu = (specific_growth_rate(t, 0.2, 1.6, -0.6 + h) - specific_growth_rate(t, 0.2, 1.6, -0.6 - h)) / (2 * h)
    assert np.allclose(gmu[:, 1], dmu_dnu, rtol=1e-4)
