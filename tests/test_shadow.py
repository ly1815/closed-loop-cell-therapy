import numpy as np
import pytest

from cltherapy.process.parameters import RUNS
from cltherapy.shadow import OnlineData, fit_realtime, fit_training
from cltherapy.sim.runner import run_published


@pytest.fixture(scope="module")
def run2():
    r = RUNS[2]
    v = run_published(r, seed=2)
    return r, v, OnlineData.from_dict(v.online_data())


def test_training_fit_recovers_published_parameters(run2):
    r, v, data = run2
    f = fit_training(r, data, r.X_I, r.X_F)  # exact counts, so only sensor noise remains
    assert f.params["Y_xs"] == pytest.approx(r.Y_xs, rel=0.05)
    assert f.params["Y_px"] == pytest.approx(r.Y_px, rel=0.08)
    assert f.params["nu"] == pytest.approx(r.nu, abs=0.05)
    assert f.params["t_lag"] - r.t_I == pytest.approx(r.lag, abs=0.03)


def test_realtime_forecast_with_2_5_days_is_within_30_percent(run2):
    """Egan et al.: with 2.5 days of online data every run was predicted within 30%."""
    r, v, data = run2
    f = fit_realtime(r, data, v.X_I_count, r.t_sed + 2.5)
    assert float(f.predict_X(r.t_F)) == pytest.approx(r.X_F, rel=0.30)
    assert f.t_now <= r.t_sed + 2.5 + 1e-9


def test_bands_contain_the_estimate_and_narrow_with_more_data(run2):
    r, v, data = run2
    t = np.array([r.t_F])
    widths = []
    for w in (1.5, 3.0):
        f = fit_realtime(r, data, v.X_I_count, r.t_sed + w)
        b = f.band(t)
        assert b["X_lo"][0] <= b["X"][0] <= b["X_hi"][0]
        widths.append(b["X_hi"][0] - b["X_lo"][0])
    assert widths[1] < widths[0]


def test_threshold_time_is_where_forecast_crosses_threshold(run2):
    r, v, data = run2
    f = fit_realtime(r, data, v.X_I_count, r.t_sed + 2.5)
    tT = f.threshold_time(10.0)
    assert float(f.predict_X(tT)) == pytest.approx(10.0, rel=1e-6)
