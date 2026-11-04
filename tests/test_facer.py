import numpy as np
import pytest

from cltherapy.control.perfusion import GlucosePerfusionPID
from cltherapy.process.facer import simulate_model
from cltherapy.process.growth import lactate_steady_state
from cltherapy.process.parameters import RUNS, TRAINING_RUNS
from cltherapy.sensors import CausalButterworth, OnlineSensors
from cltherapy.sim.runner import run_published


@pytest.mark.parametrize("n", TRAINING_RUNS)
def test_model_reaches_published_final_concentration_and_set_point(n):
    r = RUNS[n]
    m = simulate_model(r)
    end = m.at(r.t_F)
    assert end["X"] == pytest.approx(r.X_F, rel=1e-4)
    assert end["S"] == pytest.approx(r.S_sp, abs=0.25)
    # eq. 18: lactate settles near Y_ps (S_I - S_sp)
    assert end["P"] == pytest.approx(lactate_steady_state(r.Y_ps, r.S_I, r.S_sp), rel=0.05)


def test_no_perfusion_during_sedimentation_and_min_rate_after():
    r = RUNS[1]
    m = simulate_model(r)
    sed = (m.t >= r.t_I) & (m.t < r.t_sed)
    assert np.all(m.D[sed] == 0.0)
    assert np.all(m.D[m.t >= r.t_sed + 1e-4] >= r.gains.D_min - 1e-12)


def test_pid_raises_perfusion_when_glucose_is_low():
    pid = GlucosePerfusionPID(S_sp=20.0)
    assert pid.update(15.0, 0.01) == 0.0  # not started: sedimentation
    pid.start()
    D = [pid.update(s, 0.01) for s in (19.0, 18.0, 18.0, 18.0)]
    assert D[-1] > D[0] > 0
    assert pid.update(None, 0.01) == D[-1]  # no reading: hold


def test_butterworth_has_unit_gain_and_slows_steps():
    f = CausalButterworth()
    assert f(5.0) == pytest.approx(5.0)
    out = [f(10.0) for _ in range(30)]  # 5 minutes after a step at 10 s sampling
    assert 5.0 < out[-1] < 10.0
    out = [f(10.0) for _ in range(600)]  # well after
    assert out[-1] == pytest.approx(10.0, abs=1e-3)


def test_bubbles_pull_readings_down():
    s = OnlineSensors(seed=0, glucose_sd=0.0, lactate_sd=0.0)
    s.bubbles.rate_per_day = 86400 / 10.0  # a bubble on every reading
    r = s.read(20.0, 10.0, 0.8)
    assert r["S_raw"] < 20.0 * 0.81 and r["P_raw"] < 10.0 * 0.81


def test_virtual_facer_matches_model_without_noise():
    r = RUNS[2]
    v = run_published(r, realistic=False)
    m = simulate_model(r)
    assert v.X == pytest.approx(r.X_F, rel=1e-3)
    late = v.series("t") > r.t_F - 1
    assert np.mean(v.series("D")[late]) == pytest.approx(np.mean(m.D[m.t > r.t_F - 1]), rel=0.05)


def test_controller_holds_glucose_despite_bubbles():
    r = RUNS[1]
    v = run_published(r, seed=3)
    t, S = v.series("t"), v.series("S")
    late = t > r.t_sed + 2.0
    assert np.median(np.abs(S[late] - r.S_sp)) < 0.2


def test_run3_pump_stop_loses_data_and_perfusion():
    r = RUNS[3]
    v = run_published(r, seed=0)
    t = v.series("t")
    stop = (t > r.t_sed + 2.55) & (t < r.t_sed + 2.95)
    assert np.all(v.series("D")[stop] == 0.0)
    assert np.all(np.isnan(v.series("S_filt")[stop]))
