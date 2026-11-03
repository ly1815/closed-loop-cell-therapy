import numpy as np
import pytest

from cltherapy.process.parameters import PATIENT_A
from cltherapy.sensors import AtLineLab, ChannelSpec, Fault, FlowCytometry, OnlineSensors, SoftSensors
from cltherapy.sim.runner import Batch


def test_online_noise_and_drift():
    spec = ChannelSpec("pH", noise_sd=0.01, drift_per_day=-0.02, period_h=0.0, resolution=0.0001)
    s = OnlineSensors(channels=(spec,), seed=3)
    early = [s.read(h, {"pH": 7.2})["pH"] for h in np.linspace(0, 1, 400)]
    late = [s.read(h, {"pH": 7.2})["pH"] for h in np.linspace(240, 241, 400)]
    assert np.std(early) == pytest.approx(0.01, rel=0.2)
    assert np.mean(late) - np.mean(early) == pytest.approx(-0.2, abs=0.01)  # 10 days of drift


def test_reading_is_held_between_updates():
    spec = ChannelSpec("DO_pct", noise_sd=1.0, drift_per_day=0.0, period_h=1.0, resolution=0.1)
    s = OnlineSensors(channels=(spec,), seed=0)
    a = s.read(0.0, {"DO_pct": 50.0})
    b = s.read(0.5, {"DO_pct": 50.0})
    c = s.read(1.0, {"DO_pct": 50.0})
    assert a["DO_pct"] == b["DO_pct"] and b["DO_pct_age_h"] == 0.5
    assert c["DO_pct_age_h"] == 0.0


def test_faults_dropout_and_stale():
    s = OnlineSensors(seed=0)
    truth = {"pH": 7.2, "DO_pct": 50.0, "weight_kg": 1.0}
    s.inject("pH", Fault("dropout", start_h=1.0))
    assert s.read(0.5, truth)["pH"] is not None
    assert s.read(1.5, truth)["pH"] is None
    s.clear("pH")
    first = s.read(2.0, truth)["DO_pct"]
    s.inject("DO_pct", Fault("stale", start_h=2.0))
    assert s.read(5.0, {**truth, "DO_pct": 20.0})["DO_pct"] == first


def test_atline_results_arrive_after_turnaround():
    lab = AtLineLab(seed=0)
    lab.submit(10.0, {"viable_density": 5.0, "viability": 0.9, "glucose": 15.0, "lactate": 5.0})
    assert lab.poll(10.1) == []
    counter = lab.poll(10.3)
    assert [r.instrument for r in counter] == ["cell_counter"]
    metab = lab.poll(10.6)
    assert [r.instrument for r in metab] == ["metabolite_analyser"]
    assert metab[0].ready_h == pytest.approx(10.5)


def test_sampling_costs_cells():
    b = Batch(patient=PATIENT_A).run(5.0)
    before = b.truth()["total_viable_cells"]
    b.take_sample(2.0)
    after = b.truth()["total_viable_cells"]
    assert b.cells_sampled == pytest.approx(before * 2.0 / 1000.0, rel=1e-6)
    assert after == pytest.approx(before - b.cells_sampled, rel=1e-6)


def test_flow_cytometry_is_delayed():
    r = FlowCytometry(seed=0).measure(216.0, {"memory_fraction": 0.13})
    assert r.ready_h == pytest.approx(228.0)
    assert abs(r.values["memory_fraction"] - 0.13) < 0.1


def test_controller_follows_the_drifting_sensor():
    """A drifting pH patch drags the true pH away: the controller only sees the reading."""
    spec = ChannelSpec("pH", noise_sd=0.0, drift_per_day=-0.02, period_h=0.0, resolution=0.0001)
    b = Batch(patient=PATIENT_A, sensors=OnlineSensors(channels=(spec,) + OnlineSensors().channels[1:]))
    b.run(6.0)
    assert b.history[-1]["pH"] - 7.2 == pytest.approx(0.12, abs=0.02)


def test_soft_oxygen_uptake_tracks_truth():
    b = Batch(patient=PATIENT_A, sensors=OnlineSensors(seed=0)).run(8.0)
    h = b.history[1:]
    soft = SoftSensors().estimate(
        [r["t_h"] for r in h], [r["meas_DO_pct"] for r in h], [r["meas_pH"] for r in h],
        [r["y_O2"] for r in h], [r["y_CO2"] for r in h], [r["rock_rpm"] for r in h],
        [r["rock_angle_deg"] for r in h], [r["perfusion_vvd"] for r in h],
    )
    t = np.array([r["t_day"] for r in h])
    late = t > 6.0
    true_our = np.array([r["OUR"] for r in h])[late]
    est = soft["OUR_hat"][late]
    assert np.mean(np.abs(est - true_our) / true_our) < 0.12
