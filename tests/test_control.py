from cltherapy.control.pid import PI
from cltherapy.control.rocker import PerfusionSchedule, RockerController


def test_pi_holds_first_order_setpoint():
    """A PI loop drives a simple first-order process to its setpoint."""
    pi = PI(kp=2.0, ti_h=0.5, setpoint=1.0, out_min=0.0, out_max=5.0)
    y, dt = 0.0, 0.01
    for _ in range(2000):
        u = pi.update(y, dt)
        y += dt * (u - y) / 0.5
    assert abs(y - 1.0) < 1e-3


def test_pi_anti_windup_recovers_quickly():
    pi = PI(kp=1.0, ti_h=0.2, setpoint=10.0, out_min=0.0, out_max=1.0)
    for _ in range(500):  # unreachable setpoint: output saturates
        pi.update(0.0, 0.01)
    pi.setpoint = 0.0
    out = pi.update(0.5, 0.01)
    assert out == 0.0  # no stored integral holding the output high


def test_perfusion_schedule_ramps():
    s = PerfusionSchedule(start_day=3.0, ramp_days=2.0, max_vvd=1.5)
    assert s.rate(2.9) == 0.0
    assert abs(s.rate(4.0) - 0.75) < 1e-9
    assert s.rate(8.0) == 1.5


def test_supervisory_setpoints_are_clipped():
    c = RockerController()
    c.set_setpoints(pH=8.0, DO=5.0)
    assert c.pH_setpoint == c.limits.pH_setpoint[1]
    assert c.DO_setpoint == c.limits.DO_setpoint[0]
