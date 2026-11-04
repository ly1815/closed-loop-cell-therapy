"""The generic PI in control/pid.py, kept as the worked example of the control refresher notes."""

from cltherapy.control.pid import PI


def test_pi_holds_first_order_setpoint():
    pi = PI(kp=2.0, ti_h=0.5, setpoint=1.0, out_min=0.0, out_max=5.0)
    y, dt = 0.0, 0.01
    for _ in range(2000):
        u = pi.update(y, dt)
        y += dt * (u - y) / 0.5
    assert abs(y - 1.0) < 1e-3


def test_pi_anti_windup_recovers_quickly():
    pi = PI(kp=1.0, ti_h=0.2, setpoint=10.0, out_min=0.0, out_max=1.0)
    for _ in range(500):
        pi.update(0.0, 0.01)
    pi.setpoint = 0.0
    assert pi.update(0.5, 0.01) == 0.0
