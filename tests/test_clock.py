import time

from cltherapy.sim.clock import SimClock


def test_clock_runs_pauses_and_jumps():
    c = SimClock(speed_h_per_s=3600.0)  # one simulated hour per millisecond
    time.sleep(0.01)
    assert c.now_h() > 5.0
    c.pause()
    t = c.now_h()
    time.sleep(0.01)
    assert c.now_h() == t
    c.jump_to(48.0)
    assert c.now_h() == 48.0
    c.resume()
    time.sleep(0.005)
    assert c.now_h() > 48.0
