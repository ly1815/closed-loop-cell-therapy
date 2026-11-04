"""Figures for the control refresher. Run from this folder:

    python make_figures.py

Needs numpy, scipy and matplotlib (already dependencies of the package).
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy import signal

OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False,
                     "axes.grid": True, "grid.alpha": 0.3, "legend.frameon": False})
C = ["#0F6E56", "#C4532E", "#2F6DB5", "#7A4FA0", "#555555"]


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=200)
    plt.close(fig)


# 1. First-order step responses for several time constants
t = np.linspace(0, 20, 500)
fig, ax = plt.subplots(figsize=(6.0, 2.8))
for tau, c in zip((1, 2, 4), C):
    ax.plot(t, 2 * (1 - np.exp(-t / tau)), color=c, label=fr"$\tau$ = {tau}")
    ax.plot(tau, 2 * (1 - np.exp(-1)), "o", color=c, ms=4)
ax.axhline(2, color="k", lw=0.7, ls="--")
ax.axhline(2 * 0.632, color="k", lw=0.5, ls=":")
ax.text(14.5, 2 * 0.632 + 0.05, "63.2% of the final change", fontsize=8)
ax.text(14.5, 2.05, "final value = K = 2", fontsize=8)
ax.set_xlabel("time")
ax.set_ylabel("output y")
ax.set_ylim(0, 2.3)
ax.legend(loc="lower right")
save(fig, "first_order.png")

# 2. Second-order step responses for several damping ratios
t = np.linspace(0, 20, 800)
fig, ax = plt.subplots(figsize=(6.0, 2.8))
for zeta, c in zip((0.2, 0.5, 1.0, 2.0), C):
    sys = signal.lti([1], [1, 2 * zeta, 1])
    _, y = signal.step(sys, T=t)
    ax.plot(t, y, color=c, label=fr"$\zeta$ = {zeta}")
ax.axhline(1, color="k", lw=0.7, ls="--")
ax.set_xlabel(r"time (in units of $\tau$)")
ax.set_ylabel("output y")
ax.legend(loc="lower right")
save(fig, "second_order.png")

# 3. P versus PI on a first-order plant: setpoint step at t=0, disturbance step at t=40
K, tau = 2.0, 5.0
dt = 0.01
t = np.arange(0, 80, dt)


def run_loop(kc, ti=None):
    y, integ, ys, us = 0.0, 0.0, [], []
    for tk in t:
        r = 1.0
        d = 0.5 if tk >= 40 else 0.0
        e = r - y
        if ti is not None:
            integ += kc * e * dt / ti
        u = kc * e + integ
        y += dt / tau * (-y + K * u + d)
        ys.append(y)
        us.append(u)
    return np.array(ys), np.array(us)


fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
for (label, kc, ti), c in zip([("P only, $K_c$ = 1", 1.0, None), ("PI, $K_c$ = 1, $\\tau_I$ = 5", 1.0, 5.0)], C):
    y, u = run_loop(kc, ti)
    axes[0].plot(t, y, color=c, label=label)
    axes[1].plot(t, u, color=c)
axes[0].axhline(1, color="k", lw=0.7, ls="--")
axes[0].annotate("offset", xy=(30, 0.667), xytext=(18, 0.35), fontsize=8,
                 arrowprops=dict(arrowstyle="->", lw=0.7))
axes[0].axvline(40, color="grey", lw=0.6, ls=":")
axes[0].text(41, 0.05, "disturbance\nstep", fontsize=7)
axes[0].set_ylabel("output y (setpoint = 1)")
axes[1].set_ylabel("controller output u")
for a in axes:
    a.set_xlabel("time")
axes[0].legend(loc="center right", fontsize=7.5)
save(fig, "p_vs_pi.png")

# 4. P control of a third-order plant: increasing the gain leads to instability
t = np.linspace(0, 30, 1500)
fig, ax = plt.subplots(figsize=(6.0, 2.8))
for kc, c in zip((1, 4, 8, 9), C):
    # closed loop: Kc / ((s+1)^3 + Kc)
    sys = signal.lti([kc], np.polyadd(np.poly([-1, -1, -1]), [kc]))
    _, y = signal.step(sys, T=t)
    ax.plot(t, y, color=c, label=fr"$K_c$ = {kc}")
ax.set_ylim(-1.5, 3.5)
ax.set_xlabel("time")
ax.set_ylabel("output y (setpoint = 1)")
ax.legend(loc="upper left", ncol=4)
save(fig, "instability.png")

# 5. Windup: PI with a saturating actuator, with and without anti-windup
K, tau = 1.0, 10.0
kc, ti = 5.0, 10.0
umin, umax = 0.0, 1.5
t = np.arange(0, 60, dt)


def run_sat(anti_windup):
    y, integ, ys, us = 0.0, 0.0, [], []
    for tk in t:
        r = 1.0 if tk < 30 else 0.3
        e = r - y
        integ += kc * e * dt / ti
        u_raw = kc * e + integ
        u = min(max(u_raw, umin), umax)
        if anti_windup:
            integ -= u_raw - u
        y += dt / tau * (-y + K * u)
        ys.append(y)
        us.append(u)
    return np.array(ys), np.array(us)


fig, axes = plt.subplots(1, 2, figsize=(7.0, 2.7))
for (label, aw), c in zip([("no anti-windup", False), ("with anti-windup (clamping)", True)], C[1::-1]):
    y, u = run_sat(aw)
    axes[0].plot(t, y, color=c, label=label)
    axes[1].plot(t, u, color=c)
axes[0].plot(t, np.where(t < 30, 1.0, 0.3), "k--", lw=0.7, label="setpoint")
axes[1].axhline(umax, color="grey", lw=0.6, ls=":")
axes[1].axhline(umin, color="grey", lw=0.6, ls=":")
axes[0].set_ylabel("output y")
axes[1].set_ylabel("actuator u (limits 0 to 1.5)")
for a in axes:
    a.set_xlabel("time")
fig.legend(*axes[0].get_legend_handles_labels(), loc="upper center", ncol=3, fontsize=7.5)
fig.tight_layout(rect=(0, 0, 1, 0.9))
fig.savefig(OUT / "windup.png", dpi=200)
plt.close(fig)

# 6. Poles in the complex plane and the responses they produce
fig, axes = plt.subplots(1, 2, figsize=(7.0, 3.1))
poles = [(-1.0, 0.0, C[0], "real, negative: smooth decay"),
         (-0.3, 1.5, C[2], "complex, negative real part: decaying oscillation"),
         (0.15, 1.0, C[1], "positive real part: growing (unstable)")]
t = np.linspace(0, 15, 600)
for re, im, c, label in poles:
    axes[0].plot(re, im, "x", color=c, ms=8, mew=2)
    if im:
        axes[0].plot(re, -im, "x", color=c, ms=8, mew=2)
        y = np.exp(re * t) * np.cos(im * t)
    else:
        y = np.exp(re * t)
    axes[1].plot(t, y, color=c, label=label)
axes[0].axvline(0, color="k", lw=0.8)
axes[0].axhline(0, color="k", lw=0.8)
axes[0].axvspan(-1.4, 0, color=C[0], alpha=0.07)
axes[0].text(-1.3, 1.75, "stable\n(left half)", fontsize=8, color=C[0])
axes[0].text(0.05, 1.75, "unstable", fontsize=8, color=C[1])
axes[0].set_xlim(-1.4, 0.6)
axes[0].set_ylim(-2.2, 2.2)
axes[0].set_xlabel("real part of the pole")
axes[0].set_ylabel("imaginary part")
axes[1].set_ylim(-2.5, 3)
axes[1].set_xlabel("time")
axes[1].set_ylabel("response term")
fig.legend(*axes[1].get_legend_handles_labels(), loc="lower center", ncol=3, fontsize=7.5)
fig.tight_layout(rect=(0, 0.08, 1, 1))
fig.savefig(OUT / "poles.png", dpi=200)
plt.close(fig)

print("figures written to", OUT)
