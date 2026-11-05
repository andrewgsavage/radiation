"""Issue #8, no plate: outermost shield peak vs total stack length (tube surface to last shield)."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from radiation import Setup, solve

s = Setup()
fig, ax = plt.subplots(figsize=(8, 5))
for sp in (0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.010, 0.015):
    L, T = [], []
    for n in range(1, 20):
        if n * sp >= 0.020:  # total length under 20 mm
            break
        L.append(n * sp * 1e3)  # n gaps: tube->shield 1, then shield to shield
        T.append(solve(s, n, plate=False, spacing=sp, convect_front=False, n_strips=40).T_plate_max - 273.15)
    ax.plot(L, T, "o-", lw=0.8, label=f"{sp*1e3:g} mm spacing")
    for l, t, n in zip(L, T, range(1, 20)):
        ax.annotate(str(n), (l, t), textcoords="offset points", xytext=(4, 3), fontsize=7)
ax.axhline(80, color="k", ls="--", lw=1, label="80 °C limit")
ax.set_xlim(0, 20)
ax.set_xlabel("Total length, tube surface to outer shield (mm)")
ax.set_ylabel("Outer shield peak (°C)")
ax.set_title("No plate, convection from outer face only (labels = number of shields)")
ax.legend(); ax.grid(alpha=0.3, )
fig.tight_layout(); fig.savefig("issue8_length_vs_temp.png", dpi=130)
