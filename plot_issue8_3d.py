"""Issue #8 in 3D: no plate, measured tube profiles, convection from the outer face only.
Outer shield peak vs total stack length (gaps plus shield thicknesses), under 20 mm."""
import time
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from radiation import Setup
from radiation3d import load_tube_profiles, solve3d

s = Setup()
profiles = load_tube_profiles("astral_temps.csv")
SPACINGS = (0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.010, 0.015)
grid = dict(n_x=10, n_z=16, n_segments=30)

fig, axes = plt.subplots(1, len(profiles), figsize=(6 * len(profiles), 5), sharey=True)
for ax, (name, prof) in zip(axes, profiles.items()):
    t0 = time.time()
    for sp in SPACINGS:
        L, T, N = [], [], []
        n = 1
        while n * (sp + s.shield_thickness) < 0.020:
            r = solve3d(s, prof, n, plate=False, spacing=sp, convect_front=False, **grid)
            L.append(n * (sp + s.shield_thickness) * 1e3); T.append(r.T_plate_max - 273.15); N.append(n)
            n += 1
        ax.plot(L, T, "o-", lw=0.8, ms=4, label=f"{sp*1e3:g} mm spacing")
        for l, t, k in zip(L, T, N):
            ax.annotate(str(k), (l, t), textcoords="offset points", xytext=(3, 3), fontsize=6)
        print(name, f"{sp*1e3:g} mm", [f"{l:.1f}:{t:.0f}" for l, t in zip(L, T)], flush=True)
    print(name, f"{time.time() - t0:.0f} s", flush=True)
    ax.axhline(80, color="k", ls="--", lw=1, label="80 °C limit")
    ax.set_xlim(0, 20)
    ax.set_title(f'"{name}" profile')
    ax.set_xlabel(f"Total length incl. {s.shield_thickness*1e3:g} mm shields (mm)")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("Outer shield peak (°C)")
axes[0].legend(fontsize=8)
fig.suptitle("3D, no plate, convection from outer face only (labels = number of shields)")
fig.tight_layout(); fig.savefig("issue8_3d_length_vs_temp.png", dpi=130)
