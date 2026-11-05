"""Issue #8 in 3D: no plate, measured tube profiles, convection from the outer face only.
Outer shield peak vs total stack length (gaps plus shield thicknesses), under 20 mm,
for each shield emissivity. Results are cached in issue8_3d_results.csv."""
import csv
import dataclasses
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from radiation import Setup
from radiation3d import load_tube_profiles, solve3d

base = Setup()
profiles = load_tube_profiles("astral_temps.csv")
EPS = (0.2, 0.4)
SPACINGS = (0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.010, 0.015)
grid = dict(n_x=10, n_z=16, n_segments=30)
CACHE = "issue8_3d_results.csv"

cache = {}
if os.path.exists(CACHE):
    with open(CACHE, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            cache[(float(r["eps_shield"]), r["profile"], float(r["spacing_mm"]), int(r["n_shields"]))] = float(r["peak_C"])


def peak(eps, name, sp, n):
    key = (eps, name, round(sp * 1e3, 6), n)
    if key not in cache:
        s = dataclasses.replace(base, eps_shield=eps)
        r = solve3d(s, profiles[name], n, plate=False, spacing=sp, convect_front=False, **grid)
        cache[key] = r.T_plate_max - 273.15
        new = not os.path.exists(CACHE)
        with open(CACHE, "a", newline="", encoding="utf-8") as f:
            w = csv.writer(f)
            if new:
                w.writerow(["eps_shield", "profile", "spacing_mm", "n_shields", "peak_C"])
            w.writerow([*key, f"{cache[key]:.2f}"])
    return cache[key]


fig, axes = plt.subplots(len(EPS), len(profiles), figsize=(6 * len(profiles), 5 * len(EPS)), sharey=True)
for i, eps in enumerate(EPS):
    for ax, name in zip(axes[i], profiles):
        for sp in SPACINGS:
            L, T, N = [], [], []
            n = 1
            while n * (sp + base.shield_thickness) < 0.020:
                L.append(n * (sp + base.shield_thickness) * 1e3); T.append(peak(eps, name, sp, n)); N.append(n)
                n += 1
            ax.plot(L, T, "o-", lw=0.8, ms=4, label=f"{sp*1e3:g} mm spacing")
            for l, t, k in zip(L, T, N):
                ax.annotate(str(k), (l, t), textcoords="offset points", xytext=(3, 3), fontsize=6)
            print(eps, name, f"{sp*1e3:g} mm", [f"{l:.1f}:{t:.0f}" for l, t in zip(L, T)], flush=True)
        ax.axhline(80, color="k", ls="--", lw=1, label="80 °C limit")
        ax.set_xlim(0, 20)
        ax.set_title(f'"{name}" profile, shield ε = {eps}')
        ax.set_xlabel(f"Total length incl. {base.shield_thickness*1e3:g} mm shields (mm)")
        ax.grid(alpha=0.3)
    axes[i][0].set_ylabel("Outer shield peak (°C)")
axes[0][0].legend(fontsize=8)
fig.suptitle("3D, no plate, convection from outer face only (labels = number of shields)")
fig.tight_layout(); fig.savefig("issue8_3d_length_vs_temp.png", dpi=110)
