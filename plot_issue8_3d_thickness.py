"""Issue #8 in 3D, hottest ("545") profile only: shield emissivity 0.2 / 0.4 against
shield thickness 1 / 0.5 / 0.25 mm. Uses the cache in plot_issue8_3d."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_issue8_3d import EPS, SPACINGS, peak

NAME = "545°C"
THICK = (0.001, 0.0005, 0.00025)
fig, axes = plt.subplots(len(EPS), len(THICK), figsize=(7 * len(THICK), 5 * len(EPS)), sharey=True)
for i, eps in enumerate(EPS):
    for ax, t in zip(axes[i], THICK):
        for sp in SPACINGS:
            L, T, N = [], [], []
            n = 1
            while n * (sp + t) < 0.020:
                L.append(n * (sp + t) * 1e3); T.append(peak(eps, NAME, sp, n, t)); N.append(n)
                n += 1
            ax.plot(L, T, "o-", lw=0.8, ms=4, label=f"{sp*1e3:g} mm spacing")
            for l, v, k in zip(L, T, N):
                ax.annotate(str(k), (l, v), textcoords="offset points", xytext=(3, 3), fontsize=6)
            print(eps, t, f"{sp*1e3:g} mm", flush=True)
        ax.axhline(80, color="k", ls="--", lw=1, label="80 °C limit")
        ax.set_xlim(0, 20)
        ax.set_title(f'"{NAME}" profile, shield ε = {eps}, {t*1e3:g} mm thick')
        ax.set_xlabel("Total length incl. shield thickness (mm)")
        ax.grid(alpha=0.3)
    axes[i][0].set_ylabel("Outer shield peak (°C)")
axes[0][0].legend(fontsize=8)
fig.suptitle("3D, no plate, convection from outer face only (labels = number of shields)")
fig.tight_layout(); fig.savefig("issue8_3d_thickness.png", dpi=110)
