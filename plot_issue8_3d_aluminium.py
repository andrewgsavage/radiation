"""Issue #8 in 3D, hottest ("545") profile: aluminium shields (k = 200 W/m.K) at
emissivity 0.1 / 0.2 / 0.4 and thickness 1 / 0.5 mm. Uses the cache in plot_issue8_3d."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_issue8_3d import SPACINGS, peak

NAME = "545°C"
K_AL = 200.0
EPS = (0.1, 0.2, 0.4)
THICK = (0.001, 0.0005)
fig, axes = plt.subplots(len(THICK), len(EPS), figsize=(6 * len(EPS), 5 * len(THICK)), sharey=True)
for i, t in enumerate(THICK):
    for ax, eps in zip(axes[i], EPS):
        for sp in SPACINGS:
            L, T, N = [], [], []
            n = 1
            while n * (sp + t) < 0.020:
                L.append(n * (sp + t) * 1e3); T.append(peak(eps, NAME, sp, n, t, K_AL)); N.append(n)
                n += 1
            ax.plot(L, T, "o-", lw=0.8, ms=4, label=f"{sp*1e3:g} mm spacing")
            for l, v, k in zip(L, T, N):
                ax.annotate(str(k), (l, v), textcoords="offset points", xytext=(3, 3), fontsize=6)
            print(t, eps, f"{sp*1e3:g} mm", flush=True)
        ax.axhline(80, color="k", ls="--", lw=1, label="80 °C limit")
        ax.set_xlim(0, 20)
        ax.set_title(f'"{NAME}", aluminium, ε = {eps}, {t*1e3:g} mm thick')
        ax.set_xlabel("Total length incl. shield thickness (mm)")
        ax.grid(alpha=0.3)
    axes[i][0].set_ylabel("Outer shield peak (°C)")
axes[0][0].legend(fontsize=8)
fig.suptitle("3D, no plate, aluminium shields, convection from outer face only (labels = number of shields)")
fig.tight_layout(); fig.savefig("issue8_3d_aluminium.png", dpi=110)
