"""Issue #8 in 3D, hottest ("545") profile: 15 um aluminium foil shields (k = 200 W/m.K)
at emissivity 0.05 / 0.1 / 0.2. Uses the cache in plot_issue8_3d. Shield count capped at
N_MAX, since foil allows very many shields in 20 mm."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from plot_issue8_3d import SPACINGS, peak

NAME = "545°C"
K_AL = 200.0
T_FOIL = 15e-6
EPS = (0.05, 0.1, 0.2)
N_MAX = 15
fig, axes = plt.subplots(1, len(EPS), figsize=(6 * len(EPS), 5), sharey=True)
for ax, eps in zip(axes, EPS):
    for sp in SPACINGS:
        L, T, N = [], [], []
        n = 1
        while n * (sp + T_FOIL) < 0.020 and n <= N_MAX:
            L.append(n * (sp + T_FOIL) * 1e3); T.append(peak(eps, NAME, sp, n, T_FOIL, K_AL)); N.append(n)
            if T[-1] < 40:
                break
            n += 1
        ax.plot(L, T, "o-", lw=0.8, ms=4, label=f"{sp*1e3:g} mm spacing")
        for l, v, k in zip(L, T, N):
            ax.annotate(str(k), (l, v), textcoords="offset points", xytext=(3, 3), fontsize=6)
        print(eps, f"{sp*1e3:g} mm", [f"{l:.1f}:{v:.0f}" for l, v in zip(L, T)], flush=True)
    ax.axhline(80, color="k", ls="--", lw=1, label="80 °C limit")
    ax.set_xlim(0, 20)
    ax.set_title(f'"{NAME}", 15 µm aluminium foil, ε = {eps}')
    ax.set_xlabel("Total length incl. foil thickness (mm)")
    ax.grid(alpha=0.3)
axes[0].set_ylabel("Outer shield peak (°C)")
axes[0].legend(fontsize=8)
fig.suptitle("3D, no plate, foil shields, convection from outer face only (labels = number of shields)")
fig.tight_layout(); fig.savefig("issue8_3d_foil.png", dpi=110)
