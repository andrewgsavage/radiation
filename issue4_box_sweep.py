"""Issue #4: plate and shield peak temperatures inside a box of still air.

The box's air and walls are swept from 150 to 350 C in 10 C steps. Prints a
Markdown table; run with `python issue4_box_sweep.py`.
"""

import numpy as np

import radiation as rad

C = 273.15
T_box = np.arange(150, 351, 10) + C
N = (1, 2, 3)
EPS = (0.2, 0.45)

rows = rad.box_sweep(T_box, n_shields=(0,) + N, eps_shield=EPS)
get = {(r["T_box"], r["n_shields"], r["eps_shield"] if r["n_shields"] else None): r for r in rows}

head = ["Box air/walls", "No shield: plate"]
head += [f"{n} shield{'s' * (n > 1)}, ε {e}: shield / plate" for e in EPS for n in N]
print("| " + " | ".join(head) + " |")
print("|" + "---|" * len(head))
for Tb in T_box:
    cells = [f"{Tb - C:.0f}", f"{get[(Tb, 0, None)]['T_plate_max'] - C:.0f}"]
    for e in EPS:
        for n in N:
            r = get[(Tb, n, e)]
            cells.append(f"{r['T_shield_max'] - C:.0f} / {r['T_plate_max'] - C:.0f}")
    print("| " + " | ".join(cells) + " |")
