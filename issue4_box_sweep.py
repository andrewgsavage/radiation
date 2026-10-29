"""Issue #4: plate and shield peak temperatures inside a box of still air.

The box's air and walls are swept from 150 to 350 C in 10 C steps, for up to ten
0.1 mm stainless shields. Prints one Markdown table per shield emissivity; run with
`python issue4_box_sweep.py`.
"""

import sys

import numpy as np

import radiation as rad

C = 273.15
T_box = np.arange(150, 351, 10) + C
N = (1, 2, 3, 4, 6, 8, 10)
EPS = (0.2, 0.45)

sys.stdout.reconfigure(encoding="utf-8")
rows = rad.box_sweep(T_box, n_shields=(0,) + N, eps_shield=EPS,
                     setup=rad.Setup(shield_thickness=0.1e-3))
get = {(r["T_box"], r["n_shields"], r["eps_shield"] if r["n_shields"] else None): r for r in rows}

for e in EPS:
    print(f"\nShield ε {e}: hottest shield / plate peak, °C\n")
    head = ["Box", "0"] + [str(n) for n in N]
    print("| " + " | ".join(head) + " |")
    print("|" + "---|" * len(head))
    for Tb in T_box:
        cells = [f"{Tb - C:.0f}", f"– / {get[(Tb, 0, None)]['T_plate_max'] - C:.0f}"]
        for n in N:
            r = get[(Tb, n, e)]
            cells.append(f"{r['T_shield_max'] - C:.0f} / {r['T_plate_max'] - C:.0f}")
        print("| " + " | ".join(cells) + " |")
