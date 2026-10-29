"""Run the FEA mesh-convergence check and case matrix; writes fea_results.json."""

import json
import time

import radiation as rad
import radiation3d as r3
import fea3d

C = 273.15
s = rad.Setup()
profiles = r3.load_tube_profiles("astral_temps.csv")
out = {"convergence": [], "cases": []}

for n_x, n_z, n_y in [(10, 15, 2), (20, 30, 2), (20, 30, 4), (30, 45, 2), (40, 60, 2)]:
    for n in (0, 1):
        if n == 1 and n_x == 40:
            continue
        t = time.time()
        f = fea3d.solve_fea(s, profiles["545°C"], n, n_x=n_x, n_z=n_z, n_y=n_y)
        row = dict(n_x=n_x, n_z=n_z, n_y=n_y, shields=n, peak=f.T_plate_max - C, mean=f.T_plate_mean - C,
                   seconds=time.time() - t)
        out["convergence"].append(row)
        print(row, flush=True)

for name, profile in profiles.items():
    for plate, counts in ((True, (0, 1, 2)), (False, (1, 2, 3))):
        for n in counts:
            f = fea3d.solve_fea(s, profile, n, plate=plate)
            p = r3.solve3d(s, profile, n, plate=plate)
            row = dict(profile=name, plate=plate, shields=n, peak=f.T_plate_max - C, mean=f.T_plate_mean - C,
                       hottest_sheet=max(T.max() for T in f.T) - C, patch_peak=p.T_plate_max - C)
            out["cases"].append(row)
            print(row, flush=True)

with open("fea_results.json", "w", encoding="utf-8") as fh:
    json.dump(out, fh, indent=1, ensure_ascii=False)
