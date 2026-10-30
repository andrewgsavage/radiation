"""FEA case matrix with the box interior (tube side) at raised temperatures and the
outermost sheet's back face open to 30 C outside air; writes fea_ambient_results.json."""

import json
from dataclasses import replace

import radiation as rad
import radiation3d as r3
import fea3d

C = 273.15
profiles = r3.load_tube_profiles("astral_temps.csv")
rows = []
for T_box in (150, 200, 250, 300):
    s = replace(rad.Setup(), T_amb=T_box + C)
    for name, profile in profiles.items():
        for plate, counts in ((True, (0, 1, 2)), (False, (1, 2, 3))):
            for n in counts:
                f = fea3d.solve_fea(s, profile, n, plate=plate, T_outside=30 + C)
                row = dict(T_box=T_box, profile=name, plate=plate, shields=n, peak=f.T_plate_max - C,
                           mean=f.T_plate_mean - C, hottest_sheet=max(T.max() for T in f.T) - C)
                rows.append(row)
                print(row, flush=True)

with open("fea_ambient_results.json", "w", encoding="utf-8") as fh:
    json.dump(rows, fh, indent=1, ensure_ascii=False)
