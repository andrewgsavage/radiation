"""Fewest shields and shortest stack (outer peak < 80 C, total length < 20 mm) per 3D case,
from the results caches. Usage: python summarise_issue8.py [csv ...]"""
import csv
import sys
from collections import defaultdict

for path in sys.argv[1:] or ["issue8_3d_results.csv"]:
    print(path)
    ok = defaultdict(list)
    for r in csv.DictReader(open(path, encoding="utf-8")):
        t, sp, n, T = float(r["shield_thickness_mm"]), float(r["spacing_mm"]), int(r["n_shields"]), float(r["peak_C"])
        key = (r["profile"], float(r["k_shield"]), t, float(r["eps_shield"]))
        ok[key]
        if T < 80 and n * (sp + t) < 20:
            ok[key].append((n, round(n * (sp + t), 2), sp, round(T, 1)))
    for key in sorted(ok):
        v = ok[key]
        print(" ", key, "fewest", min(v) if v else None, "shortest", min(v, key=lambda x: x[1]) if v else None)
