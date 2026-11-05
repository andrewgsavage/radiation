"""Issue #8: 2D model, constant tube temperature, convection only from the outermost face.
Shields needed to keep the outermost sheet below 150 C, and how sheet spacing affects it."""
from radiation import Setup, solve

s = Setup()
LIMIT = 150 + 273.15
kw = dict(convect_front=False, n_strips=40)

print("Fixed 10 mm gap, shields evenly spaced inside it (peak of outermost sheet, C)")
for plate in (True, False):
    row = [f"{solve(s, n, plate=plate, **kw).T_plate_max - 273.15:5.0f}" for n in range(0 if plate else 1, 7)]
    print(" plate" if plate else " no plate", row)

print("\nFixed pitch: every gap = spacing (stack grows with shields)")
for plate in (True, False):
    print(" plate" if plate else " no plate")
    for sp in (0.002, 0.005, 0.010, 0.020, 0.040):
        peaks, need = [], None
        for n in range(0 if plate else 1, 9):
            T = solve(s, n, plate=plate, spacing=sp, **kw).T_plate_max
            peaks.append(f"{T - 273.15:4.0f}")
            if need is None and T < LIMIT:
                need = n
        print(f"  {sp*1e3:4.0f} mm  needed={need}  peaks={peaks}")
