# radiation

Radiation heat transfer from a 650 °C tube to a nearby vertical plate, and how many
heat shields keep the plate below 150 °C ([issue #1](https://github.com/andrewgsavage/radiation/issues/1)),
and how many keep the outermost shield below 150 °C when there's no plate
([issue #3](https://github.com/andrewgsavage/radiation/issues/3)).

- [`issue1_plate_temperature.ipynb`](issue1_plate_temperature.ipynb): the worked problem,
  with geometry drawings, results and a sensitivity study
- [`radiation.py`](radiation.py): view factors, the radiation and conduction network, and
  natural-convection correlations
- [`radiation3d.py`](radiation3d.py): the 3D version, for a tube whose temperature varies
  along its length; sheets are split into x–z patches
- [`astral_temps.csv`](astral_temps.csv): measured temperatures (°C) against position x (mm),
  one column per case. The column headers (545, 500 and 300 °C) are external pyrometer
  readings of the tube, which carry an emissivity calibration error, so they sit below the
  measured values. Every column is readings on the tube itself, symmetric about the tube
  centre at x = 200 mm. Blank cells are positions with no reading.
- [`tests/`](tests/): checks against textbook view factors and the 1/(N + 1) shield result

## Results (base case)

| | Plate peak | Plate mean |
|---|---|---|
| No shields | 188 °C | 119 °C |
| One ε = 0.2 shield | 61 °C | 52 °C |

**One shield is needed.** The plate's hottest band is directly opposite the tube. A single
plate-temperature estimate (about 122 °C) hides it. Assumptions are listed in the notebook:
tube and plate ε = 0.8, carbon-steel plate, 25 °C still air, and heat carried by hot air
rising off the tube ignored.

## No plate (issue #3)

| Shields | Outer shield peak |
|---|---|
| 1 | 131 °C |
| 2 | 83 °C |

**One shield** is enough in the base case, but only just: a black tube, ε = 0.3 shields, or
0.1 mm shields each push a lone shield over 150 °C. Two shields give a safe margin.

## 3D, measured tube temperatures

Uses the profiles in `astral_temps.csv`, with the plate centred on the tube (peak, °C):

| Case | Plate, no shields | Plate, 1 shield | No plate, 1 shield | No plate, 2 shields |
|---|---|---|---|---|
| "545" | 130 | 45 | 98 | 61 |
| "500" | 119 | 43 | 91 | 57 |
| "300" | 64 | 33 | 51 | 39 |

The measured tube is only near its peak over the middle ~50 mm, so **no shields are needed
for the plate, and with no plate one shield is enough**. The plate's margin is thin in the
"545" case: a black tube or 50 °C surroundings bring it to just under 150 °C, and a black
stainless plate needs one shield. With no plate, an oxidised (ε = 0.45) 0.1 mm shield needs
a second one.

## Outer-face convection only (issue #8)

2D model, tube at a constant 650 °C, and only the face furthest from the tube loses heat
to room air by convection (`solve(..., convect_front=False)`). Run `python study_issue8.py`.
Peak of the outermost sheet, °C:

| Every gap (tube→sheet, sheet→sheet) | Plate, 0 / 1 / 2 shields | No plate, 1 / 2 / 3 shields | Shields needed (plate / no plate) |
|---|---|---|---|
| 2 mm | 244 / 95 / 81 | 205 / 142 / 112 | 1 / 2 |
| 5 mm | 230 / 82 / 64 | 189 / 114 / 82 | 1 / 2 |
| 10 mm | 210 / 69 / 50 | 168 / 89 / 60 | 1 / 2 |
| 20 mm | 181 / 56 / 38 | 139 / 65 / 43 | 1 / 1 |
| 40 mm | 142 / 42 / 30 | 105 / 46 / 32 | 0 / 1 |

**With a plate, one shield; with no plate, two shields** at spacings up to 10 mm. Losing
the tube-side convection raises temperatures (plate with no shields: 210 °C, against 188 °C
in issue #1). Wider spacing helps because the sheets are finite (140 × 220 mm): more of
each sheet's radiation escapes past the edges instead of reaching the next sheet. Squeezing
more shields into a fixed 10 mm gap gives diminishing returns (no plate: 168, 114, 96, 87 °C
for 1–4 shields), since the gaps thin out and air conduction across them takes over.

With an **80 °C limit** and the stack kept under 20 mm from the tube surface, counting
the 0.5 mm thickness of each stainless shield (`python plot_issue8.py`): the fewest shields
is 4 at 4 mm spacing (18 mm, 72 °C), and the shortest stack is 9 at 1 mm spacing
(13.5 mm, 77 °C). Packing shields closer than 1 mm makes the stack longer, because
the shields' own thickness then outweighs the gaps.

### Issue #8 in 3D, measured tube temperatures

Same setup (no plate, outer-face convection only, 80 °C limit, under 20 mm including 0.5 mm
shields), with the `astral_temps.csv` profiles (`python plot_issue8_3d.py`, ~10 min).
Because the tube is only near its peak over the middle ~50 mm, far fewer shields are needed than in 2D:

| Profile | Fewest shields | Shortest stack |
|---|---|---|
| "545" | 2 at 7.5 mm (16 mm, 71 °C); 2 at 5 mm is marginal (79 °C) | 5 at 0.5 mm (5 mm, 77 °C) |
| "500" | 2 at 4 mm (9 mm, 78 °C) | 4 at 0.5 mm (4 mm, 78 °C) |
| "300" | 1 at any spacing (70 °C or less) | 1 at 0.5 mm (1 mm, 70 °C) |

With **ε = 0.4 shields** (oxidised stainless), the second row of the plot:

| Profile | Fewest shields | Shortest stack |
|---|---|---|
| "545" | 3 at 4 mm (13.5 mm, 77 °C) | 8 at 0.5 mm (8 mm, 78 °C) |
| "500" | 3 at 3 mm (10.5 mm, 79 °C) | 7 at 0.5 mm (7 mm, 77 °C) |
| "300" | 1 at 10 mm (10.5 mm, 78 °C) | 2 at 0.5 mm (2 mm, 75 °C) |

Every solved point is in `issue8_3d_results.csv`.

**Half-thickness (0.25 mm) shields**, "545" profile only (`python plot_issue8_3d_thickness.py`):

| Shield ε | Thickness | Fewest shields | Shortest stack |
|---|---|---|---|
| 0.2 | 0.5 mm | 2 at 5 mm (11 mm, 79 °C) | 5 at 0.5 mm (5 mm, 77 °C) |
| 0.2 | 0.25 mm | 3 at 4 mm (12.75 mm, 74 °C) | 7 at 0.5 mm (5.25 mm, 78 °C) |
| 0.4 | 0.5 mm | 3 at 4 mm (13.5 mm, 77 °C) | 8 at 0.5 mm (8 mm, 78 °C) |
| 0.4 | 0.25 mm | 4 at 3 mm (13 mm, 78 °C) | 10 at 0.5 mm (7.5 mm, 78 °C) |

Thinner shields run hotter: they conduct less heat sideways, so the hot spot opposite the
tube's peak spreads less. This costs about one extra shield, and saves very little length.

## Running

```bash
python -m venv .venv
```

```bash
.venv/Scripts/python -m pip install -r requirements.txt
```

```bash
.venv/Scripts/python -m pytest tests
```
