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

No plate, and only the face furthest from the tube loses heat to room air by convection
(`convect_front=False`). Shields are 0.5 mm stainless (k = 16 W/m·K) unless stated. The limit is
**80 °C on the outermost shield**, and the total length (tube surface to outer face, gaps plus
shield thickness) is kept **under 20 mm**. Air in the gaps between shields is still air
(conduction only), with k_air at each gap's own mean temperature; the gap between the tube
and the first shield is radiation only. Every 3D point solved is in `issue8_3d_results.csv`
(`issue8_3d_results_fixed_kair.csv` holds the earlier runs with k_air fixed at 200 °C);
`python summarise_issue8.py` prints the tables below from it.

**2D, tube at a constant 650 °C** (`python plot_issue8.py`): fewest shields 4 at 3 mm
(14 mm, 79 °C); shortest stack 8 at 1 mm (12 mm, 79 °C). `python study_issue8.py` also has
the plate cases and wider spacings.

**3D, measured tube profiles**, 0.5 mm stainless (`python plot_issue8_3d.py`):

| Profile | Shield ε | Fewest shields | Shortest stack |
|---|---|---|---|
| "545" | 0.2 | 2 at 5 mm (11 mm, 77 °C) | 5 at 0.5 mm (5 mm, 75 °C) |
| "545" | 0.4 | 3 at 4 mm (13.5 mm, 74 °C) | 7 at 0.5 mm (7 mm, 79 °C) |
| "500" | 0.2 | 2 at 4 mm (9 mm, 76 °C) | 4 at 0.5 mm (4 mm, 76 °C) |
| "500" | 0.4 | 3 at 3 mm (10.5 mm, 76 °C) | 6 at 0.5 mm (6 mm, 80 °C) |
| "300" | 0.2 | 1 at any spacing | 1 at 0.5 mm (1 mm, 70 °C) |
| "300" | 0.4 | 1 at 10 mm (10.5 mm, 78 °C) | 2 at 0.5 mm (2 mm, 74 °C) |

**Shield material and thickness**, "545" profile only (`plot_issue8_3d_thickness.py`,
`plot_issue8_3d_aluminium.py`, `plot_issue8_3d_foil.py`; aluminium k = 200 W/m·K):

| Shield | ε | Fewest shields | Shortest stack |
|---|---|---|---|
| Stainless 1 mm | 0.2 | 2 at 2 mm (6 mm, 79 °C) | 3 at 0.5 mm (4.5 mm, 78 °C) |
| Stainless 1 mm | 0.4 | 2 at 7.5 mm (17 mm, 76 °C) | 4 at 1 mm (8 mm, 78 °C) |
| Stainless 0.25 mm | 0.2 | 3 at 3 mm (9.75 mm, 78 °C) | 7 at 0.5 mm (5.25 mm, 75 °C) |
| Stainless 0.25 mm | 0.4 | 4 at 3 mm (13 mm, 76 °C) | 9 at 0.5 mm (6.75 mm, 80 °C) |
| Aluminium 1 mm | 0.1 | 1 (1.5 mm, 57 °C) | same |
| Aluminium 1 mm | 0.2 | 1 (1.5 mm, 74 °C) | same |
| Aluminium 1 mm | 0.4 | 2 at 1 mm (4 mm, 80 °C) | same |
| Aluminium 0.5 mm | 0.1 | 1 (1 mm, 60 °C) | same |
| Aluminium 0.5 mm | 0.2 | 1 at 1 mm (1.5 mm, 80 °C, marginal) | same |
| Aluminium 0.5 mm | 0.4 | 2 at 2 mm (5 mm, 77 °C) | 3 at 0.5 mm (3 mm, 80 °C) |
| Aluminium foil 15 µm | 0.05 | 1 at 5 mm (5 mm, 79 °C) | 2 at 0.5 mm (1 mm, 72 °C) |
| Aluminium foil 15 µm | 0.1 | 2 at 4 mm (8 mm, 79 °C) | 4 at 0.5 mm (2 mm, 78 °C) |
| Aluminium foil 15 µm | 0.2 | 3 at 4 mm (12 mm, 76 °C) | 8 at 0.5 mm (4.1 mm, 75 °C) |

What matters, in order: **emissivity**, then **sideways conduction** (k × thickness), which
spreads the hot spot opposite the tube's peak over the sheet. A single 1 mm aluminium shield
at ε ≤ 0.2 is enough. For a given length, closer spacing with more shields runs cooler,
until each shield's own thickness outweighs the gaps.

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
