# radiation

Radiation heat transfer from a 650 °C tube to a nearby vertical plate, and how many
heat shields keep the plate below 150 °C ([issue #1](https://github.com/andrewgsavage/radiation/issues/1)),
and how many keep the outermost shield below 150 °C when there's no plate
([issue #3](https://github.com/andrewgsavage/radiation/issues/3)).

- [`issue1_plate_temperature.ipynb`](issue1_plate_temperature.ipynb): the worked problem,
  with geometry drawings, results and a sensitivity study
- [`radiation.py`](radiation.py): view factors, the radiation and conduction network, and
  natural-convection correlations
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

## 3D FEA with the measured tube temperatures

[`fea3d.py`](fea3d.py) is a scikit-fem finite-element model: solid hex meshes for the plate
and shields, and radiation from the tube split into 40 segments at the
[`astral_temps.csv`](astral_temps.csv) temperatures. Results and 3D contour plots are in
[`fea3d_astral.ipynb`](fea3d_astral.ipynb), and `python run_fea_cases.py` regenerates
`fea_results.json`. It matches the x-z patch model ([`radiation3d.py`](radiation3d.py)) to
within about 2 °C, and the mesh is converged to 0.2 °C.

| Case | Plate, no shield | Plate, 1 shield | No plate, 1 shield | No plate, 2 shields |
|---|---|---|---|---|
| 545 °C | 130 °C | 46 °C | 99 °C | 62 °C |
| 500 °C | 120 °C | 43 °C | 92 °C | 59 °C |
| 300 °C | 64 °C | 33 °C | 51 °C | 40 °C |

Peaks of the outermost sheet. With the measured tube temperatures, **the plate stays below
150 °C with no shields**, and with no plate **one shield** is enough. The tube is only near its
peak temperature close to x = 200 mm, so these runs are much cooler than the uniform 650 °C case.

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
