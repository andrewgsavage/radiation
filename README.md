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
