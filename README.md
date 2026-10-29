# radiation

Radiation heat transfer from a 650 °C tube to a nearby vertical plate, and how many
heat shields keep the plate below 150 °C ([issue #1](https://github.com/andrewgsavage/radiation/issues/1)).

- [`issue1_plate_temperature.ipynb`](issue1_plate_temperature.ipynb): the worked problem,
  with geometry drawings, results and a sensitivity study
- [`radiation.py`](radiation.py): view factors, the radiation and conduction network, and
  natural-convection correlations
- [`issue4_box_sweep.py`](issue4_box_sweep.py): plate and shield temperatures inside a
  hot box of still air (issue #4)
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

## Inside a box of still air ([issue #4](https://github.com/andrewgsavage/radiation/issues/4))

[`issue4_box_sweep.py`](issue4_box_sweep.py) puts the tube, shields and plate inside a
box whose air and walls are both at a set temperature, swept from 150 to 350 °C. The
sheets radiate to the box walls, treated as black because a box much larger than the
sheets acts as a cavity, and lose heat to the box air by natural convection. The tube
is held at 650 °C. Shields are 0.1 mm stainless at ε 0.2 (bright) or 0.45 (oxidized).
The table gives peak temperatures in °C.

| Box air/walls | No shield: plate | 1 shield, ε 0.2: shield / plate | 2 shields, ε 0.2: shield / plate | 3 shields, ε 0.2: shield / plate | 1 shield, ε 0.45: shield / plate | 2 shields, ε 0.45: shield / plate | 3 shields, ε 0.45: shield / plate |
|---|---|---|---|---|---|---|---|
| 150 | 266 | 260 / 175 | 272 / 168 | 276 / 165 | 316 / 197 | 338 / 183 | 346 / 176 |
| 160 | 272 | 267 / 184 | 280 / 177 | 284 / 174 | 322 / 206 | 343 / 191 | 351 / 185 |
| 170 | 279 | 275 / 193 | 288 / 186 | 292 / 183 | 327 / 214 | 349 / 200 | 357 / 194 |
| 180 | 285 | 283 / 202 | 295 / 195 | 299 / 193 | 333 / 223 | 354 / 209 | 362 / 203 |
| 190 | 292 | 290 / 211 | 303 / 205 | 307 / 202 | 339 / 231 | 360 / 218 | 368 / 212 |
| 200 | 298 | 298 / 220 | 311 / 214 | 315 / 211 | 344 / 240 | 365 / 226 | 373 / 221 |
| 210 | 305 | 305 / 230 | 318 / 223 | 322 / 221 | 350 / 248 | 371 / 235 | 378 / 230 |
| 220 | 312 | 313 / 239 | 325 / 233 | 329 / 230 | 356 / 257 | 376 / 244 | 384 / 239 |
| 230 | 318 | 320 / 248 | 332 / 242 | 336 / 240 | 361 / 265 | 381 / 254 | 389 / 248 |
| 240 | 325 | 327 / 258 | 340 / 252 | 344 / 249 | 367 / 274 | 387 / 263 | 394 / 257 |
| 250 | 332 | 335 / 267 | 347 / 261 | 351 / 259 | 373 / 283 | 392 / 272 | 400 / 267 |
| 260 | 339 | 342 / 276 | 354 / 271 | 358 / 268 | 379 / 292 | 398 / 281 | 405 / 276 |
| 270 | 346 | 350 / 286 | 362 / 280 | 365 / 278 | 385 / 300 | 403 / 290 | 411 / 285 |
| 280 | 353 | 357 / 295 | 369 / 290 | 373 / 288 | 391 / 309 | 409 / 299 | 416 / 294 |
| 290 | 360 | 365 / 304 | 376 / 299 | 380 / 297 | 397 / 318 | 415 / 308 | 422 / 304 |
| 300 | 368 | 372 / 314 | 384 / 309 | 387 / 307 | 403 / 327 | 420 / 318 | 427 / 313 |
| 310 | 375 | 380 / 323 | 391 / 319 | 395 / 317 | 409 / 336 | 426 / 327 | 433 / 323 |
| 320 | 382 | 387 / 333 | 398 / 328 | 402 / 326 | 416 / 345 | 432 / 336 | 439 / 332 |
| 330 | 390 | 395 / 342 | 406 / 338 | 409 / 336 | 422 / 354 | 438 / 345 | 444 / 341 |
| 340 | 397 | 403 / 352 | 413 / 347 | 417 / 346 | 428 / 363 | 444 / 355 | 450 / 351 |
| 350 | 405 | 410 / 361 | 420 / 357 | 424 / 355 | 435 / 372 | 450 / 364 | 456 / 360 |

**No number of shields keeps the plate below 150 °C once the box is at 150 °C or more.**
The plate always runs hotter than the air around it, and shields only cut its rise above
the box. At a 150 °C box the plate's peak rise is +116 °C with no shield, then +25, +18 and
+15 °C with one, two and three ε 0.2 shields, or +47, +33 and +26 °C at ε 0.45. At 350 °C
the rises shrink to +55 °C unshielded and +11 to +5 °C shielded. Beyond one shield the
gain is small, because heat also crosses the gaps by conduction through the air. To keep
the plate below 150 °C, the box must be cooled or vented to below about 125 °C. The
preliminary table in issue #4 gives 147 °C at a 125 °C box with two ε 0.2 shields.

The hottest shield, the one nearest the tube, runs well above the box. With ε 0.2 it sits
+110 to +126 °C above a 150 °C box and +60 to +74 °C above a 350 °C box. With ε 0.45 the
figures are +166 to +196 °C and +85 to +106 °C. Across the sweep that is 260–456 °C, well
inside the scaling limits of annealed 316 and 310 stainless. Its peak rises slightly as
shields are added, because the shields behind it block its view of the cooler box walls.

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
