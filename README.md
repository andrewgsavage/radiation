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
is held at 650 °C. Up to ten 0.1 mm stainless shields are spaced evenly in the 10 mm gap,
at ε 0.2 (bright) or 0.45 (oxidized).

**Shields ε 0.2** (hottest shield / plate peak, °C; columns are the number of shields)

| Box | 0 | 1 | 2 | 3 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|
| 150 | – / 266 | 292 / 177 | 315 / 170 | 324 / 167 | 328 / 166 | 331 / 164 | 331 / 163 | 330 / 162 |
| 160 | – / 272 | 299 / 186 | 322 / 179 | 331 / 177 | 335 / 175 | 338 / 173 | 338 / 172 | 337 / 171 |
| 170 | – / 279 | 306 / 195 | 329 / 189 | 338 / 186 | 342 / 184 | 345 / 182 | 345 / 181 | 344 / 181 |
| 180 | – / 285 | 313 / 204 | 336 / 198 | 344 / 195 | 348 / 193 | 351 / 192 | 352 / 191 | 351 / 190 |
| 190 | – / 292 | 320 / 213 | 342 / 207 | 351 / 204 | 355 / 203 | 358 / 201 | 358 / 200 | 358 / 199 |
| 200 | – / 298 | 327 / 223 | 349 / 216 | 357 / 214 | 361 / 212 | 365 / 210 | 365 / 209 | 364 / 209 |
| 210 | – / 305 | 333 / 232 | 355 / 226 | 363 / 223 | 368 / 222 | 371 / 220 | 371 / 219 | 370 / 218 |
| 220 | – / 312 | 340 / 241 | 361 / 235 | 370 / 232 | 374 / 231 | 377 / 229 | 377 / 228 | 377 / 228 |
| 230 | – / 318 | 346 / 250 | 367 / 244 | 376 / 242 | 380 / 241 | 383 / 239 | 384 / 238 | 383 / 237 |
| 240 | – / 325 | 353 / 260 | 374 / 254 | 382 / 251 | 386 / 250 | 389 / 249 | 390 / 248 | 389 / 247 |
| 250 | – / 332 | 360 / 269 | 380 / 263 | 388 / 261 | 392 / 260 | 395 / 258 | 396 / 257 | 395 / 257 |
| 260 | – / 339 | 366 / 278 | 386 / 273 | 394 / 270 | 398 / 269 | 401 / 268 | 402 / 267 | 402 / 266 |
| 270 | – / 346 | 373 / 287 | 392 / 282 | 400 / 280 | 404 / 279 | 407 / 277 | 408 / 277 | 408 / 276 |
| 280 | – / 353 | 379 / 297 | 399 / 292 | 406 / 290 | 410 / 288 | 414 / 287 | 414 / 286 | 414 / 286 |
| 290 | – / 360 | 386 / 306 | 405 / 301 | 413 / 299 | 416 / 298 | 420 / 297 | 420 / 296 | 420 / 295 |
| 300 | – / 368 | 393 / 315 | 411 / 311 | 419 / 309 | 423 / 308 | 426 / 306 | 427 / 306 | 426 / 305 |
| 310 | – / 375 | 400 / 325 | 418 / 320 | 425 / 318 | 429 / 317 | 432 / 316 | 433 / 315 | 432 / 315 |
| 320 | – / 382 | 406 / 334 | 424 / 330 | 431 / 328 | 435 / 327 | 438 / 326 | 439 / 325 | 439 / 325 |
| 330 | – / 390 | 413 / 344 | 430 / 339 | 437 / 337 | 441 / 336 | 444 / 335 | 445 / 335 | 445 / 334 |
| 340 | – / 397 | 420 / 353 | 437 / 349 | 444 / 347 | 447 / 346 | 450 / 345 | 451 / 344 | 451 / 344 |
| 350 | – / 405 | 427 / 362 | 443 / 358 | 450 / 357 | 453 / 356 | 457 / 355 | 457 / 354 | 457 / 354 |

**Shields ε 0.45**

| Box | 0 | 1 | 2 | 3 | 4 | 6 | 8 | 10 |
|---|---|---|---|---|---|---|---|---|
| 150 | – / 266 | 350 / 201 | 384 / 187 | 399 / 181 | 406 / 177 | 414 / 172 | 417 / 170 | 418 / 169 |
| 160 | – / 272 | 355 / 209 | 388 / 196 | 403 / 189 | 411 / 186 | 418 / 181 | 421 / 179 | 422 / 177 |
| 170 | – / 279 | 360 / 217 | 393 / 204 | 407 / 198 | 415 / 194 | 422 / 190 | 425 / 188 | 426 / 186 |
| 180 | – / 285 | 364 / 226 | 397 / 213 | 411 / 207 | 419 / 203 | 426 / 199 | 429 / 197 | 430 / 195 |
| 190 | – / 292 | 369 / 234 | 401 / 222 | 415 / 216 | 423 / 212 | 431 / 208 | 434 / 206 | 435 / 205 |
| 200 | – / 298 | 374 / 243 | 406 / 230 | 420 / 224 | 427 / 221 | 435 / 217 | 438 / 215 | 439 / 214 |
| 210 | – / 305 | 379 / 251 | 410 / 239 | 424 / 233 | 431 / 230 | 439 / 226 | 442 / 224 | 443 / 223 |
| 220 | – / 312 | 383 / 260 | 414 / 248 | 428 / 243 | 435 / 239 | 442 / 236 | 445 / 234 | 447 / 232 |
| 230 | – / 318 | 388 / 268 | 418 / 257 | 432 / 252 | 439 / 249 | 446 / 245 | 449 / 243 | 450 / 242 |
| 240 | – / 325 | 393 / 277 | 423 / 266 | 436 / 261 | 443 / 258 | 450 / 254 | 453 / 252 | 454 / 251 |
| 250 | – / 332 | 398 / 286 | 427 / 275 | 440 / 270 | 447 / 267 | 454 / 264 | 457 / 262 | 458 / 261 |
| 260 | – / 339 | 403 / 294 | 431 / 284 | 444 / 279 | 451 / 276 | 458 / 273 | 461 / 271 | 462 / 270 |
| 270 | – / 346 | 408 / 303 | 436 / 293 | 448 / 288 | 455 / 286 | 462 / 282 | 465 / 281 | 467 / 280 |
| 280 | – / 353 | 413 / 312 | 440 / 302 | 453 / 297 | 459 / 295 | 466 / 292 | 469 / 290 | 471 / 289 |
| 290 | – / 360 | 418 / 320 | 445 / 311 | 457 / 307 | 464 / 304 | 471 / 301 | 474 / 300 | 475 / 299 |
| 300 | – / 368 | 424 / 329 | 449 / 320 | 461 / 316 | 468 / 313 | 475 / 311 | 478 / 309 | 479 / 308 |
| 310 | – / 375 | 429 / 338 | 454 / 329 | 466 / 325 | 472 / 323 | 479 / 320 | 482 / 319 | 483 / 318 |
| 320 | – / 382 | 434 / 347 | 459 / 338 | 470 / 335 | 477 / 332 | 483 / 330 | 486 / 328 | 487 / 327 |
| 330 | – / 390 | 440 / 356 | 464 / 348 | 475 / 344 | 481 / 342 | 488 / 339 | 490 / 338 | 492 / 337 |
| 340 | – / 397 | 445 / 365 | 469 / 357 | 479 / 353 | 486 / 351 | 492 / 349 | 495 / 347 | 496 / 346 |
| 350 | – / 405 | 451 / 373 | 474 / 366 | 484 / 363 | 490 / 361 | 496 / 358 | 499 / 357 | 500 / 356 |

**No number of shields keeps the plate below 150 °C once the box is at 150 °C or more.**
The plate always runs hotter than the air around it, and shields only cut its rise above
the box. In a 150 °C box the plate peaks at 266 °C unshielded and 177 °C behind one ε 0.2
shield. Ten shields only bring it to 162 °C. At ε 0.45 the plate peaks at 201 °C behind one
shield and 169 °C behind ten. More shields stop helping after three or four: they block
radiation, but conduction through the air across the 10 mm gap stays about the same
however finely the gap is split. To keep the plate below 150 °C, the box must be cooled
or vented to below about 125 °C. The preliminary table in issue #4 gives 147 °C at a
125 °C box with two ε 0.2 shields.

The hottest shield, the one nearest the tube, runs 290–460 °C with ε 0.2 and 350–500 °C
with ε 0.45 across the sweep. That is well inside the scaling limits of annealed 316 and
310 stainless. Its peak rises as shields are added, because the shields behind it are
hot and block its view of the cooler box walls.

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
