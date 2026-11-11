# bie: boundary-integral solver

Solves the mark of [MODEL.md](../MODEL.md) with boundary integral equations. Needs numpy and scipy; `selftest.py` also needs mpmath. It shares no code with `fem/`: it reads the material files itself and builds the profile itself.

```python
from bie.solver import solve
sol = solve(marks, channels)     # marks: 12 numbers, or one row of 12 per mark; channels: [(wavelength_nm, "TE" or "TM"), ...]
sol["Ap"], sol["Am"]             # complex amplitudes of the reflected +1 and -1 orders in air
sol["ep"], sol["em"]             # their diffraction efficiencies
```

Any wavelength from 400 to 1100 nm works. One mark in one channel takes 5 to 15 s on a laptop.

## Method

- **Formulation.** Only the silicon/oxide profile of one period is discretised. The oxide surface is built into the Green's function of the oxide (the periodic free-space function plus the smooth part reflected by the plane `z = h + t`), so the air and the silicon half-space are exact. The unknowns are the field along the lines and its weighted normal derivative on the profile, in a Müller-type system of the second kind; the reflected orders follow from the representation formula.
- **Periodic Green's functions.** Ewald splitting with a frequency-adapted parameter, for the complex wavenumber of silicon as well.
- **Discretisation.** Nyström on composite 16-point Gauss-Legendre panels, graded towards every breakpoint; logarithmic product weights on a panel's own interactions, dyadically graded rules on its neighbours'.
- **Corners.** Recursively compressed inverse preconditioning (Helsing) at every breakpoint, 55 dyadic levels by default, so the field singularity at a silicon corner in TM costs no accuracy.

Wavelengths within about 1e-8 nm of an oxide Rayleigh wavelength (`2000 n_oxide / m`: 419.47, 487.68, 583.44, 727.35, 967.23 nm), where the periodic Green's function is singular, are refused.

## How accurate

- `python bie/selftest.py` (half a minute): the Ewald sums against brute-force sums of 60001 orders, the quadrature weights against mpmath, and a flat interface against the thin-film formula.
- Every one of the 51 pairs was solved at two resolutions (the default, and panels of at most 40 nm with 70 corner levels). Each record in `records/` gives the largest change; over all pairs and all eight channels it stays below 1e-11 nm in the reading.
- Against `fem/`, on every pair: within 7e-9 nm in the reading, which is the finite elements' own resolution (`python summary.py`).
