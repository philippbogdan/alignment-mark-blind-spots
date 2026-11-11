# fem: finite-element solver

Solves the mark of [MODEL.md](../MODEL.md) with high-order finite elements. Needs [NGSolve](https://ngsolve.org) (`pip install ngsolve`), numpy and scipy. It shares no code with `bie/`.

```python
from fem.solver import amplitudes
sol = amplitudes(marks, channels, order=8)    # or order=10; same keys as bie.solver.solve
```

- **Geometry.** The exact profile, straight edges and circular arcs, meshed with curved elements on a cell periodic in x, with geometric (hp) refinement at every silicon corner.
- **Radiation.** Exact Dirichlet-to-Neumann conditions in Fourier modes on a plane 100 nm above the line and on a plane in the silicon.
- **Oxide surface.** The solver computes the reflection matrix of the grating under an oxide half-space; the oxide/air surface is then added by the exact multiple-reflection formula, so the oxide thickness enters analytically.

Two resolutions are used in `records/`: order 8 (a minute per pair on a laptop, about 1.5 GB) and order 10 (a few minutes, about 6 GB).

`geometry.py` also holds the admissibility rules of the model, which `check.py` applies to every pair.
