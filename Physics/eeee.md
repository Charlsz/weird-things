blackhole_simulator.py
~~~~~~~~~~~~~~~~~~~~~~

A minimalist “black‑hole simulator” that solves the Schwarzschild
geodesic equations for massive particles and photons.

Features
--------
* Exact Schwarzschild metric (G, c kept explicit – SI units).
* Closed‑form expressions for the conserved energy (E) and angular momentum (L)
  of a circular orbit.
* First‑integral radial equation derived from the normalization of the
  4‑velocity (timelike) or 4‑momentum (null).
* ODE integration in proper time τ with SciPy’s solve_ivp.
* Simple polar → Cartesian visualisation, horizon plotted as a filled circle.

Author   : ChatGPT (OpenAI) – 2024‑06
License  : MIT (feel free to adapt)