# Guenthner inversion runs (March 2026)

Runs conducted by Willy Guenthner in late March 2026: he committed the script
and mrm.csv on March 30, 2026 (bea3d63) and the three output folders on
April 2, 2026 (558a073). The script here is the exact committed state, which
corresponds to the "all constraints" run (`name = "MRM_jacobsville"` with the
constrained block active). The other two runs used the same script with the
constraint block comment-toggled and the base name edited:

- `unconstrained_output/` — `constraint = Constraint()`, `name = "MRM"`
  (output prefix `MRM_`)
- `no_jacobsville_box_output/` — MCR burial + Cambrian constraints
  (output prefix `MRM_constrained_`)
- `jacobsville_box_output/` — all three constraints
  (output prefix `MRM_jacobsville_constrained_`)

Constraint priors in these runs used **Normal age distributions**:
`Normal(505, 20)`, `Normal(990, 10)`, `Normal(1075, 25)` Ma, with Uniform
temperature priors (20–50, 20–50, 250–400 °C).

No RNG seed was set, so the runs are not exactly
reproducible.
