# Swanson-Hysell inversion runs (April 9, 2026)

Runs conducted by NSH on April 9, 2026 with Thermochron.jl v0.16.13 on Julia
1.12. The script here is the comment-toggle version used for these runs,
archived in the state with Option 1 (unconstrained) active; the
`no_jacobsville` and `all_constraints` runs were made by toggling the
constraint blocks.

Key difference from the Guenthner March 2026 runs (`../guenthner_runs_2026-03/`):
constraint age priors were changed from Normal distributions to **Uniform (box)
constraints** — `Uniform(485, 505)`, `Uniform(980, 1000)`, `Uniform(1050, 1100)`
Ma — a significant change for how the constraints condition the inversion.
Temperature priors are unchanged.

No RNG seed was set, so these runs are not exactly reproducible. Note that the
unconstrained run here recovered a substantially different t-T history than
Guenthner's unconstrained run at comparable log-likelihood (mean ca. −477 vs
−470), with a lower acceptance rate (ca. 0.30 vs 0.47) — consistent with a
multimodal posterior sampled by chains that did not mix between modes.

Outputs: `MRM_unconstrained_*`, `MRM_no_jacobsville_*`, `MRM_all_constraints_*`
(compressed posterior distributions in the `.jld2` files, diagnostic and result
plots in the PDFs).
