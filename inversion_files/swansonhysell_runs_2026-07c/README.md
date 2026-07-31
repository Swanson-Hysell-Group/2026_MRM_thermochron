# Inversion runs, July 2026 — batch c (additional seeds)

Additional Markov chains run to characterize chain-to-chain variability of the
Thermochron.jl time-temperature inversions, for the seed-comparison figure in the
supporting information. Same data (`../mrm.csv`), same script
(`../tTinversion_mrm.jl`), and the same pinned Julia environment (Thermochron.jl
v0.16.13, Julia v1.12.5) as batches a and b — these chains differ only in their
random seed.

Background and the two-seed results that motivated these runs are in
`../swansonhysell_runs_2026-07a/README.md` and
`../swansonhysell_runs_2026-07b/README.md`.

## Why additional seeds

The unconstrained prior volume is large (0 to 400 °C across 2700 to 0 Ma), so a
single transdimensional chain of 10⁶ steps explores only part of the accessible
t-T space, and independent chains can settle in different regions at different
likelihoods. Replicate chains at independent seeds make that variability
explicit, and identify which features of the posterior are properties of the data
rather than of an individual chain.

Seeds are concentrated on the two configurations shown in the main text
(`unconstrained` and `no_jacobsville`), since it is their reproducibility that
bears on the interpretation.

## Contents

- `Project.toml`, `Manifest.toml` — pinned Julia environment (identical to
  batches a and b)
- `*_seed<N>.log` — run logs with configuration/seed header, package versions,
  post-burn-in summary statistics, and runtime
- `MRM_*_seed<N>_*` — posterior t-T path density maps (`_tT.pdf`), observed
  versus predicted date-eU plots (`_Age-eU.pdf`), kinetic parameter posteriors
  (`_kinetics.pdf`), and MCMC diagnostic traces (`_lldist.png`,
  `_acceptance.png`, `_ndist.png`; rasterized, see batch a README)

Full posterior distributions (`.jld2`) are not archived; each chain is exactly
reproducible from its seed and the pinned environment.

## Reproducing these runs

From the `inversion_files/` directory, passing configuration, seed, and output
directory:

```
julia --project=swansonhysell_runs_2026-07c tTinversion_mrm.jl unconstrained 777 swansonhysell_runs_2026-07c
julia --project=swansonhysell_runs_2026-07c tTinversion_mrm.jl no_jacobsville 777 swansonhysell_runs_2026-07c
```

## Runs

| Configuration | Seed | Mean log-likelihood | Mean acceptance | Status |
| --- | --- | --- | --- | --- |
| `unconstrained` | 777 | | | launched 2026-07-14 |
| `no_jacobsville` | 777 | | | launched 2026-07-14 |

Results will be recorded here on completion and combined with batches a and b for
the supporting information seed-comparison figure.
