# Inversion runs, July 2026 — batch a (seed 20260710)

> **Superseded.** These runs used a version of `../mrm.csv` in which the
> `raw He age [Ma]` column held alpha-ejection-corrected dates rather than raw
> dates. Because Thermochron.jl applies its own alpha-ejection correction
> internally, the input dates were too old by a grain-size-dependent factor
> (median 1.35 for the zircons, 1.53 for the apatites). The error was found in
> August 2026 and is described in `../swansonhysell_runs_2026-08a/README.md`,
> which supersedes this batch. The t-T histories and the convergence behaviour
> discussed below should not be interpreted; the seed-dependence documented
> here does not persist once the input data are corrected.

Thermochron.jl time-temperature inversions of the combined zircon and apatite
(U-Th)/He dataset from crystalline basement of the Montreal River monocline
(`../mrm.csv`: 54 ZHe and 13 AHe analyses, 67 chronometers total). Runs
conducted by N. Swanson-Hysell, completed July 11, 2026.

These are the first runs in this project made with a **fixed random seed and a
pinned software environment**, and are therefore exactly reproducible. Earlier
runs (`../guenthner_runs_2026-03/`, `../swansonhysell_runs_2026-04/`) set no
seed and are not.

## Contents

- `tTinversion_mrm.jl` — the exact inversion script used for these runs
- `Project.toml`, `Manifest.toml` — pinned Julia environment (Thermochron.jl
  v0.16.13, Julia v1.12.5, and the exact version of every dependency)
- `unconstrained.log`, `no_jacobsville.log` — run logs, including the
  configuration/seed header, package versions, post-burn-in summary statistics,
  and runtime (progress-bar output stripped)
- `MRM_unconstrained_*`, `MRM_no_jacobsville_*` — posterior t-T path density
  maps (`_tT.pdf`), observed versus predicted date-eU plots for each mineral
  (`_Age-eU.pdf`), kinetic parameter posteriors (`_kinetics.pdf`), and MCMC
  diagnostic traces of log-likelihood, acceptance rate, and t-T node count
  (`_lldist.png`, `_acceptance.png`, `_ndist.png`)

The diagnostic traces are archived as 150 dpi rasterizations because the vector
originals plot all 10⁶ steps and run to ~15 MB each. The full posterior
distributions (`.jld2`, ~800 MB each) are not archived, as each chain can be
regenerated exactly from its seed and the pinned environment.

## Reproducing these runs

From the `inversion_files/` directory:

```
julia --project=swansonhysell_runs_2026-07a tTinversion_mrm.jl unconstrained
julia --project=swansonhysell_runs_2026-07a tTinversion_mrm.jl no_jacobsville
```

The script takes the constraint configuration as its first argument
(`unconstrained`, `no_jacobsville`, or `all_constraints`) and an optional random
seed as its second; with no seed argument it uses the default seed 20260710 that
produced the results here. Each run takes roughly 20 hours.

## Configuration

Both runs used 1,000,000 MCMC steps (500,000 burn-in), an 8 Ma timestep over a
2704 Ma model span, 15 to 50 t-T nodes, and the `MCMC_varkinetics` routine, in
which RDAAM (apatite) and ZRDAAM (zircon) damage-diffusivity parameters are
resampled during the inversion. Geological constraints are applied as uniform
(box) priors:

| Configuration | Constraints applied |
| --- | --- |
| `unconstrained` | none |
| `no_jacobsville` | Midcontinent Rift burial (250–400 °C, 1050–1100 Ma); Cambrian near-surface (20–50 °C, 485–505 Ma) |

## Results

| Run | Seed | Mean log-likelihood | Mean acceptance | Mean nodes |
| --- | --- | --- | --- | --- |
| `unconstrained` | 20260710 | −492.3 | 0.41 | 20.1 |
| `no_jacobsville` | 20260710 | −481.5 | 0.33 | 19.7 |

Statistics are means over the post-burn-in portion of each chain.

**The unconstrained posterior is multimodal.** This chain recovered a t-T
history distinct from both the April 2026 run (mean log-likelihood −477) and the
Guenthner March 2026 run (−470), despite the April run using the same package
version, data, and physics and differing only in chain stochasticity. Comparable
likelihoods across visibly different histories indicate that a single 10⁶-step
chain does not fully explore the posterior of this deep-time inversion.
Interpretation should accordingly be focused on features that reproduce
across independent chains. A second batch with an independent seed, and adding
the `all_constraints` configuration, is in `../swansonhysell_runs_2026-07b/`.

Features that reproduce across all chains run to date include entry into
Midcontinent Rift burial temperatures culminating ca. 1050–1000 Ma followed by
pronounced cooling, and post-Jacobsville residence above ~80–100 °C at some
point prior to final cooling. The timing of Neoproterozoic cooling does **not**
reproduce across chains.
