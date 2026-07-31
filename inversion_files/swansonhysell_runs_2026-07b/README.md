# Inversion runs, July 2026 — batch b (seed 12345)

Thermochron.jl time-temperature inversions of the combined zircon and apatite
(U-Th)/He dataset from crystalline basement of the Montreal River monocline
(`../mrm.csv`: 54 ZHe and 13 AHe analyses, 67 chronometers total). Runs
conducted by N. Swanson-Hysell, launched July 12, 2026.

This batch is a **seed-sensitivity test**. It repeats the two configurations of
batch a (`../swansonhysell_runs_2026-07a/`) with an independent random seed, and
adds the `all_constraints` configuration. Because the batch a unconstrained chain
recovered a t-T history distinct from two earlier unconstrained runs at
comparable likelihood, the posterior of this deep-time inversion appears to be
multimodal. Running independent chains of the same configuration establishes
which features are properties of the data and which are artifacts of an
individual chain's incomplete exploration of the posterior.

## Contents

- `Project.toml`, `Manifest.toml` — pinned Julia environment (Thermochron.jl
  v0.16.13, Julia v1.12.5, and the exact version of every dependency), identical
  to batch a
- `unconstrained.log`, `no_jacobsville.log`, `all_constraints.log` — run logs,
  including the configuration/seed header, package versions, post-burn-in
  summary statistics, and runtime (progress-bar output stripped)
- `MRM_*_seed12345_*` — posterior t-T path density maps (`_tT.pdf`), observed
  versus predicted date-eU plots for each mineral (`_Age-eU.pdf`), kinetic
  parameter posteriors (`_kinetics.pdf`), and MCMC diagnostic traces of
  log-likelihood, acceptance rate, and t-T node count (`_lldist.png`,
  `_acceptance.png`, `_ndist.png`)

The diagnostic traces are archived as 150 dpi rasterizations because the vector
originals plot all 10⁶ steps and run to ~15 MB each. The inversion script is the
one at `../tTinversion_mrm.jl`, identical to the copy archived in batch a except
for the output directory. The full posterior distributions (`.jld2`, ~800 MB
each) are not archived, as each chain can be regenerated exactly from its seed
and the pinned environment.

## Reproducing these runs

From the `inversion_files/` directory:

```
julia --project=swansonhysell_runs_2026-07b tTinversion_mrm.jl unconstrained 12345
julia --project=swansonhysell_runs_2026-07b tTinversion_mrm.jl no_jacobsville 12345
julia --project=swansonhysell_runs_2026-07b tTinversion_mrm.jl all_constraints 12345
```

The script takes the constraint configuration as its first argument and the
random seed as its second. Each run takes roughly 20 to 30 hours; the three were
run concurrently.

## Configuration

All runs used 1,000,000 MCMC steps (500,000 burn-in), an 8 Ma timestep over a
2704 Ma model span, 15 to 50 t-T nodes, and the `MCMC_varkinetics` routine, in
which RDAAM (apatite) and ZRDAAM (zircon) damage-diffusivity parameters are
resampled during the inversion. Geological constraints are applied as uniform
(box) priors:

| Configuration | Constraints applied |
| --- | --- |
| `unconstrained` | none |
| `no_jacobsville` | Midcontinent Rift burial (250–400 °C, 1050–1100 Ma); Cambrian near-surface (20–50 °C, 485–505 Ma) |
| `all_constraints` | the above, plus Jacobsville near-surface (20–50 °C, 980–1000 Ma) |

Comparison of `no_jacobsville` and `all_constraints` tests the extent to which
imposing near-surface conditions at the time of Jacobsville Formation deposition
determines the inferred Neoproterozoic thermal history, as opposed to that
history being required by the (U-Th)/He data themselves.

## Results

All three runs completed July 13, 2026.

| Run | Seed | Mean log-likelihood | Mean acceptance | Mean nodes |
| --- | --- | --- | --- | --- |
| `unconstrained` | 12345 | −492.9 | 0.43 | 19.8 |
| `no_jacobsville` | 12345 | −471.3 | 0.32 | 21.2 |
| `all_constraints` | 12345 | −482.7 | 0.33 | 19.0 |

Statistics are means over the post-burn-in portion of each chain.

### Cross-seed comparison with batch a

**Unconstrained: the two seeded chains agree with each other.** Seed 12345
(−492.9) and batch a's seed 20260710 (−492.3) recover closely similar t-T path
distributions — warm Proterozoic residence, an excursion to ~250–300 °C
culminating ca. 1050–1000 Ma, rapid cooling to near-surface conditions after
ca. 1000 Ma, and moderate Neoproterozoic–Phanerozoic temperatures. Given a fixed
software environment, independent seeds therefore reproduce.

**However, both seeded chains sit at appreciably lower likelihood than the two
earlier unseeded runs** (April 2026: −477; Guenthner March 2026: −470), which
each found different histories again. A gap of 15 to 22 log-likelihood units is
not "comparable" — the earlier chains fit the data *better* than either seeded
chain. The unconstrained posterior therefore has at least three distinct modes
across four chains, and 10⁶ steps is evidently not sufficient for a chain to
locate the highest-likelihood mode. Neither seeded chain should be treated as a
converged characterization of the unconstrained posterior.

Note that the Guenthner runs used Thermochron.jl v0.15.xx rather than the
v0.16.13 used here (see `../guenthner_runs_2026-03/`), so that chain is not
strictly comparable. The April 2026 run, however, used v0.16.13 — the same
version as these — so the mode difference between it and the seeded chains is
not attributable to package version.

W. Guenthner's assessment (July 2026) is that the unconstrained prior volume is
simply very large (0 to 400 °C over 2700 to 0 Ma), so a chain has an enormous
amount of t-T space to search and different chains landing in different regions
is expected behavior rather than a defect; he anticipates that chains of order
10⁷ iterations would be needed for convergence and agreement. He emphasizes that
the pronounced cooling episode ca. 1000 Ma emerges from every chain regardless,
indicating that it is a strong feature of the data rather than of any individual
inversion. This is consistent with the cross-chain-robust features noted below.

A practical constraint on chain length: the recorded posterior arrays scale
linearly with step count and total ~1.3 GB at 10⁶ steps, so a 10⁷-step chain
would require ~13 GB of the 16 GB available on the machine used here (before
runtime and save-time overhead). Chains beyond ~5 × 10⁶ steps will need a
larger-memory machine.

**Rift burial + Cambrian: the seeds do not agree.** Seed 12345 (−471.3) and
batch a's seed 20260710 (−481.5) differ by 10 log-likelihood units and place
different amounts of Neoproterozoic path density: seed 12345 keeps paths at
moderate temperatures (~50–150 °C) through the Neoproterozoic, while seed
20260710 cools to near-surface ca. 950 Ma and then reheats to ~200 °C ca. 700 to
600 Ma. The timing of Neoproterozoic cooling is thus seed-dependent even with
the geological constraints imposed.

Features that are common to all chains: entry into the Midcontinent Rift burial
box with cooling culminating ca. 1050–1000 Ma, and post-Jacobsville residence
above ~80 to 100 °C at some point prior to final cooling.

A second `all_constraints` chain (seed 20260710) is running in
`../swansonhysell_runs_2026-07a/` to complete a three-configuration by two-seed
matrix.
