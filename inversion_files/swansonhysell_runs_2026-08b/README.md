# Inversion runs, August 2026 — batch b (grain-exclusion sensitivity test)

Thermochron.jl time-temperature inversions run on the batch 08a dataset with
two zircon analyses removed, to test how much of the recovered t-T history
depends on them. Input file `../mrm_no_MRM4_z3_z5.csv`: 53 ZHe and 13 AHe
analyses, 66 chronometers, identical to `../mrm.csv` apart from the two
excluded rows.

Script, seeds, constraint configurations, and the pinned Julia environment are
the same as `../swansonhysell_runs_2026-08a/`, so each run here pairs with an
08a run that differs **only** in those two grains.

## The two grains

`TAI_MRM4_z3` (eU 903, raw date 746 Ma) and `TAI_MRM4_z5` (eU 721, raw date
1288 Ma) sit far above the Archean basement date-eU array. Grains of comparable
eU give 89 to 313 Ma, so they are 3 to 14 times older than their neighbours.
The anomaly is in the helium rather than the parent: their He concentrations
(3880 and 5654 nmol/g) are the third-highest and highest in the Archean array
against a median of 1057 nmol/g, and high He concentration is otherwise a
property of the high-eU grains. Within sample MRM4, `z1` has essentially the
same eU as `z5` (622 versus 721) but a third of the helium and half the date,
so a thermal-history or radiation-damage explanation is not available — every
grain in a sample shares one t-T path and eU is the damage proxy.

Checks against the laboratory files in the HAL repository
(https://github.com/wrguenthner/HAL_data, `projects_2023/TA23_MRM_MCR/`) rule
out several artifacts: reported U and Th agree with the raw ICP-MS intensities
at the session-median calibration for both grains, He re-extracts are normal
(3.4 and 3.5 per cent), neither grain is affected by the TA_/TAI_ prefix
collision (sample MRM4 has only TAI_ analyses), and grain dimensions and F_T
are unremarkable. Zonation cannot account for it either, since a low-eU core
lowers rather than raises bulk He concentration. What remains is either parent
that was never dissolved and so never measured — an incompletely dissolved
U-Th-rich inclusion whose helium was still released on degassing — or helium
implanted from a U-Th-rich neighbouring phase. Th/U for the two grains (0.63
and 0.65) is two to four times that of the other MRM4 analyses, which is
consistent with a Th-rich accessory phase being involved, though Th/U near 0.6
is not extreme across the dataset as a whole.

These grains violate no criterion applied to the dataset so far — unlike the
MRM34 analyses excluded in August 2026, their dates do not exceed the
crystallization age — and they are retained in `../mrm.csv` and in batch 08a.
This batch tests their influence rather than presuming it.

## Why it matters beyond the two dates

`tTinversion_mrm.jl` calls `empiricaluncertainty!` to estimate the model
uncertainty (`σcalc`) for each mineral from the excess dispersion of the data
themselves. Two grains lying far off the date-eU array therefore inflate the
uncertainty assigned to *every* zircon, not only their own, and so loosen the
constraint the whole array places on the inversion. The paired comparison
measures that effect as well as the direct pull of the two dates.

## Contents

- `Project.toml`, `Manifest.toml` — pinned Julia environment (Thermochron.jl
  v0.16.13, Julia v1.12.5), identical to batches 07a/b/c and 08a
- `run_batch.sh` — batch runner; refuses to start if the inversion script does
  not accept an input-file argument, if the reduced input file is missing, or
  if another batch is still running
- `<config>_seed<N>.log` — run logs, with the configuration, seed, input file
  and analysis count in the header
- `MRM_*_tT_density.png`, `MRM_*_tT_density.json` — the posterior path-density
  image of each chain exported from its `.jld2` posterior by
  `../export_tT_density.jl` (8-bit, 1200 × 1800, normalised exactly as in the
  `_tT.pdf` plot; ranges and chain statistics in the `.json`). These are the
  inputs to `../../code/make_inversion_figures.py`, which assembles the
  composite manuscript figures in `../../figures/`
- `MRM_<config>_seed<N>_no_MRM4_z3_z5_*` — posterior t-T path density maps,
  date-eU plots, kinetic parameter posteriors, and MCMC diagnostic traces

## Running

From the `inversion_files/` directory, after batch 08a has finished:

```
nohup bash swansonhysell_runs_2026-08b/run_batch.sh > swansonhysell_runs_2026-08b/batch.log 2>&1 &
```

This requires the fourth command line argument (input data file) added to
`../tTinversion_mrm.jl`; the runner checks for it and stops if it is absent.
The default is a single `no_jacobsville` chain at seed 20260710, the
configuration carried in the main text, paired against
`../swansonhysell_runs_2026-08a/no_jacobsville_seed20260710.log`. Add entries
to `CONFIGS` and `SEEDS` in the runner to widen the test; each addition needs
an 08a counterpart at the same seed for the pairing to hold.

## Results

| Configuration | Seed | Mean log-likelihood | Mean acceptance | Mean nodes | Status |
| --- | --- | --- | --- | --- | --- |
| `no_jacobsville` | 20260710 | -423.1 | 0.30 | 20.1 | completed 2026-08-29 (18.1 h) |

Log-likelihoods are not comparable to the 08a runs, since the number of
chronometers being fit differs.
