"""Forward-model test of the burial-temperature bound from pre-depositional Freda ZHe dates.

Detrital zircon grains in the Freda Formation retain (U-Th)/He dates that
pre-date deposition (ca. 1080-1050 Ma), which requires that burial heating
after deposition never fully reset the ZHe system in these grains. This
script uses pythermo (https://github.com/OpenThermochronology/PyThermo) with
the zircon radiation damage and annealing model of Guenthner et al. (2013)
to ask: at what sustained Neoproterozoic burial temperature would the
pre-depositional He in these grains have been erased?

Model setup:
    - Grains: the three clearly pre-depositional grains from the Freda
      samples, with measured U, Th, and radius (from data/UTh_He_zircon_data.csv):
          TAI_LRA_67_z1  1159 +/- 20 Ma   U 141, Th 48 ppm,  radius 63 um (eU 152)
          TAI_BRB_z5     1192 +/- 33 Ma   U 84,  Th 58 ppm,  radius 30 um (eU 98)
          TAI_BRB_z7     1176 +/- 29 Ma   U 247, Th 162 ppm, radius 37 um (eU 285)
    - Pre-depositional history: crystallization at 1200 Ma with cooling to
      the surface by 1160 Ma. This is a deliberately simple, minimally
      damaged (and therefore maximally He-retentive) end member: a grain
      with more inherited radiation damage would be reset at lower
      temperatures, so the thresholds found here are upper limits.
    - Post-depositional history (modeled for the upper Freda position of
      LRa-67.0; the BRb grains additionally experienced syn-depositional
      burial beneath ~3.7 km of Freda that is neglected here, again in the
      conservative direction): near-surface residence through monocline
      tilting, burial heating to a peak temperature T_hold representing
      residence in the Jacobsville/Bayfield basin, Neoproterozoic cooling
      to the surface, surface exposure at the 505 Ma Sauk transgression,
      and Paleozoic burial peaking at 70 degrees C at 250 Ma (matching
      the Phanerozoic history used in the manuscript forward models).
    - T_hold is varied from 35 (no reheating) to 200 degrees C.
    - Two burial-duration scenarios are run, since the erasure threshold
      depends on how long the grains stay hot: a peak hold from 900 to
      850 Ma reached by gradual burial from 990 Ma, with cooling complete
      by ca. 690 Ma (consistent with the Cryogenian AHe date-eU pediment
      at 688 +/- 22 Ma), and a protracted hold from 950 to 650 Ma
      following rapid burial, with cooling complete by 600 Ma.

Because each detrital grain carries an unknown pre-depositional history,
the ZHe date a grain had at deposition is unknown and partial He loss
cannot be assessed from any individual date. What can be assessed is
complete loss: each grain is therefore also modeled with no inherited He
or damage (crystallization at the 1052 Ma deposition age), and the
temperature at which the two models converge marks complete loss of the
pre-depositional He signal --- a threshold that is independent of the
assumed pre-depositional history. Above that temperature no grain could
retain a date older than deposition, in conflict with the measured dates.

Because He diffusion is linear in concentration, a grain with a larger
(older) initial He reservoir retains a proportionally larger date residual
at a given temperature, and inherited damage below the percolation
threshold increases retentivity. The complete-loss temperature is
therefore also evaluated against a worst-case end member that retained He
(and accumulated damage) since ca. 1830 Ma: the temperature at which even
this grain cannot keep a date older than deposition (ca. 170 degrees C)
is the robust upper bound quoted in the manuscript.

A separate syn-depositional model applies to the BRb grains alone. Sample
BRb-1.2 sits ~150 m above the base of the Freda Formation and was
deposited ca. 1075 Ma, then progressively buried beneath the ~4 km of
Freda strata exposed above it in the monocline until the end of Oronto
Group deposition ca. 1050 Ma, followed by rapid exhumation during
contractional inversion (near-surface by ca. 1030 Ma). The retention of
pre-depositional dates through this brief but potentially deep burial
caps the peak temperature reached at ca. 1050 Ma, which translates the
observed 4 km of overburden into a maximum syn-rift geothermal gradient
(peak T = 25 C surface temperature + 4 km x gradient), or equivalently
caps any additional burial beyond the exposed stratigraphic thickness.
The same complete-loss logic applies, with the no-inheritance comparison
grain crystallizing at the 1075 Ma deposition age of BRb. Because the
Mesoproterozoic and Neoproterozoic burial episodes draw on the same He
budget (and syn-rift heating partially anneals radiation damage, changing
later retentivity), the syn-Freda sweep is run at two fixed Neoproterozoic
burial levels spanning the AHe-permitted range (80 and 140 degrees C),
which shows how the two bounds trade off against each other.

Usage:
    pip install pythermo
    python code/freda_zhe_burial_bound.py

Output:
    - printed tables of modeled present-day dates for each grain and T_hold
    - figures/freda_zhe_burial_bound.pdf (Neoproterozoic burial sweep)
    - figures/brb_syn_freda_burial.pdf (BRb syn-depositional burial sweep)
"""

from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pythermo as pyt

# figures are written to the repository-root figures/ directory regardless of
# the working directory this script is run from
FIGURES_DIR = Path(__file__).resolve().parent.parent / "figures"

# measured pre-depositional Freda grains: name, U ppm, Th ppm, radius um, observed date, 1s
GRAINS = [
    ("LRa-67.0 z1 (eU 152)", 140.91, 48.07, 62.89, 1159.1, 19.8),
    ("BRb-1.2 z5 (eU 98)", 84.38, 57.57, 29.54, 1192.4, 32.9),
    ("BRb-1.2 z7 (eU 285)", 247.39, 161.76, 37.27, 1175.7, 28.9),
]

# Neoproterozoic (Jacobsville/Bayfield) burial hold temperatures to test (degrees C)
HOLD_TEMPS = [35, 80, 100, 120, 130, 140, 150, 160, 170, 175, 180, 200]

FREDA_DEPOSITION_AGE = 1052  # Ma, maximum depositional age of Hodgin et al. (2024)
LOG2_NODES = 7  # finite-difference node control, as in MRM_forward_models.ipynb


def build_tT_path(t_hold, scenario="protracted"):
    """Build the time-temperature path for a given Neoproterozoic hold temperature.

    Args:
        t_hold: Peak burial temperature (degrees C) representing residence
            beneath Jacobsville/Bayfield basin fill.
        scenario: "peak_900_850" ramps up gradually from 990 Ma, resides
            at t_hold from 900 to 850 Ma, and cools to the surface by
            690 Ma (consistent with the Cryogenian AHe date-eU pediment);
            "protracted" holds t_hold from 950 to 650 Ma following rapid
            burial at 990 Ma, with cooling complete by 600 Ma.

    Returns:
        2D numpy array of [time (Ma), temperature (C)] rows ordered young to
        old, as required by pythermo.tT_path.
    """
    if scenario == "peak_900_850":
        burial_nodes = [
            [690, 35],      # Neoproterozoic cooling complete (AHe pediment timing)
            [850, t_hold],  # end of burial peak
            [900, t_hold],  # start of burial peak; gradual ramp from 990 Ma
        ]
    elif scenario == "protracted":
        burial_nodes = [
            [600, 35],      # Neoproterozoic cooling complete
            [650, t_hold],  # end of Jacobsville/Bayfield burial hold
            [950, t_hold],  # start of burial hold
        ]
    else:
        raise ValueError(f"unknown scenario: {scenario}")

    return np.array([
        [0, 10],       # present-day surface
        [250, 70],     # Paleozoic burial peak (manuscript forward models)
        [505, 35],     # surface at Sauk transgression
        *burial_nodes,
        [990, 35],     # near-surface at onset of Jacobsville deposition
        [1010, 35],    # near-surface following monocline tilting
        [1051, 35],    # shallow burial within the Freda Formation
        [1052, 20],    # deposition at the surface
        [1160, 20],    # source region cooled to surface
        [1200, 180],   # crystallization in source region
    ])


def model_date(u_ppm, th_ppm, radius, t_hold, scenario="protracted", inherited=True):
    """Compute an alpha-ejection-corrected ZHe date for one grain and path.

    Args:
        u_ppm: Uranium concentration (ppm).
        th_ppm: Thorium concentration (ppm).
        radius: Equivalent spherical grain radius (micrometers).
        t_hold: Peak Neoproterozoic burial temperature (degrees C).
        scenario: Burial-duration scenario passed to build_tT_path
            ("peak_900_850" or "protracted").
        inherited: If True, include the pre-depositional (source region)
            history. If False, start the grain at the 1052 Ma deposition
            age with no inherited He or radiation damage, for the
            complete-loss comparison.

    Returns:
        Modeled present-day ZHe date (Ma) using the damage-diffusivity and
        annealing kinetics of Guenthner et al. (2013).
    """
    tT = build_tT_path(t_hold, scenario)
    if not inherited:
        tT = tT[tT[:, 0] <= FREDA_DEPOSITION_AGE]
    path = pyt.tT_path(tT)
    rho_r_array, relevant_tT = path.guenthner_anneal()
    grain = pyt.zircon(radius, LOG2_NODES, relevant_tT, rho_r_array, u_ppm, th_ppm)
    return grain.zirc_date("guenthner", "guenthner")[0]


SCENARIOS = {
    "peak_900_850": "gradual burial with peak temperature 900-850 Ma",
    "protracted": "protracted hold at peak temperature, 950-650 Ma",
}

# BRb syn-depositional burial model ------------------------------------------

BRB_DEPOSITION_AGE = 1075  # Ma, deposition of the lower Freda at BRb
BRB_OVERBURDEN_KM = 4.0    # exposed Freda stratigraphic thickness above BRb
SURFACE_TEMP = 25          # degrees C
# independent upper-bound geotherm from the preservation of ca. 1700 Ma Rb-Sr
# biotite ages in the structurally highest basement rocks through ~15 km of rift
# burial (biotite closure ~400 degrees C taken as a conservative limit; Cannon 1993)
RBSR_GEOTHERM_MAX = 25     # degrees C/km
# syn-Freda peak temperatures; at 4 km of burial these correspond to
# geothermal gradients of 10 to 50 C/km
SYN_FREDA_PEAK_TEMPS = [65, 85, 105, 125, 145, 165, 185, 205, 225]
# fixed Neoproterozoic burial levels spanning the AHe-permitted range,
# showing how the Mesoproterozoic and Neoproterozoic bounds trade off
NEO_LEVELS = [80, 140]


def build_brb_syn_freda_path(t_peak, t_neo):
    """Build the BRb path with syn-Freda burial and fixed Neoproterozoic burial.

    Args:
        t_peak: Peak temperature (degrees C) at ca. 1050 Ma from progressive
            burial beneath the Freda Formation, followed by rapid exhumation
            during contractional inversion.
        t_neo: Fixed Neoproterozoic (Jacobsville/Bayfield) burial temperature
            (degrees C), held 900 to 850 Ma as in the peak_900_850 scenario.

    Returns:
        2D numpy array of [time (Ma), temperature (C)] rows ordered young to
        old, as required by pythermo.tT_path.
    """
    return np.array([
        [0, 10],        # present-day surface
        [250, 70],      # Paleozoic burial peak (manuscript forward models)
        [505, 35],      # surface at Sauk transgression
        [690, 35],      # Neoproterozoic cooling complete (AHe pediment timing)
        [850, t_neo],   # end of Jacobsville/Bayfield burial peak
        [900, t_neo],   # start of burial peak; gradual ramp from 990 Ma
        [990, 35],      # near-surface at onset of Jacobsville deposition
        [1010, 35],     # near-surface following monocline tilting
        [1030, 35],     # exhumed to near-surface by rapid inversion
        [1050, t_peak],  # maximum syn-Freda burial at end of Oronto deposition
        [1075, SURFACE_TEMP],  # deposition of BRb at the surface
        [1160, 20],     # source region cooled to surface
        [1200, 180],    # crystallization in source region
    ])


def model_brb_date(u_ppm, th_ppm, radius, t_peak, t_neo, inherited=True):
    """Compute a ZHe date for a BRb grain through the syn-Freda burial path.

    Args:
        u_ppm: Uranium concentration (ppm).
        th_ppm: Thorium concentration (ppm).
        radius: Equivalent spherical grain radius (micrometers).
        t_peak: Peak syn-Freda burial temperature (degrees C) at ca. 1050 Ma.
        t_neo: Fixed Neoproterozoic burial temperature (degrees C).
        inherited: If False, start the grain at the 1075 Ma BRb deposition
            age with no inherited He or damage, for the complete-loss
            comparison.

    Returns:
        Modeled present-day ZHe date (Ma).
    """
    tT = build_brb_syn_freda_path(t_peak, t_neo)
    if not inherited:
        tT = tT[tT[:, 0] <= BRB_DEPOSITION_AGE]
    path = pyt.tT_path(tT)
    rho_r_array, relevant_tT = path.guenthner_anneal()
    grain = pyt.zircon(radius, LOG2_NODES, relevant_tT, rho_r_array, u_ppm, th_ppm)
    return grain.zirc_date("guenthner", "guenthner")[0]


def main():
    results = {}
    residual = {}
    for scenario in SCENARIOS:
        res = np.zeros((len(HOLD_TEMPS), len(GRAINS)))
        resid = np.zeros((len(HOLD_TEMPS), len(GRAINS)))
        for j, (name, u_ppm, th_ppm, radius, obs_date, obs_1s) in enumerate(GRAINS):
            for i, t_hold in enumerate(HOLD_TEMPS):
                res[i, j] = model_date(u_ppm, th_ppm, radius, t_hold, scenario)
                no_inherit = model_date(u_ppm, th_ppm, radius, t_hold, scenario,
                                        inherited=False)
                resid[i, j] = res[i, j] - no_inherit
        results[scenario] = res
        residual[scenario] = resid

    # print results tables
    for scenario, res in results.items():
        print(f"\nModeled present-day ZHe date (Ma): {SCENARIOS[scenario]}")
        print("(residual: date difference vs. the same grain with no inherited He;")
        print(" a residual near zero marks complete loss of pre-depositional He)")
        print(f"{'T_peak (C)':>10}" + "".join(f"{name:>22}{'residual':>10}" for name, *_ in GRAINS))
        for i, t_hold in enumerate(HOLD_TEMPS):
            print(f"{t_hold:>10}" + "".join(
                f"{res[i, j]:>22.0f}{residual[scenario][i, j]:>10.0f}"
                for j in range(len(GRAINS))))
    print("\nObserved dates: " + ", ".join(
        f"{name}: {obs:.0f} +/- {err:.0f} Ma" for name, _, _, _, obs, err in GRAINS))

    # worst-case initial condition: He and damage retained since ca. 1830 Ma
    # (maximal initial reservoir; damage below the percolation threshold
    # increases retentivity for the moderate-eU grains)
    print("\nWorst-case end member: He and damage retained since ca. 1830 Ma "
          "(peak_900_850 scenario).\nA modeled date above the ca. 1052 Ma "
          "depositional age means a pre-depositional date could\nsurvive for "
          "some initial condition; complete loss for ANY initial condition is "
          "reached at\nthe temperature where all dates fall below it.")
    print(f"{'T_peak (C)':>10}" + "".join(f"{name:>24}" for name, *_ in GRAINS))
    for t_hold in [150, 160, 170, 175, 180]:
        row = f"{t_hold:>10}"
        for name, u_ppm, th_ppm, radius, obs_date, obs_1s in GRAINS:
            tT = build_tT_path(t_hold, "peak_900_850")
            tT = tT[tT[:, 0] <= FREDA_DEPOSITION_AGE]
            tT = np.vstack([tT, np.array([[1830, 20], [1850, 180]])])
            path = pyt.tT_path(tT)
            rho_r_array, relevant_tT = path.guenthner_anneal()
            grain = pyt.zircon(radius, LOG2_NODES, relevant_tT, rho_r_array,
                               u_ppm, th_ppm)
            row += f"{grain.zirc_date('guenthner', 'guenthner')[0]:>24.0f}"
        print(row)
    print("The measured dates older than the ca. 1052 Ma depositional age require "
          "that complete\nloss of pre-depositional He did not occur: peak "
          "Neoproterozoic burial temperatures\nremained below ca. 170 C.")

    # figure: modeled date (upper) and pre-depositional He residual (lower)
    fig, (ax, ax_r) = plt.subplots(2, 1, figsize=(7, 8.5), sharex=True,
                                   height_ratios=[1, 1])
    for j, (name, _, _, _, obs_date, obs_1s) in enumerate(GRAINS):
        line, = ax.plot(HOLD_TEMPS, results["peak_900_850"][:, j], marker="o", label=name)
        ax.plot(HOLD_TEMPS, results["protracted"][:, j], marker="o", markersize=3,
                linestyle="--", color=line.get_color())
        ax.axhline(obs_date, color=line.get_color(), linestyle=":", linewidth=0.8)
        ax_r.plot(HOLD_TEMPS, residual["peak_900_850"][:, j], marker="o",
                  color=line.get_color(), label=name)
        ax_r.plot(HOLD_TEMPS, residual["protracted"][:, j], marker="o", markersize=3,
                  linestyle="--", color=line.get_color())
    ax.axhspan(1050, 1080, color="#db806a", alpha=0.4)
    ax.text(38, 1065, "Freda deposition", fontsize=9, color="#b3563e", va="center")
    ax.set_ylabel("modeled present-day ZHe date (Ma)")
    ax.set_title("Survival of pre-depositional He in Freda detrital zircon\n"
                 "(solid: peak 900–850 Ma; dashed: protracted 950–650 Ma hold;\n"
                 "dotted lines = observed dates)",
                 fontsize=10)
    ax.legend(fontsize=9, loc="lower left")

    ax_r.axhline(0, color="black", linewidth=0.8)
    ax_r.axvline(170, color="gray", linestyle="--", linewidth=1)
    ax_r.text(171, ax_r.get_ylim()[1] * 0.6,
              "complete loss\nfor any initial He\nby ca. 170 °C",
              fontsize=9, color="gray", va="top")
    ax_r.set_xlabel("peak Neoproterozoic burial temperature (°C)")
    ax_r.set_ylabel("pre-depositional He signal (Myr)\n"
                    "(date difference vs. no inherited He)")
    fig.tight_layout()
    fig.savefig(FIGURES_DIR / "freda_zhe_burial_bound.pdf", bbox_inches="tight")
    print(f"\nFigure saved to {FIGURES_DIR / 'freda_zhe_burial_bound.pdf'}")

    # BRb syn-depositional burial sweep
    brb_grains = [g for g in GRAINS if g[0].startswith("BRb")]
    brb_results = {}
    brb_residual = {}
    for t_neo in NEO_LEVELS:
        res = np.zeros((len(SYN_FREDA_PEAK_TEMPS), len(brb_grains)))
        resid = np.zeros_like(res)
        for j, (name, u_ppm, th_ppm, radius, obs_date, obs_1s) in enumerate(brb_grains):
            for i, t_peak in enumerate(SYN_FREDA_PEAK_TEMPS):
                res[i, j] = model_brb_date(u_ppm, th_ppm, radius, t_peak, t_neo)
                no_inherit = model_brb_date(u_ppm, th_ppm, radius, t_peak, t_neo,
                                            inherited=False)
                resid[i, j] = res[i, j] - no_inherit
        brb_results[t_neo] = res
        brb_residual[t_neo] = resid

    print("\n\nBRb syn-depositional burial model: progressive burial beneath "
          f"{BRB_OVERBURDEN_KM:.0f} km of Freda strata\nfrom 1075 to 1050 Ma, "
          "rapid exhumation by 1030 Ma. Peak temperature at 1050 Ma implies\n"
          f"a geothermal gradient of (T - {SURFACE_TEMP})/{BRB_OVERBURDEN_KM:.0f} "
          "C/km for the known overburden.")
    for t_neo in NEO_LEVELS:
        print(f"\nNeoproterozoic burial fixed at {t_neo} C (900-850 Ma):")
        print(f"{'T_peak (C)':>10}{'gradient':>10}" + "".join(
            f"{name:>22}{'residual':>10}" for name, *_ in brb_grains))
        for i, t_peak in enumerate(SYN_FREDA_PEAK_TEMPS):
            gradient = (t_peak - SURFACE_TEMP) / BRB_OVERBURDEN_KM
            print(f"{t_peak:>10}{gradient:>10.0f}" + "".join(
                f"{brb_results[t_neo][i, j]:>22.0f}{brb_residual[t_neo][i, j]:>10.0f}"
                for j in range(len(brb_grains))))

    # figure: BRb syn-Freda burial sweep (date and residual panels)
    fig2, (bx, bx_r) = plt.subplots(2, 1, figsize=(7, 8.5), sharex=True)
    for j, (name, _, _, _, obs_date, obs_1s) in enumerate(brb_grains):
        line, = bx.plot(SYN_FREDA_PEAK_TEMPS, brb_results[NEO_LEVELS[0]][:, j],
                        marker="o", label=name)
        bx.plot(SYN_FREDA_PEAK_TEMPS, brb_results[NEO_LEVELS[1]][:, j], marker="o",
                markersize=3, linestyle="--", color=line.get_color())
        bx.axhline(obs_date, color=line.get_color(), linestyle=":", linewidth=0.8)
        bx_r.plot(SYN_FREDA_PEAK_TEMPS, brb_residual[NEO_LEVELS[0]][:, j],
                  marker="o", color=line.get_color(), label=name)
        bx_r.plot(SYN_FREDA_PEAK_TEMPS, brb_residual[NEO_LEVELS[1]][:, j],
                  marker="o", markersize=3, linestyle="--", color=line.get_color())
    bx.axhspan(1050, 1080, color="#db806a", alpha=0.4)
    bx.set_ylabel("modeled present-day ZHe date (Ma)")
    bx.set_title("BRb-1.2: syn-Freda burial (1075–1050 Ma) beneath 4 km of strata\n"
                 f"(solid: Neoproterozoic burial {NEO_LEVELS[0]} °C; "
                 f"dashed: {NEO_LEVELS[1]} °C;\ndotted lines = observed dates)",
                 fontsize=10)
    bx.legend(fontsize=9, loc="lower left")
    grad_axis = bx.secondary_xaxis(
        "top", functions=(lambda t: (t - SURFACE_TEMP) / BRB_OVERBURDEN_KM,
                          lambda g: SURFACE_TEMP + BRB_OVERBURDEN_KM * g))
    grad_axis.set_xlabel("implied geothermal gradient for 4 km of burial (°C/km)")
    # independent upper bound on the geotherm from basement Rb-Sr preservation,
    # marked on the shared gradient axis (T = surface + 4 km x gradient)
    rbsr_T = SURFACE_TEMP + BRB_OVERBURDEN_KM * RBSR_GEOTHERM_MAX
    for ax in (bx, bx_r):
        ax.axvspan(rbsr_T, SYN_FREDA_PEAK_TEMPS[-1] + 5, color="gray", alpha=0.12)
        ax.axvline(rbsr_T, color="#2e6f40", linestyle="--", linewidth=1.3)
    _ylo, _yhi = bx.get_ylim()
    bx.text(rbsr_T - 4, _ylo + 0.62 * (_yhi - _ylo),
            f"basement Rb-Sr:\ngeotherm ≲ {RBSR_GEOTHERM_MAX:.0f} °C/km",
            fontsize=8.5, color="#2e6f40", ha="right", va="center")
    bx_r.axhline(0, color="black", linewidth=0.8)
    bx_r.set_xlabel("peak syn-Freda burial temperature at ca. 1050 Ma (°C)")
    bx_r.set_ylabel("pre-depositional He signal (Myr)\n"
                    "(date difference vs. no inherited He)")
    fig2.tight_layout()
    fig2.savefig(FIGURES_DIR / "brb_syn_freda_burial.pdf", bbox_inches="tight")
    print(f"\nFigure saved to {FIGURES_DIR / 'brb_syn_freda_burial.pdf'}")


if __name__ == "__main__":
    main()
