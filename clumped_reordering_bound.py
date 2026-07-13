"""Clumped-isotope reordering bound on Neoproterozoic burial temperature.

Vein calcite sample AMF6 from the Montreal River region (U-Pb age 1049 +/- 13
Ma) yields a Delta_47-based (apparent) temperature of 54.4 +/- 1.7 degrees C.
\\citeA{Hodgin2024a} interpret this as a maximum temperature of formation at
1052-1036 Ma, since low-temperature solid-state reordering raises the apparent
temperature without resetting the U-Pb age (Passey and Henkes, 2012; Stolper
and Eiler, 2015). The apparent temperature therefore caps peak Neoproterozoic
burial: burial hot enough to reorder the calcite would have driven its apparent
temperature above the observed value.

The bound is evaluated conservatively along three axes:

  1. Formation temperature. The calcite may have precipitated near 54.4 degrees
     C with negligible reordering, or precipitated cold (near-surface, following
     structural uplift) and been reordered upward during burial. A colder
     formation retains its ordered signal to a higher burial temperature, so the
     coldest plausible formation sets the loosest bound. A family of formation
     temperatures (20-54.4 degrees C) is swept.

  2. Reordering kinetics. Only the Passey & Henkes (2012) calcite model is used.
     The Stolper & Eiler (2015) routine in the ClumpyCool distribution fails an
     isothermal-equilibration sanity check in this implementation (an isothermal
     hold at 160 degrees C for 1 Gyr yields only ~74 degrees C apparent, whereas
     the signal must equilibrate to the hold temperature), so it is excluded.
     Passey & Henkes also reorders less readily in the relevant range and thus
     gives the more conservative (higher) bound.

  3. Thermal-history scenario. The burial/exhumation history is unknown, so two
     scenarios bracketing slow and rapid exhumation are modeled (see SCENARIOS).
     A slowly cooled sample partially re-equilibrates its clumped signal downward
     during cooling and so tolerates a higher peak temperature before its
     apparent temperature exceeds the measured value; the slow-exhumation
     scenario therefore sets the conservative (highest) bound.

For each scenario the calcite starts at its formation Delta_47 and reordering is
forward modeled along the burial path with a swept peak temperature. The peak
temperature at which the modeled present-day apparent temperature rises above
the measured 54.4 +/- 1.7 degrees C is the upper bound on Neoproterozoic burial.

Dependency: ClumpyCool (Lloyd, 2020; doi:10.17605/OSF.IO/JYHSW). A clone of
https://github.com/maxmansaxman/Clumpy_reordering_model_distribution in this
directory (as Clumpy_reordering_model_distribution/) is found automatically;
alternatively place ClumpyCool_func.py anywhere on the Python path.

Usage:
    python clumped_reordering_bound.py

Output:
    - printed tables of modeled apparent temperature and crossing temperatures
    - figures/clumped_reordering_bound.pdf (4 thermal-history scenarios)
"""

import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import PchipInterpolator

try:
    import ClumpyCool_func as cc
except ImportError:  # fall back to a local clone of the ClumpyCool distribution
    sys.path.append(str(Path(__file__).parent / "Clumpy_reordering_model_distribution"))
    import ClumpyCool_func as cc

MYR = 1e6

# AMF6 clumped isotope constraint (Hodgin et al., 2024)
FORMATION_AGE = 1049      # Ma, U-Pb age
MEASURED_TEMP = 54.4      # degrees C, apparent Delta_47 temperature
MEASURED_1S = 1.7         # degrees C, 1 sigma

KINETIC_MODEL = "p"       # Passey & Henkes (2012); see module docstring

FORMATION_TEMPS = [20, 30, 40, 54.4]   # degrees C
PEAK_TEMPS = [35, 60, 80, 90, 100, 105, 110, 115, 120, 125, 130, 135, 140, 150, 160]

# Four Neoproterozoic burial/exhumation scenarios. Each: ordered (age_Ma,
# temperature_C) nodes with None at the swept peak node(s). Calcite forms at
# 1049 Ma; ages are converted to elapsed time (years since formation).
SCENARIOS = {
    "slow exhumation\n(peak 980 Ma, cool to 505 Ma)": {
        "ages": [1049, 1010, 980, 505, 250, 0],
        "temps": [35, 35, None, 35, 50, 10],
    },
    "rapid exhumation\n(hold 950-750 Ma, exhumed by 700 Ma)": {
        "ages": [1049, 1010, 950, 750, 700, 505, 250, 0],
        "temps": [35, 35, None, None, 20, 35, 50, 10],
    },
}


def apparent_temperature(scenario, t_form, t_peak):
    """Model the present-day apparent Delta_47 temperature for a burial path.

    Args:
        scenario: Key into SCENARIOS defining the burial-path node structure.
        t_form: Calcite formation temperature (degrees C), setting initial Delta_47.
        t_peak: Peak Neoproterozoic burial temperature (degrees C).

    Returns:
        Modeled apparent (Delta_47-based) present-day temperature (degrees C).
    """
    spec = SCENARIOS[scenario]
    temps = [t_peak if t is None else t for t in spec["temps"]]
    elapsed_yrs = (FORMATION_AGE - np.array(spec["ages"], dtype=float)) * MYR
    tT_fn = PchipInterpolator(elapsed_yrs, np.array(temps, dtype=float))
    time_array = np.linspace(0, elapsed_yrs.max(), 5000)
    d47_init = cc.temp_to_D47_ARF(t_form)
    pred, _ = cc.predict_D47(time_array, tT_fn, d47_init,
                             model=KINETIC_MODEL, mineral="calcite_arf", sigma_choice=1)
    return pred["Temp_pred"].iloc[-1]


def crossing(apparent, threshold):
    """Linearly interpolate the peak temperature at which apparent = threshold."""
    a = np.array(apparent)
    for i in range(1, len(a)):
        if a[i - 1] < threshold <= a[i]:
            f = (threshold - a[i - 1]) / (a[i] - a[i - 1])
            return PEAK_TEMPS[i - 1] + f * (PEAK_TEMPS[i] - PEAK_TEMPS[i - 1])
    return None


def main():
    results = {s: {t: [apparent_temperature(s, t, tp) for tp in PEAK_TEMPS]
                   for t in FORMATION_TEMPS}
               for s in SCENARIOS}

    print(f"\nAMF6 measured Delta_47 temperature: {MEASURED_TEMP} +/- {MEASURED_1S} "
          f"degrees C (1 sigma); U-Pb age {FORMATION_AGE} Ma (Passey & Henkes kinetics)")
    conservative = 0
    for s in SCENARIOS:
        label = s.replace("\n", " ")
        print(f"\n  scenario: {label}")
        print("    peak burial T where apparent exceeds 54.4 + 1.7 C (1 sigma):")
        for t in FORMATION_TEMPS:
            c = crossing(results[s][t], MEASURED_TEMP + MEASURED_1S)
            print(f"      formation {t:>6} C:  ~{c:.0f} C" if c
                  else f"      formation {t:>6} C:  not reached")
            if c:
                conservative = max(conservative, c)
    print(f"\nConservative upper bound (highest crossing across scenarios and "
          f"formation temperatures): ~{conservative:.0f} degrees C")

    # figure: one panel per thermal-history scenario
    n = len(SCENARIOS)
    fig, axes = plt.subplots(1, n, figsize=(6 * n, 5.5), sharex=True, sharey=True)
    axes = np.atleast_1d(axes)
    cmap = plt.get_cmap("viridis")
    for ax, s in zip(axes, SCENARIOS):
        for k, t in enumerate(FORMATION_TEMPS):
            ax.plot(PEAK_TEMPS, results[s][t], marker="o", markersize=4,
                    color=cmap(k / max(1, len(FORMATION_TEMPS) - 1)),
                    label=f"{t} °C")
        ax.axhspan(MEASURED_TEMP - MEASURED_1S, MEASURED_TEMP + MEASURED_1S,
                   color="gray", alpha=0.3)
        ax.axhline(MEASURED_TEMP, color="gray", linewidth=1)
        ax.annotate(f"measured Δ$_{{47}}$ T = {MEASURED_TEMP} ± {MEASURED_1S} °C\n"
                    "(Hodgin et al., 2024)",
                    xy=(PEAK_TEMPS[0] + 3, MEASURED_TEMP + MEASURED_1S),
                    xytext=(PEAK_TEMPS[0] + 3, 95), fontsize=8.5, color="dimgray",
                    arrowprops=dict(arrowstyle="->", color="gray", lw=0.8))
        ax.set_title(s, fontsize=10)
        ax.set_xlabel("peak Neoproterozoic burial temperature (°C)")
        ax.legend(fontsize=8, loc="lower right", title="calcite\nformation T",
                  title_fontsize=8)
    axes[0].set_ylabel("modeled apparent Δ$_{47}$ temperature, present day (°C)")
    fig.suptitle("Clumped-isotope reordering limit on Neoproterozoic burial "
                 "(AMF6 vein calcite; Passey & Henkes 2012 kinetics)", fontsize=12)
    fig.tight_layout()
    fig.savefig("figures/clumped_reordering_bound.pdf", bbox_inches="tight")
    print("\nFigure saved to figures/clumped_reordering_bound.pdf")


if __name__ == "__main__":
    main()
