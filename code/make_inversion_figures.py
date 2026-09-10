"""Assemble the composite time-temperature inversion figures for the manuscript.

The Thermochron.jl inversions in ``inversion_files/`` write one posterior
path-density map per chain (``<prefix>_tT.pdf``). Those single-chain plots
are archived with each run, but the manuscript figures place several chains
side by side on common axes with a shared colour scale. This script builds
those composites from the density arrays exported by
``inversion_files/export_tT_density.jl``:

* ``<prefix>_tT_density.npy`` -- 8-bit density image (1200 x 1800), rows from
  0 to 400 °C (top to bottom), columns from 0 to 2704 Ma (left to right),
  normalised to the 98.5th percentile of the image exactly as in
  ``tTinversion_mrm.jl``. These are intermediate files (git-ignored).
* ``<prefix>_tT_density.json`` -- time and temperature ranges and the
  post-burn-in chain statistics.

On its first run the script converts each ``.npy`` to a lossless 8-bit PNG
(``<prefix>_tT_density.png``) next to it; the PNGs are what is archived in
the repository, and later runs read them directly so that the figures can be
regenerated without the (unarchived) ``.jld2`` posteriors.

Figures written to ``figures/``:

* ``inversion_tT_1500Ma.png`` -- main-text figure: the three constraint
  configurations at the reference seed (20260710), one row, shown from
  1500 Ma to the present where the data have resolving power.
* ``inversion_tT_1500Ma_vertical.png`` -- the same three panels stacked,
  for a single-column layout.
* ``inversion_tT.png`` and ``inversion_tT_vertical.png`` -- the same panels
  over the full model span (2704 Ma to the present), for the supporting
  information.
* ``inversion_seed_comparison.png`` -- supporting information: the three
  configurations (rows) at the three seeds (columns).
* ``inversion_grain_sensitivity.png`` -- supporting information: the rift
  burial + Cambrian configuration at seed 20260710 with the full dataset and
  with zircon analyses TAI_MRM4_z3 and TAI_MRM4_z5 removed
  (``inversion_files/swansonhysell_runs_2026-08b/``).

The colour scale is the ``ylcn`` colormap of StatGeochemBase.jl used by
Thermochron.jl, transcribed here so the composites match the single-chain
plots. Constraint boxes are drawn from the same uniform priors defined in
``tTinversion_mrm.jl``.

Run from anywhere::

    python code/make_inversion_figures.py
"""

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Rectangle
import numpy as np
from PIL import Image

REPO = Path(__file__).resolve().parent.parent
INVERSION = REPO / "inversion_files"
FIGURES = REPO / "figures"

BATCH_A = INVERSION / "swansonhysell_runs_2026-08a"
BATCH_B = INVERSION / "swansonhysell_runs_2026-08b"

SEEDS = [20260710, 12345, 777]
CONFIGS = ["unconstrained", "no_jacobsville", "all_constraints"]
CONFIG_TITLES = {
    "unconstrained": "Unconstrained",
    "no_jacobsville": "Rift burial + Cambrian constraints",
    "all_constraints": "Rift burial + Jacobsville + Cambrian constraints",
}

# Uniform (box) priors from tTinversion_mrm.jl: (tmin, tmax, Tmin, Tmax)
BOXES = {
    "mcr_burial": (1050, 1100, 250, 400),
    "jacobsville": (980, 1000, 20, 50),
    "cambrian": (485, 505, 20, 50),
}
CONFIG_BOXES = {
    "unconstrained": [],
    "no_jacobsville": ["mcr_burial", "cambrian"],
    "all_constraints": ["mcr_burial", "jacobsville", "cambrian"],
}

# StatGeochemBase.jl `ylcn` colormap (256 colours), as used by Thermochron.jl
YLCN_HEX = (
    "7AFEFF 78FCFF 76FAFF 74F8FF 72F6FF 71F4FF 6FF2FF 6DF0FF 6BEFFF 69EDFF 67EBFF "
    "65E9FF 63E7FF 61E5FF 5FE3FF 5DE1FF 5BDFFF 5ADDFF 58DBFF 56D9FF 54D7FF 52D5FF "
    "50D3FF 4ED1FF 4CCFFF 4ACDFF 48CBFF 46C9FF 44C7FF 43C5FF 41C3FF 3FC1FF 3DBFFF "
    "3BBDFF 39BBFF 37BAFF 35B8FF 33B6FF 31B4FF 2FB2FF 2DB0FF 2CAEFF 2AACFF 28AAFF "
    "26A8FF 24A6FF 22A4FF 20A2FF 1EA0FF 1C9EFF 1A9CFF 189AFF 1798FF 1596FF 1394FF "
    "1192FF 0F90FF 0D8EFF 0B8CFF 098AFF 0788FF 0586FF 0384FF 0183FF 0181FF 037FFD "
    "057CFB 077AF9 0978F8 0B76F6 0D74F4 0F72F2 1170F0 136EEE 156CEC 176AEA 1968E8 "
    "1B66E6 1D64E4 1F62E2 2160E0 235EDE 255CDC 275ADA 2958D9 2C56D7 2E54D5 3052D3 "
    "3250D1 344ECF 364CCD 384ACB 3A48C9 3C46C7 3E44C5 4042C3 4240C1 443EBF 463CBD "
    "483ABB 4A38B9 4C36B8 4E34B6 5032B4 5230B2 552EB0 572CAE 5929AC 5B27AA 5D25A8 "
    "5F23A6 6121A4 631FA2 651DA0 671B9E 69199C 6B179A 6D1598 6F1397 711195 730F93 "
    "750D91 770B8F 79098D 7B078B 7D0589 800387 820185 840183 860381 88057F 89077D "
    "8B097B 8D0C79 8F0E77 911075 931273 951471 97166F 99186C 9B1A6A 9D1C68 9F1E66 "
    "A12164 A32362 A52560 A7275E A9295C AB2B5A AD2D58 AF2F56 B13154 B33352 B53550 "
    "B7384E B93A4B BB3C49 BD3E47 BE4045 C04243 C24441 C4463F C6483D C84A3B CA4D39 "
    "CC4F37 CE5135 D05333 D25531 D4572F D6592D D85B2A DA5D28 DC5F26 DE6224 E06422 "
    "E26620 E4681E E66A1C E86C1A EA6E18 EC7016 EE7214 F07412 F27610 F3790E F57B0C "
    "F77D0A F97F07 FB8105 FD8303 FF8501 FF8702 FF8904 FF8B06 FF8D08 FF8F0A FF910C "
    "FF930E FF9510 FF9612 FF9814 FF9A16 FF9C18 FF9E1A FFA01C FFA21E FFA420 FFA622 "
    "FFA824 FFAA26 FFAC28 FFAD2A FFAF2C FFB12E FFB330 FFB532 FFB734 FFB936 FFBB38 "
    "FFBD3A FFBF3C FFC13E FFC240 FFC442 FFC644 FFC846 FFCA48 FFCC4A FFCE4C FFD04E "
    "FFD250 FFD452 FFD654 FFD856 FFD958 FFDB5A FFDD5C FFDF5E FFE160 FFE362 FFE564 "
    "FFE766 FFE968 FFEB6A FFED6C FFEF6E FFF070 FFF272 FFF474 FFF676 FFF878 FFFA7A "
    "FFFC7C FFFE7E FFFF80"
).split()
YLCN = LinearSegmentedColormap.from_list("ylcn", ["#" + h for h in YLCN_HEX], N=256)


def prefix(config, seed, batch=BATCH_A, suffix=""):
    """Output-file prefix of a chain, e.g. .../MRM_all_constraints_seed777."""
    return batch / f"MRM_{config}_seed{seed}{suffix}"


def load_density(prefix_path):
    """Load a chain's density map and metadata, archiving the PNG if needed.

    Args:
        prefix_path: Path prefix of the chain's output files.

    Returns:
        Tuple of (density array scaled 0-1, metadata dict).
    """
    png = Path(str(prefix_path) + "_tT_density.png")
    npy = Path(str(prefix_path) + "_tT_density.npy")
    meta = json.loads(Path(str(prefix_path) + "_tT_density.json").read_text())
    if npy.exists():
        A = np.load(npy)
        assert A.dtype == np.uint8 and list(A.shape) == meta["shape"]
        if not png.exists():
            Image.fromarray(A, mode="L").save(png, optimize=True)
            print(f"archived {png.relative_to(REPO)}")
    elif png.exists():
        A = np.asarray(Image.open(png))
    else:
        raise FileNotFoundError(f"no density map for {prefix_path}")
    return A.astype(float) / 255.0, meta


def draw_panel(ax, density, meta, boxes, title=None, label=None,
               xlabel=True, ylabel=True, title_pad=None, tmax=None):
    """Draw one posterior path-density map on an axis.

    Args:
        ax: Matplotlib axis.
        density: 2-D array scaled 0-1, rows 0 to Tinit, columns 0 to tinit.
        meta: Metadata dict with ``time_Ma`` and ``temperature_C`` ranges.
        boxes: Names of the constraint boxes to draw (keys of ``BOXES``).
        title: Optional panel title.
        label: Optional bold panel letter placed above the axis.
        xlabel, ylabel: Whether to label the axes.
        title_pad: Title offset in points (default clears the panel label).
        tmax: Oldest time (Ma) to show; default is the full model span.
    """
    t0, t1 = meta["time_Ma"]
    T0, T1 = meta["temperature_C"]
    ax.imshow(density, cmap=YLCN, vmin=0, vmax=1, aspect="auto",
              origin="upper", extent=(t0, t1, T1, T0),
              interpolation="nearest")
    for name in boxes:
        box_t0, box_t1, box_T0, box_T1 = BOXES[name]
        ax.add_patch(Rectangle((box_t0, box_T0), box_t1 - box_t0, box_T1 - box_T0,
                               fill=False, edgecolor="black", linewidth=1.2))
    if tmax is None:
        tmax = t1
    ax.set_xlim(tmax, t0)
    ax.set_ylim(T1, T0)
    ax.set_xticks(np.arange(0, tmax + 1, 500 if tmax > 1500 else 200))
    ax.set_yticks(np.arange(T0, T1 + 1, 50))
    ax.minorticks_on()
    ax.tick_params(direction="out", which="both", top=True, right=True,
                   labelsize=8)
    if xlabel:
        ax.set_xlabel("Time (Ma)", fontsize=9)
    if ylabel:
        ax.set_ylabel("Temperature (°C)", fontsize=9)
    if title:
        if title_pad is None:
            title_pad = 18 if label else 6
        ax.set_title(title, fontsize=10, pad=title_pad)
    if label:
        ax.text(0, 1.03, label, transform=ax.transAxes, fontsize=10,
                fontweight="bold", ha="left", va="bottom")


def add_colorbar(fig, axes, label="Normalized path density"):
    """Add a shared colour bar to the right of a set of axes."""
    cbar = fig.colorbar(plt.cm.ScalarMappable(cmap=YLCN),
                        ax=axes, fraction=0.03, pad=0.02, shrink=0.9)
    cbar.set_label(label, fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    return cbar


def stats_line(meta):
    """One-line summary of a chain's post-burn-in statistics."""
    return (f"mean log-likelihood {meta['mean_loglikelihood']:.1f}, "
            f"acceptance {meta['mean_acceptance']:.2f}, "
            f"nodes {meta['mean_npoints']:.1f}")


def main_text_figure(seed=SEEDS[0], vertical=False, tmax=None):
    """Three constraint configurations at one seed (Figure: inversion_tT).

    Args:
        seed: Random seed of the chains to show.
        vertical: Stack the panels in one column instead of one row.
        tmax: Oldest time (Ma) to show; None for the full model span.
    """
    if vertical:
        fig, axes = plt.subplots(3, 1, figsize=(5.4, 10.8))
        fig.subplots_adjust(hspace=0.22)
    else:
        fig, axes = plt.subplots(1, 3, figsize=(13.5, 3.6))
    for ax, config, letter in zip(axes, CONFIGS, "abc"):
        density, meta = load_density(prefix(config, seed))
        draw_panel(ax, density, meta, CONFIG_BOXES[config],
                   label=f"({letter}) {CONFIG_TITLES[config]}",
                   xlabel=(not vertical or ax is axes[-1]),
                   ylabel=(vertical or ax is axes[0]), tmax=tmax)
        print(f"  {config} seed {seed}: {stats_line(meta)}")
    add_colorbar(fig, list(axes))
    name = "inversion_tT" + (f"_{tmax:.0f}Ma" if tmax else "") + ("_vertical" if vertical else "")
    out = FIGURES / f"{name}.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.relative_to(REPO)}")


def seed_comparison_figure():
    """Three configurations (rows) at three seeds (columns) (Figure S5)."""
    fig, axes = plt.subplots(3, 3, figsize=(13.5, 10.2))
    for i, (config, letter) in enumerate(zip(CONFIGS, "abc")):
        for j, seed in enumerate(SEEDS):
            ax = axes[i, j]
            density, meta = load_density(prefix(config, seed))
            draw_panel(ax, density, meta, CONFIG_BOXES[config],
                       title=(f"seed {seed}" if i == 0 else None),
                       title_pad=18,
                       label=(f"({letter}) {CONFIG_TITLES[config]}"
                              if j == 0 else None),
                       xlabel=(i == 2), ylabel=(j == 0))
            print(f"  {config} seed {seed}: {stats_line(meta)}")
    add_colorbar(fig, axes.ravel().tolist())
    out = FIGURES / "inversion_seed_comparison.png"
    fig.savefig(out, dpi=250, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.relative_to(REPO)}")


def grain_sensitivity_figure(seed=SEEDS[0]):
    """Rift burial + Cambrian chain with and without TAI_MRM4_z3 and _z5."""
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    panels = [
        (prefix("no_jacobsville", seed), "All 55 zircon and 13 apatite analyses"),
        (prefix("no_jacobsville", seed, BATCH_B, "_no_MRM4_z3_z5"),
         "TAI_MRM4_z3 and TAI_MRM4_z5 removed"),
    ]
    for ax, (pfx, title), letter in zip(axes, panels, "ab"):
        density, meta = load_density(pfx)
        draw_panel(ax, density, meta, CONFIG_BOXES["no_jacobsville"],
                   label=f"({letter}) {title}", ylabel=(ax is axes[0]))
        print(f"  {pfx.name}: {stats_line(meta)}")
    add_colorbar(fig, list(axes))
    out = FIGURES / "inversion_grain_sensitivity.png"
    fig.savefig(out, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {out.relative_to(REPO)}")


if __name__ == "__main__":
    FIGURES.mkdir(exist_ok=True)
    print("main-text figure (1500 Ma to present)")
    main_text_figure(tmax=1500)
    main_text_figure(vertical=True, tmax=1500)
    print("full-span three-panel figure")
    main_text_figure()
    main_text_figure(vertical=True)
    print("seed comparison figure")
    seed_comparison_figure()
    print("grain sensitivity figure")
    grain_sensitivity_figure()
