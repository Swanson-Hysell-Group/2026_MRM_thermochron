"""Build LaTeX supplementary data tables from the (U-Th)/He data files.

This script converts the reduced single-grain data in

    data/UTh_He_apatite_data.csv
    data/UTh_He_zircon_data.csv

into two ``longtable`` LaTeX files that are input by the supporting information
document si_materials/si_materials.tex:

    si_materials/si_UTh_He_apatite_table.tex
    si_materials/si_UTh_He_zircon_table.tex

The data treatment applied here mirrors that of code/MRM_UTh_He_data.ipynb and
code/MRM_forward_models.ipynb:

* Analysis names are split into a sample name and a grain name. Sample ``397``
  is the Mellen Complex sample MC4, sample ``BRB`` is BRb-1.2 and sample ``LRA``
  (analysis names ``TAI_LRA_67_z*``) is LRa-67.0.
* Samples are assigned to the groups defined in the notebooks: Archean basement,
  Midcontinent Rift, and Freda Formation detrital zircon.
* Grains that the notebooks filter out as spurious (TA_MRM20_z3 and
  TAI_MRM34_z2) are retained in the table but flagged as excluded, with the
  reason given, so that the reported dataset is complete.
* Analyses of samples belonging to other projects that share the same
  laboratory report files (APR, DKR1, DTR, MG, MG18, MA5 zircon; MG_G1 apatite)
  are dropped. Note that TAI_MA5_z1 and TAI_MA5_z2 are named in the notebook
  filtering step but MA5 is not one of the samples of this study, so those
  grains do not appear here either way. Two further apatite analyses
  (TAI_MG_G1_a2 and TAI_BRB_a1) were removed at the export step in
  MRM_UTh_He_data.ipynb and are therefore not present in the csv file.
* Fish Canyon Tuff zircon (FCT) and Durango apatite (DUR) aliquots run as
  secondary standards in the same analytical sessions are reported in a final
  block of each table.
* Grains flagged as ``not included in the inversion dataset`` are determined by
  comparing against inversion_files/mrm.csv, the input file of the
  Thermochron.jl inversions, rather than being hard coded.

Notes on individual analyses are carried as superscript footnote symbols in the
last column, with the symbols defined in a legend beneath each table, so that no
row has to wrap onto a second line.

Several zircon grains have two rows in the csv file: the He measurement, mass,
radius and F_T are identical but the U-Th determination differs between
laboratory report files, so the two rows give different dates. Both rows are
carried through the notebooks (both are in the inversion input file), so both
are reported here, distinguished by the report file they come from and flagged
with a dagger.

Run from anywhere:

    python code/make_UTh_He_tables.py
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import string

REPO = Path(__file__).resolve().parent.parent
DATA = REPO / "data"
APATITE_CSV = DATA / "UTh_He_apatite_data.csv"
ZIRCON_CSV = DATA / "UTh_He_zircon_data.csv"
INVERSION_CSV = REPO / "inversion_files" / "mrm.csv"
SI_DIR = REPO / "si_materials"
APATITE_TEX = SI_DIR / "si_UTh_He_apatite_table.tex"
ZIRCON_TEX = SI_DIR / "si_UTh_He_zircon_table.tex"

# analysis-name sample codes mapped to the field sample names used in the paper
SAMPLE_DISPLAY = {
    "397": "MC4",
    "BRB": "BRb-1.2",
    "LRA_67": "LRa-67.0",
    "MG_G1": "MG-G1",
}

# sample groups following MRM_UTh_He_data.ipynb and MRM_forward_models.ipynb
ZIRCON_GROUPS = [
    ("Archean basement",
     ["MRM1", "MRM2", "MRM4", "MRM5", "MRM6", "MRM7", "MRM8", "MRM16",
      "MRM17", "MRM20", "MRM21", "MRM24", "MRM35", "MRM36"]),
    ("Midcontinent Rift", ["MRM9", "MRM34", "397"]),
    ("Freda Formation detrital zircon", ["LRA_67", "BRB"]),
    ("Fish Canyon Tuff zircon standard", ["FCT"]),
]

APATITE_GROUPS = [
    ("Archean basement",
     ["MRM6", "MRM8", "MRM16", "MRM20", "MRM21", "MRM24", "MRM35", "MRM36"]),
    ("Midcontinent Rift (Mellen Complex)", ["397"]),
    ("Durango apatite standard", ["DUR"]),
]

# grains filtered out in the notebooks, with the stated reason
EXCLUDED = {
    "TA_MRM20_z3": "excluded: anomalous U--Th determination",
    #"TAI_MRM34_z2": "excluded: date exceeds crystallization age",
}

# groups whose grains are expected in the Thermochron.jl inversion input, so
# that a grain that is absent can be flagged
INVERSION_GROUPS = {"Archean basement", "Midcontinent Rift (Mellen Complex)"}

# morphology shorthand in the morph_notes column that describes eU rather than
# grain morphology and so is not carried into the table
EU_SHORTHAND = re.compile(r"^(very |low|mod|high|to |-|,|\s)*(eu)?[\s,-]*$", re.I)

# near-duplicate wordings of the same laboratory morphology note, mapped onto a
# single wording so that they share one table footnote
NOTE_NORMALIZATION = {
    "one tip rounded": "one tip is rounded",
}

# notes assigned by this script; laboratory notes read from the data files are
# escaped for LaTeX where they are read
REPEAT_NOTE = (r"two U--Th determinations reported for the same He aliquot, "
               r"listed here as (a) and (b)")
NOT_MODELED_NOTE = r"not included in the Thermochron.jl inversion dataset"

# footnote symbols, assigned in the order in which the notes first appear
SYMBOLS = ["*", r"\dag", r"\ddag", r"\S", r"\P", r"\#",
           "**", r"\dag\dag", r"\ddag\ddag", r"\S\S", r"\P\P", r"\#\#"]

EMDASH = "---"


def parse_analysis_name(name):
    """Split a laboratory analysis name into sample and grain names.

    Analysis names take the form ``TAI_MRM35_a4`` or ``TAI_LRA_67_z1``, where
    the leading token is the analytical session prefix, the trailing token is
    the grain, and everything in between is the sample.

    Args:
        name (str): Analysis name from the ``Unnamed: 0`` column.

    Returns:
        tuple: (sample code, grain name) with the grain name lower cased.
    """
    tokens = name.split("_")
    if tokens[0] in ("TAI", "TA"):
        tokens = tokens[1:]
    grain = tokens[-1].lower()
    sample = "_".join(tokens[:-1])
    return sample, grain

def undupe_grain_names(df, sample_col="sample", grain_col="grain"):
    """Give grains unique names when the same sample/grain combo shows up
    more than once (e.g. from TAI_ and TA_ prefixed duplicates).

    Args:
        df (pandas.DataFrame): DataFrame containing the sample and grain
            columns produced by "parse_analysis_name".
        sample_col (str): Name of the column holding sample codes.
            Defaults to "sample".
        grain_col (str): Name of the column holding grain names. Defaults
            to "grain".

    Returns:
        pandas.DataFrame: The same DataFrame with "grain_col" updated so
        that repeated sample/grain combinations get a lower-case letter
        suffix (first occurrence unchanged, then 'a', 'b', 'c', ...).
    """
    combo = df[sample_col] + "_" + df[grain_col]

    # 0 = first occurrence, 1 = second, 2 = third, ...
    occurrence = combo.groupby(combo).cumcount()

    # 0 -> '', 1 -> 'a', 2 -> 'b', ...
    suffix = occurrence.map(lambda n: "" if n == 0 else string.ascii_lowercase[n - 1])

    df[grain_col] = df[grain_col] + suffix
    return df

def report_number(url):
    """Return the HAL report file number that an analysis was taken from."""
    if not isinstance(url, str):
        return None
    match = re.search(r"MCR_(\d+)_report", url)
    return int(match.group(1)) if match else None


def latex_escape(text):
    """Escape LaTeX special characters in a note string."""
    if not isinstance(text, str):
        return ""
    replacements = {
        "\\": r"\textbackslash{}", "&": r"\&", "%": r"\%", "$": r"\$",
        "#": r"\#", "_": r"\_", "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for old, new in replacements.items():
        text = text.replace(old, new)
    # laboratory notes write temperatures as e.g. 750C
    text = re.sub(r"(\d+)\s?C\b", r"\1\\,\\textdegree C", text)
    return text


def sig(value, digits=4):
    """Format a number to a fixed number of significant figures for LaTeX.

    Values that are very large or very small are written in scientific
    notation so that the anomalous U-Th determinations remain legible.

    Args:
        value (float): Value to format.
        digits (int): Number of significant figures.

    Returns:
        str: LaTeX-ready string, or an em dash if the value is missing.
    """
    if value is None or (isinstance(value, float) and not np.isfinite(value)):
        return EMDASH
    if value == 0:
        return "0"
    if abs(value) >= 1e5 or abs(value) < 1e-3:
        mantissa, exponent = f"{value:.{digits - 1}e}".split("e")
        return f"${mantissa}\\times10^{{{int(exponent)}}}$"
    return f"{value:.{digits}g}"


def fmt_percent(percent):
    """Format a relative uncertainty in percent with a sensible precision."""
    if percent >= 10:
        return f"{percent:.0f}"
    if percent >= 1:
        return f"{percent:.1f}"
    return f"{percent:.2f}"


def fmt_amount(value, sigma, digits=4):
    """Format a measured amount with its 1 sigma uncertainty as a percentage.

    Args:
        value (float): Measured amount (pmol or ng).
        sigma (float): Absolute 1 sigma uncertainty on the amount.
        digits (int): Significant figures on the amount.

    Returns:
        str: e.g. ``0.2045 (0.52)`` where the value in parentheses is 1 sigma
        as a percentage of the amount.
    """
    text = sig(value, digits)
    if text == EMDASH or value is None or value == 0:
        return text
    if sigma is None or not np.isfinite(sigma):
        return text
    return f"{text} ({fmt_percent(sigma / abs(value) * 100)})"


def fmt_date(value, sigma):
    """Format a date and its 1 sigma uncertainty in Ma to one decimal place."""
    if value is None or not np.isfinite(value):
        return EMDASH
    return f"{value:.1f} $\\pm$ {sigma:.1f}"


def fmt_fixed(value, decimals=1):
    """Format a number to a fixed number of decimal places."""
    if value is None or not np.isfinite(value):
        return EMDASH
    return f"{value:.{decimals}f}"


def clean_morph_note(note):
    """Return morphology notes, dropping the laboratory eU shorthand."""
    if not isinstance(note, str) or not note.strip():
        return ""
    if EU_SHORTHAND.match(note.strip()):
        return ""
    # strip any leading eU shorthand from a compound note
    parts = [p.strip() for p in note.split(",")]
    kept = [NOTE_NORMALIZATION.get(p, p)
            for p in parts if p and not EU_SHORTHAND.match(p)]
    return "; ".join(kept)


def prepare(csv_path, groups, inversion_names):
    """Read a data csv file and return the analyses of this study, annotated.

    Args:
        csv_path (Path): Path to the reduced data csv file.
        groups (list): List of (group name, sample code list) tuples defining
            which samples belong to this study and how they are grouped.
        inversion_names (set): Analysis names present in the Thermochron.jl
            inversion input file.

    Returns:
        pandas.DataFrame: The analyses of this study with sample, grain, group,
        report number and note columns added, sorted for presentation.
    """
    df = pd.read_csv(csv_path)
    df["csv_order"] = np.arange(len(df))
    parsed = df["Unnamed: 0"].apply(parse_analysis_name)
    df["sample_code"] = [p[0] for p in parsed]
    df["grain"] = [p[1] for p in parsed]
    df = undupe_grain_names(df, sample_col="sample_code", grain_col="grain")
    df["sample"] = df["sample_code"].map(lambda s: SAMPLE_DISPLAY.get(s, s))
    df["report"] = df["source_url"].apply(report_number) if "source_url" in df else None

    group_of = {}
    order_of = {}
    for group_index, (group_name, samples) in enumerate(groups):
        for sample_index, sample in enumerate(samples):
            group_of[sample] = group_name
            order_of[sample] = (group_index, sample_index)
    df = df[df["sample_code"].isin(group_of)].copy()
    df["group"] = df["sample_code"].map(group_of)
    df["group_order"] = df["sample_code"].map(lambda s: order_of[s][0])
    df["sample_order"] = df["sample_code"].map(lambda s: order_of[s][1])
    df["grain_order"] = df["grain"].str.extract(r"(\d+)").astype(float)

    # grains with more than one U-Th determination on the same He aliquot
    counts = df["Unnamed: 0"].value_counts()
    df["repeat_UTh"] = df["Unnamed: 0"].map(counts) > 1

    notes = []
    for _, row in df.iterrows():
        parts = []
        if row["Unnamed: 0"] in EXCLUDED:
            parts.append(EXCLUDED[row["Unnamed: 0"]])
        elif (row["group"] in INVERSION_GROUPS
              and row["Unnamed: 0"] not in inversion_names):
            parts.append(NOT_MODELED_NOTE)
        if row["repeat_UTh"]:
            parts.append(REPEAT_NOTE)
        morph = clean_morph_note(row.get("morph_notes"))
        if morph:
            parts.append(latex_escape(morph))
        analysis_note = row.get("analysis_notes")
        if isinstance(analysis_note, str) and analysis_note.strip():
            # quoted verbatim so that the laboratory shorthand is not read as
            # our own wording
            parts.append("laboratory note: ``"
                         + latex_escape(analysis_note.strip()) + "''")
        notes.append(parts)
    df["notes"] = notes

    df = df.sort_values(
        ["group_order", "sample_order", "grain_order", "grain", "report",
         "csv_order"],
        na_position="last",
    ).reset_index(drop=True)

    # label repeat U-Th determinations (a), (b) in the order in which they
    # appear in the archived csv file
    letters = df.groupby("Unnamed: 0").cumcount().map(
        lambda i: chr(ord("a") + i))
    df["repeat_label"] = np.where(df["repeat_UTh"], letters, "")
    return df


def table_rows(df, mineral, n_columns):
    """Build the LaTeX body rows of a data table.

    Notes are carried as superscript footnote symbols in the last column so
    that no row has to wrap; the symbols are assigned in the order in which the
    notes first appear and are defined in a legend beneath the table.

    Args:
        df (pandas.DataFrame): Prepared analyses from :func:`prepare`.
        mineral (str): ``'zircon'`` or ``'apatite'``; apatite tables carry Sm
            columns and zircon tables carry a report-file column.
        n_columns (int): Number of columns, used for the group header rows.

    Returns:
        tuple: (list of LaTeX source lines for the body of the longtable,
        list of (symbol, note text) pairs in the order the symbols were
        assigned).
    """
    is_apatite = mineral == "apatite"
    assigned = {}

    def symbol_for(note):
        """Return the footnote symbol for a note, assigning one if needed."""
        if note not in assigned:
            if len(assigned) >= len(SYMBOLS):
                raise ValueError("ran out of footnote symbols")
            assigned[note] = SYMBOLS[len(assigned)]
        return assigned[note]

    lines = []
    previous_group = None
    previous_sample = None
    for _, row in df.iterrows():
        if row["group"] != previous_group:
            if previous_group is not None:
                lines.append(r"\addlinespace[2pt]")
                lines.append(r"\midrule[0.1pt]")
            lines.append(
                f"\\multicolumn{{{n_columns}}}{{@{{}}l}}{{\\textit{{{row['group']}}}}}"
                r" \tabularnewline"
            )
            previous_group = row["group"]
            previous_sample = None
        elif row["sample"] != previous_sample and previous_sample is not None:
            lines.append(r"\addlinespace[2pt]")
        previous_sample = row["sample"]

        grain = row["grain"]
        if row["repeat_UTh"]:
            grain += f" ({row['repeat_label']})"

        cells = [row["sample"], grain]
        if not is_apatite:
            cells.append("" if pd.isna(row["report"]) else f"{int(row['report'])}")
        cells += [
            sig(row["mass_morph_g"] * 1e6, 3) if np.isfinite(row["mass_morph_g"]) else EMDASH,
            fmt_fixed(row["radius"], 1),
            fmt_fixed(row["Ft"], 3),
            fmt_amount(row["He_pmol"], row["He_pmol_1s"]),
            fmt_amount(row["U_ng"], row["U_ng_1s"]),
            fmt_amount(row["Th_ng"], row["Th_ng_1s"]),
        ]
        if is_apatite:
            cells.append(fmt_amount(row["Sm_ng"], row["Sm_ng_1s"]))
        cells += [
            sig(row["U_ppm"], 4) if row["U_ppm"] else EMDASH,
            sig(row["Th_ppm"], 4) if row["Th_ppm"] else EMDASH,
            sig(row["eU"], 4) if row["eU"] else EMDASH,
            fmt_date(row["uncorrected_date"], row["uncorrected_1s"]),
            fmt_date(row["corrected_date"], row["corrected_1s"]),
            (r"\textsuperscript{"
             + ",".join(symbol_for(n) for n in row["notes"]) + "}"
             if row["notes"] else ""),
        ]
        lines.append(" & ".join(cells) + r" \tabularnewline")
    legend = [(symbol, note) for note, symbol in assigned.items()]
    return lines, legend


def standard_summary(df, sample_code, label):
    """Summarize the secondary standard analyses for the table caption.

    Args:
        df (pandas.DataFrame): Prepared analyses.
        sample_code (str): Sample code of the standard (``FCT`` or ``DUR``).
        label (str): Name of the standard as it should read in the caption.

    Returns:
        str: Sentence giving n, the mean date and the standard deviation.
    """
    dates = df[df["sample_code"] == sample_code]["corrected_date"]
    if dates.empty:
        return ""
    return (f"The {label} aliquots reported here give a mean date of "
            f"{dates.mean():.1f}~Ma with a standard deviation of "
            f"{dates.std():.1f}~Ma ($n$~=~{len(dates)}).")


def write_table(df, mineral, path, table_number):
    """Write a longtable LaTeX file for one mineral.

    Args:
        df (pandas.DataFrame): Prepared analyses from :func:`prepare`.
        mineral (str): ``'zircon'`` or ``'apatite'``.
        path (Path): Output ``.tex`` path.
        table_number (str): Supporting information table number, e.g. ``S2``.
    """
    is_apatite = mineral == "apatite"
    source_csv = f"data/UTh\\_He\\_{mineral}\\_data.csv"

    if is_apatite:
        column_spec = "@{}l l r r r r r r r r r r r r c@{}"
        header = [
            r"\textbf{Sample}", r"\textbf{Grain}",
            r"\shortstack{\textbf{Mass}\\\textbf{($\mu$g)}}",
            r"\shortstack{\textbf{R}\\\textbf{($\mu$m)}}",
            r"\shortstack{\textbf{$F_T$}\\{}}",
            r"\shortstack{\textbf{He}\\\textbf{(pmol)}}",
            r"\shortstack{\textbf{U}\\\textbf{(ng)}}",
            r"\shortstack{\textbf{Th}\\\textbf{(ng)}}",
            r"\shortstack{\textbf{Sm}\\\textbf{(ng)}}",
            r"\shortstack{\textbf{U}\\\textbf{(ppm)}}",
            r"\shortstack{\textbf{Th}\\\textbf{(ppm)}}",
            r"\shortstack{\textbf{eU}\\\textbf{(ppm)}}",
            r"\shortstack{\textbf{Raw date}\\\textbf{(Ma)}}",
            r"\shortstack{\textbf{Corrected date}\\\textbf{(Ma)}}",
            r"\textbf{Notes}",
        ]
    else:
        column_spec = "@{}l l c r r r r r r r r r r r c@{}"
        header = [
            r"\textbf{Sample}", r"\textbf{Grain}", r"\textbf{Rpt.}",
            r"\shortstack{\textbf{Mass}\\\textbf{($\mu$g)}}",
            r"\shortstack{\textbf{R}\\\textbf{($\mu$m)}}",
            r"\shortstack{\textbf{$F_T$}\\{}}",
            r"\shortstack{\textbf{He}\\\textbf{(pmol)}}",
            r"\shortstack{\textbf{U}\\\textbf{(ng)}}",
            r"\shortstack{\textbf{Th}\\\textbf{(ng)}}",
            r"\shortstack{\textbf{U}\\\textbf{(ppm)}}",
            r"\shortstack{\textbf{Th}\\\textbf{(ppm)}}",
            r"\shortstack{\textbf{eU}\\\textbf{(ppm)}}",
            r"\shortstack{\textbf{Raw date}\\\textbf{(Ma)}}",
            r"\shortstack{\textbf{Corrected date}\\\textbf{(Ma)}}",
            r"\textbf{Notes}",
        ]
    n_columns = len(header)
    header_row = " & ".join(header) + r" \tabularnewline"

    mineral_label = "apatite" if is_apatite else "zircon"
    standards = (standard_summary(df, "DUR", "Durango apatite") if is_apatite
                 else standard_summary(df, "FCT", "Fish Canyon Tuff zircon"))

    caption_parts = [
        f"Single-grain {mineral_label} (U--Th)/He data for the Montreal River "
        f"monocline, listed by sample group.",
        f"The reduced data are archived as {source_csv} in the project "
        r"repository; the underlying laboratory report files are archived in "
        r"the \texttt{projects\_2023/TA23\_MRM\_MCR} directory of the "
        r"\texttt{wrguenthner/HAL\_data} repository.",
    ]
    if not is_apatite:
        caption_parts.append(
            r"``Rpt.'' identifies the laboratory report file "
            r"(\texttt{TA23\_MRM\_MCR\_$n$\_report.csv}) that an analysis is "
            r"taken from.")
    caption_parts += [
        r"Mass is the aliquot mass from the measured grain dimensions, $R$ is "
        r"the grain radius used in the $\alpha$-ejection correction and in the "
        r"diffusion calculations, and $F_T$ is the $\alpha$-ejection "
        r"correction.",
        (r"He, U, Th, and Sm amounts" if is_apatite
         else r"He, U, and Th amounts")
        + r" are followed in parentheses by their 1$\sigma$ analytical "
        r"uncertainties expressed as a percentage of the amount.",
        r"Concentrations are calculated from these amounts and the aliquot "
        r"mass, and eU~=~[U]~+~0.235[Th].",
        r"``Raw date'' is the date before the $\alpha$-ejection correction and "
        r"``corrected date'' is the date after it; both are given with the "
        r"1$\sigma$ uncertainties returned by the laboratory data reduction, "
        r"which propagate analytical uncertainties only and do not include "
        r"uncertainty in $F_T$ or dispersion between grains.",
    ]
    if is_apatite:
        caption_parts.append(
            r"The Durango apatite aliquots were not measured for grain "
            r"dimensions, so mass, $R$, and concentrations are not reported "
            r"for them and $F_T$ is set to unity.")
    else:
        caption_parts.append(
            r"Sm was not determined for the zircon aliquots.")
        repeats = df[df["repeat_UTh"]]["Unnamed: 0"].nunique()
        caption_parts.append(
            f"{repeats} grains have two U--Th determinations reported for the "
            r"same He aliquot: the He measurement, mass, radius, and $F_T$ of "
            r"the two rows are identical and only the U--Th determination, "
            r"and hence the date, differs. Both determinations are listed, "
            r"labeled (a) and (b) in the order in which they appear in the "
            r"archived data file, and both are carried through the analyses "
            r"in the project repository.")
        caption_parts.append(
            r"Grains flagged as excluded are the analyses filtered out as "
            r"spurious in the project repository.")
    caption_parts.append(
        r"Superscript symbols in the ``Notes'' column are defined beneath the "
        r"table.")
    if standards:
        caption_parts.append(standards)
    caption = " ".join(caption_parts)

    body, legend = table_rows(df, mineral, n_columns)
    legend_rows = []
    for index, (symbol, note) in enumerate(legend):
        if index == 0:
            legend_rows.append(r"\addlinespace[3pt]")
        legend_rows.append(
            f"\\multicolumn{{{n_columns}}}{{@{{}}l}}{{\\textsuperscript{{{symbol}}}"
            f"~{note}}} \\tabularnewline")

    lines = [
        "% Auto-generated by code/make_UTh_He_tables.py",
        f"% Source: data/UTh_He_{mineral}_data.csv",
        "% Do not edit by hand -- rerun the script instead.",
        "%",
        "% Requires: longtable, booktabs, pdflscape (for the landscape environment).",
        f"% Numbered as table {table_number} in the supporting information; in the AGU",
        "% template set the number with \\settablenum or the table counter before",
        "% including this file.",
        "% The table is wide; it is set in landscape here and the landscape",
        "% environment can be removed if the supporting information is",
        "% typeset differently.",
        r"\begin{landscape}",
        r"\begingroup",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{2.5pt}",
        r"\setlength{\LTcapwidth}{\textwidth}",
        r"\renewcommand{\arraystretch}{1.15}",
        f"\\begin{{longtable}}{{{column_spec}}}",
        f"\\caption{{{caption}}}\\label{{tab:{mineral}_UThHe}} \\\\",
        r"\toprule",
        header_row,
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        header_row,
        r"\midrule",
        r"\endhead",
        r"\midrule",
        f"\\multicolumn{{{n_columns}}}{{r}}{{\\textit{{continued on next page}}}}"
        r" \tabularnewline",
        r"\endfoot",
        r"\bottomrule",
    ]
    lines += legend_rows
    lines += [
        r"\endlastfoot",
    ]
    lines += body
    lines += [
        r"\end{longtable}",
        r"\endgroup",
        r"\end{landscape}",
        "",
    ]
    path.write_text("\n".join(lines))
    print(f"wrote {path.relative_to(REPO)}: {len(df)} analyses, "
          f"table {table_number}")


def main():
    inversion = pd.read_csv(INVERSION_CSV)
    inversion_names = set(inversion["grain name"])

    zircon = prepare(ZIRCON_CSV, ZIRCON_GROUPS, inversion_names)
    apatite = prepare(APATITE_CSV, APATITE_GROUPS, inversion_names)

    write_table(zircon, "zircon", ZIRCON_TEX, "S2")
    write_table(apatite, "apatite", APATITE_TEX, "S3")

    for name, df in (("zircon", zircon), ("apatite", apatite)):
        print(f"\n{name} analyses by group:")
        for group, count in df.groupby("group", sort=False).size().items():
            print(f"  {group}: {count}")
        repeats = df[df["repeat_UTh"]]["Unnamed: 0"].nunique()
        if repeats:
            print(f"  grains with repeat U-Th determinations: {repeats}")
        flagged = df["notes"].map(
            lambda notes: any(n.startswith("excluded") or n == NOT_MODELED_NOTE
                              for n in notes))
        if flagged.sum():
            print(f"  analyses flagged as excluded or not modeled: "
                  f"{flagged.sum()}")


if __name__ == "__main__":
    main()
