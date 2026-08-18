#!/usr/bin/env python3
"""Summarize and visualize the completed sample-level FAPROTAX LOSO analysis."""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont


STUDY_ORDER = [
    "PRJEB110492",
    "PRJEB35843",
    "PRJEB98254",
    "PRJNA1156347",
    "PRJNA1253324",
    "PRJNA577971_PRJNA578171",
    "PRJNA835240",
]

STUDY_LABELS = {
    "PRJEB110492": "PRJEB110492",
    "PRJEB35843": "PRJEB35843",
    "PRJEB98254": "PRJEB98254",
    "PRJNA1156347": "PRJNA1156347",
    "PRJNA1253324": "PRJNA1253324",
    "PRJNA577971_PRJNA578171": "PRJNA577971+\nPRJNA578171",
    "PRJNA835240": "PRJNA835240",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, type=Path)
    parser.add_argument("--output-dir", type=Path)
    return parser.parse_args()


def read_tsv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    return pd.read_csv(path, sep="\t")


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, lineterminator="\n")


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont:
    candidates = [
        Path("C:/Windows/Fonts/arialbd.ttf" if bold else "C:/Windows/Fonts/arial.ttf"),
        Path("C:/Windows/Fonts/calibrib.ttf" if bold else "C:/Windows/Fonts/calibri.ttf"),
    ]
    for candidate in candidates:
        if candidate.exists():
            return ImageFont.truetype(str(candidate), size=size)
    return ImageFont.load_default()


def quantile(values: np.ndarray, q: float) -> float:
    return float(np.quantile(values.astype(float), q))


def coverage_summary(coverage: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for study in STUDY_ORDER:
        values = coverage.loc[
            coverage["study"].eq(study), "assigned_read_fraction"
        ].to_numpy(float)
        rows.append(
            {
                "study": study,
                "n_samples": len(values),
                "mean_assigned_read_fraction": float(values.mean()),
                "median_assigned_read_fraction": float(np.median(values)),
                "q1_assigned_read_fraction": quantile(values, 0.25),
                "q3_assigned_read_fraction": quantile(values, 0.75),
                "minimum_assigned_read_fraction": float(values.min()),
                "maximum_assigned_read_fraction": float(values.max()),
            }
        )
    return pd.DataFrame(rows)


def function_prevalence(
    abundance: pd.DataFrame, stable_functions: list[str]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    function_columns = [
        column for column in abundance.columns if column not in {"sample_uid", "study"}
    ]
    rows = []
    stable_abundance_rows = []
    for study in STUDY_ORDER:
        cohort = abundance.loc[abundance["study"].eq(study)]
        for function_name in function_columns:
            values = cohort[function_name].to_numpy(float)
            rows.append(
                {
                    "study": study,
                    "function": function_name,
                    "prevalence": float((values > 0).mean()),
                }
            )
            if function_name in stable_functions:
                stable_abundance_rows.append(
                    {
                        "study": study,
                        "function": function_name,
                        "n_samples": len(values),
                        "prevalence": float((values > 0).mean()),
                        "median_fraction_of_16S": float(np.median(values)),
                        "q1_fraction_of_16S": quantile(values, 0.25),
                        "q3_fraction_of_16S": quantile(values, 0.75),
                    }
                )
    return pd.DataFrame(rows), pd.DataFrame(stable_abundance_rows)


def make_figure(
    output_path_png: Path,
    output_path_tiff: Path,
    coverage: pd.DataFrame,
    stable: pd.DataFrame,
    prevalence: pd.DataFrame,
    folds: pd.DataFrame,
) -> None:
    width, height = 5000, 3200
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)
    title_font = font(66, bold=True)
    panel_font = font(58, bold=True)
    axis_font = font(34)
    small_font = font(29)
    tiny_font = font(25)

    draw.text(
        (150, 70),
        "Sample-level FAPROTAX annotation and seven-cohort LOSO validation",
        fill="#111111",
        font=title_font,
    )

    # Panel A: annotation coverage by cohort.
    ax_a = (170, 310, 1400, 2400)
    draw.text((ax_a[0], 215), "A", fill="#111111", font=panel_font)
    draw.text(
        (ax_a[0] + 95, 225),
        "FAPROTAX annotation coverage",
        fill="#111111",
        font=font(44, bold=True),
    )
    left, top, right, bottom = ax_a
    ymax = 0.60
    for tick in np.arange(0, ymax + 0.001, 0.1):
        y = bottom - (tick / ymax) * (bottom - top)
        draw.line((left, y, right, y), fill="#DDDDDD", width=3)
        draw.text((left - 120, y - 18), f"{tick:.1f}", fill="#333333", font=axis_font)
    draw.line((left, top, left, bottom), fill="#222222", width=5)
    draw.line((left, bottom, right, bottom), fill="#222222", width=5)
    box_width = 78
    x_positions = np.linspace(left + 110, right - 80, len(STUDY_ORDER))
    rng = np.random.default_rng(20260729)
    for x, study in zip(x_positions, STUDY_ORDER):
        vals = coverage.loc[
            coverage["study"].eq(study), "assigned_read_fraction"
        ].to_numpy(float)
        q1, med, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
        vmin, vmax = vals.min(), vals.max()
        scale_y = lambda value: bottom - (value / ymax) * (bottom - top)
        draw.line((x, scale_y(vmin), x, scale_y(vmax)), fill="#333333", width=4)
        draw.rectangle(
            (x - box_width / 2, scale_y(q3), x + box_width / 2, scale_y(q1)),
            fill="#73A9C2",
            outline="#1E526A",
            width=4,
        )
        draw.line(
            (x - box_width / 2, scale_y(med), x + box_width / 2, scale_y(med)),
            fill="#111111",
            width=6,
        )
        for value in vals:
            jitter = float(rng.uniform(-box_width * 0.35, box_width * 0.35))
            yy = scale_y(value)
            draw.ellipse(
                (x + jitter - 5, yy - 5, x + jitter + 5, yy + 5),
                fill="#1E526A",
            )
        label = STUDY_LABELS[study]
        label_img = Image.new("RGBA", (360, 180), (255, 255, 255, 0))
        label_draw = ImageDraw.Draw(label_img)
        label_draw.multiline_text(
            (0, 0), label, fill="#222222", font=tiny_font, spacing=2, align="center"
        )
        label_img = label_img.rotate(55, expand=True, resample=Image.Resampling.BICUBIC)
        image.alpha_composite(label_img, (int(x - 35), bottom + 20)) if image.mode == "RGBA" else image.paste(
            label_img, (int(x - 35), bottom + 20), label_img
        )
    y_label = Image.new("RGBA", (650, 70), (255, 255, 255, 0))
    ImageDraw.Draw(y_label).text(
        (0, 0), "Assigned fraction of original 16S reads", fill="#222222", font=axis_font
    )
    y_label = y_label.rotate(90, expand=True)
    image.paste(y_label, (18, 810), y_label)

    # Panel B: prevalence heatmap for functions stable in all seven folds.
    heat_left, heat_top, heat_right, heat_bottom = 2300, 310, 3650, 2400
    draw.text((1510, 215), "B", fill="#111111", font=panel_font)
    draw.text(
        (1605, 225),
        "Held-out prevalence of LOSO-stable functions",
        fill="#111111",
        font=font(44, bold=True),
    )
    stable_names = stable.loc[stable["stable_all_heldouts"], "function"].tolist()
    cell_w = (heat_right - heat_left) / len(STUDY_ORDER)
    cell_h = (heat_bottom - heat_top) / len(stable_names)
    matrix = prevalence.pivot(index="function", columns="study", values="prevalence")

    def heat_color(value: float) -> tuple[int, int, int]:
        low = np.array([244, 239, 226])
        high = np.array([24, 118, 120])
        blend = np.clip((value - 0.45) / 0.55, 0, 1)
        return tuple((low * (1 - blend) + high * blend).astype(int))

    for row_idx, function_name in enumerate(stable_names):
        y0 = heat_top + row_idx * cell_h
        y1 = heat_top + (row_idx + 1) * cell_h
        draw.text(
            (heat_left - 770, y0 + cell_h * 0.28),
            function_name,
            fill="#222222",
            font=small_font,
        )
        for col_idx, study in enumerate(STUDY_ORDER):
            value = float(matrix.loc[function_name, study])
            x0 = heat_left + col_idx * cell_w
            x1 = heat_left + (col_idx + 1) * cell_w
            draw.rectangle((x0, y0, x1, y1), fill=heat_color(value), outline="white", width=3)
            text_fill = "white" if value >= 0.76 else "#222222"
            value_text = f"{value:.2f}"
            bbox = draw.textbbox((0, 0), value_text, font=tiny_font)
            draw.text(
                (
                    (x0 + x1 - (bbox[2] - bbox[0])) / 2,
                    (y0 + y1 - (bbox[3] - bbox[1])) / 2 - 3,
                ),
                value_text,
                fill=text_fill,
                font=tiny_font,
            )
    for col_idx, study in enumerate(STUDY_ORDER):
        x = heat_left + (col_idx + 0.5) * cell_w
        label_img = Image.new("RGBA", (350, 180), (255, 255, 255, 0))
        ImageDraw.Draw(label_img).multiline_text(
            (0, 0),
            STUDY_LABELS[study],
            fill="#222222",
            font=tiny_font,
            spacing=2,
            align="center",
        )
        label_img = label_img.rotate(55, expand=True, resample=Image.Resampling.BICUBIC)
        image.paste(label_img, (int(x - 30), heat_bottom + 20), label_img)

    # Panel C: fold validation rates against the matched null.
    c_left, c_top, c_right, c_bottom = 4000, 310, 4850, 2400
    draw.text((3830, 215), "C", fill="#111111", font=panel_font)
    draw.text(
        (3925, 225),
        "LOSO validation",
        fill="#111111",
        font=font(44, bold=True),
    )
    for tick in np.arange(0, 1.01, 0.2):
        y = c_bottom - tick * (c_bottom - c_top)
        draw.line((c_left, y, c_right, y), fill="#DDDDDD", width=3)
        draw.text((c_left - 95, y - 18), f"{tick:.1f}", fill="#333333", font=axis_font)
    draw.line((c_left, c_top, c_left, c_bottom), fill="#222222", width=5)
    draw.line((c_left, c_bottom, c_right, c_bottom), fill="#222222", width=5)
    c_x = np.linspace(c_left + 70, c_right - 55, len(STUDY_ORDER))
    bar_w = 35
    folds_idx = folds.set_index("heldout_study")
    for x, study in zip(c_x, STUDY_ORDER):
        obs = float(folds_idx.loc[study, "validation_rate"])
        null_mean = float(folds_idx.loc[study, "matched_null_mean"])
        obs_y = c_bottom - obs * (c_bottom - c_top)
        null_y = c_bottom - null_mean * (c_bottom - c_top)
        draw.rectangle(
            (x - bar_w - 3, obs_y, x - 3, c_bottom),
            fill="#D45D3A",
            outline="#8C2F18",
            width=3,
        )
        draw.rectangle(
            (x + 3, null_y, x + bar_w + 3, c_bottom),
            fill="#9C9C9C",
            outline="#555555",
            width=3,
        )
        label_img = Image.new("RGBA", (350, 180), (255, 255, 255, 0))
        ImageDraw.Draw(label_img).multiline_text(
            (0, 0),
            STUDY_LABELS[study],
            fill="#222222",
            font=tiny_font,
            spacing=2,
            align="center",
        )
        label_img = label_img.rotate(55, expand=True, resample=Image.Resampling.BICUBIC)
        image.paste(label_img, (int(x - 28), c_bottom + 20), label_img)
    draw.rectangle((4020, 2860, 4070, 2910), fill="#D45D3A")
    draw.text((4090, 2860), "Observed", fill="#222222", font=small_font)
    draw.rectangle((4350, 2860, 4400, 2910), fill="#9C9C9C")
    draw.text((4420, 2860), "Matched-null mean", fill="#222222", font=small_font)
    draw.text(
        (3820, 2990),
        "Global macro-average: observed 0.944; null 0.910; empirical P = 0.001",
        fill="#222222",
        font=tiny_font,
    )

    image.save(output_path_png, dpi=(300, 300), optimize=True)
    image.save(output_path_tiff, dpi=(300, 300), compression="tiff_lzw")


def main() -> None:
    args = parse_args()
    input_dir = args.input_dir.resolve()
    output_dir = (args.output_dir or input_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    coverage = read_tsv(input_dir / "sample_annotation_coverage.tsv")
    abundance = read_tsv(input_dir / "sample_function_abundance_fraction_of_16S.tsv")
    folds = read_tsv(input_dir / "faprotax_loso_fold_summary.tsv")
    details = read_tsv(input_dir / "faprotax_loso_function_details.tsv")
    stable = read_tsv(input_dir / "faprotax_stable_functions.tsv")
    null = read_tsv(input_dir / "faprotax_loso_prevalence_matched_null.tsv")

    stable["stable_all_heldouts"] = stable["stable_all_heldouts"].astype(str).eq("True")
    stable_functions = stable.loc[stable["stable_all_heldouts"], "function"].tolist()

    if coverage["sample_uid"].nunique() != 265:
        raise ValueError("Expected 265 unique samples")
    if set(coverage["study"]) != set(STUDY_ORDER):
        raise ValueError("Unexpected study identifiers")
    if abundance["sample_uid"].tolist() != coverage["sample_uid"].tolist():
        raise ValueError("Sample order differs between coverage and abundance outputs")
    if len(stable_functions) != 12:
        raise ValueError(f"Expected 12 all-heldout stable functions, found {len(stable_functions)}")
    function_columns = [
        column for column in abundance.columns if column not in {"sample_uid", "study"}
    ]
    if len(function_columns) != 44:
        raise ValueError(f"Expected 44 FAPROTAX functions, found {len(function_columns)}")
    if abundance[function_columns].lt(0).any().any():
        raise ValueError("Negative function abundance detected")

    coverage_by_study = coverage_summary(coverage)
    prevalence, stable_abundance = function_prevalence(abundance, stable_functions)
    macro_observed = float(folds["validation_rate"].mean())
    micro_observed = float(
        folds["validated_core_functions"].sum()
        / folds["discovery_core_functions"].sum()
    )
    null_macro = null.groupby("permutation", as_index=False)["null_validation_rate"].mean()
    global_p = float(
        (1 + (null_macro["null_validation_rate"] >= macro_observed).sum())
        / (len(null_macro) + 1)
    )
    global_summary = pd.DataFrame(
        [
            {
                "samples": coverage["sample_uid"].nunique(),
                "studies": coverage["study"].nunique(),
                "functions": len(function_columns),
                "stable_all_heldouts": len(stable_functions),
                "candidate_fold_predictions": int(folds["discovery_core_functions"].sum()),
                "validated_fold_predictions": int(folds["validated_core_functions"].sum()),
                "micro_validation_rate": micro_observed,
                "macro_validation_rate": macro_observed,
                "matched_null_macro_mean": float(null_macro["null_validation_rate"].mean()),
                "matched_null_macro_q025": float(
                    null_macro["null_validation_rate"].quantile(0.025)
                ),
                "matched_null_macro_q975": float(
                    null_macro["null_validation_rate"].quantile(0.975)
                ),
                "matched_null_global_empirical_p": global_p,
                "mean_annotation_coverage": float(
                    coverage["assigned_read_fraction"].mean()
                ),
                "median_annotation_coverage": float(
                    coverage["assigned_read_fraction"].median()
                ),
                "minimum_annotation_coverage": float(
                    coverage["assigned_read_fraction"].min()
                ),
                "maximum_annotation_coverage": float(
                    coverage["assigned_read_fraction"].max()
                ),
            }
        ]
    )

    write_tsv(coverage_by_study, output_dir / "faprotax_annotation_coverage_by_study.tsv")
    write_tsv(prevalence, output_dir / "faprotax_function_prevalence_by_study.tsv")
    write_tsv(
        stable_abundance,
        output_dir / "faprotax_stable_function_abundance_by_study.tsv",
    )
    write_tsv(global_summary, output_dir / "faprotax_loso_global_summary.tsv")

    make_figure(
        output_dir / "Supplementary_Figure_FAPROTAX_LOSO.png",
        output_dir / "Supplementary_Figure_FAPROTAX_LOSO.tiff",
        coverage,
        stable,
        prevalence,
        folds,
    )

    failed = details.loc[details["validated"].astype(str).eq("False")].copy()
    stable_lines = "\n".join(
        f"- `{row.function}`: minimum held-out prevalence "
        f"{float(row.minimum_heldout_prevalence):.3f}; mean "
        f"{float(row.mean_heldout_prevalence):.3f}."
        for row in stable.loc[stable["stable_all_heldouts"]].itertuples()
    )
    failed_lines = "\n".join(
        f"- {row.heldout_study}: `{row.function}` "
        f"(held-out prevalence {float(row.heldout_prevalence):.3f})."
        for row in failed.itertuples()
    )
    report = f"""# Sample-level FAPROTAX and seven-cohort LOSO validation

## Analysis status

The analysis used the official FAPROTAX 1.2.12 database without manual additions.
The input comprised 265 bacterial 16S samples from seven independent studies and
744 genus-level taxonomic rows. FAPROTAX returned 44 represented functional
groups.

## Strict LOSO definition

For each held-out study, a candidate function was discovered using only the six
training studies. A function qualified when it occurred in at least 50% of
samples within every training study. It validated when its prevalence was also
at least 50% in the untouched held-out study. The final all-heldout set contains
only functions that qualified and validated in all seven folds.

## Main results

- Twelve functions qualified and validated in all seven held-out studies.
- Across the seven folds, 84 of 90 candidate-function predictions validated
  (micro-average {micro_observed:.3f}).
- The mean fold validation rate was {macro_observed:.3f}.
- Against 999 training-prevalence-decile-matched null draws, the global null
  mean was {float(null_macro["null_validation_rate"].mean()):.3f}, with an
  empirical global P value of {global_p:.3f}.
- The mean fraction of original 16S reads assigned to at least one FAPROTAX
  function was {float(coverage["assigned_read_fraction"].mean()):.3f}; the
  sample range was {float(coverage["assigned_read_fraction"].min()):.3f} to
  {float(coverage["assigned_read_fraction"].max()):.3f}.

## Functions stable in every held-out study

{stable_lines}

The 12 labels are not 12 statistically independent processes. FAPROTAX contains
nested terms, including aerobic chemoheterotrophy within chemoheterotrophy,
aerobic ammonia oxidation within nitrification, and photoheterotrophy within
phototrophy. They are therefore best described as 12 database labels spanning
approximately nine nonredundant functional families.

## Fold-specific failures

{failed_lines}

The combined PRJNA577971/PRJNA578171 cohort had the lowest annotation coverage
and accounted for four of the six failed fold-level predictions. PRJNA835240
accounted for the remaining two.

## Interpretation for the manuscript

These data support cross-study recurrence of FAPROTAX-predicted ecological
functions, not direct measurement of genes, transcripts, metabolites, or
biochemical rates. Because mean annotation coverage was only about 27%, the
result should be presented as a complementary, hypothesis-generating functional
layer. It should not replace the paper's primary contribution: the
study-aware taxonomic transferability framework.

## Manuscript-ready Methods text

“We inferred sample-level ecological functions from the harmonized bacterial
genus count matrix using the unmodified FAPROTAX v1.2.12 database. Counts
assigned to each represented functional group were retained both as counts and
as fractions of the original 16S library; annotation coverage was calculated
from the FAPROTAX unassigned-taxa output. Functional stability was evaluated by
seven-fold leave-one-study-out validation. In each fold, a function was
discovered solely from the six training studies if its within-study prevalence
was at least 50% in every training study, and it was considered validated if its
prevalence was at least 50% in the untouched held-out study. We additionally
compared the observed fold validation rate with 999 random function sets matched
to the training prevalence deciles of the discovered functions (seed
20260728).”

## Manuscript-ready Results text

“Official FAPROTAX annotation yielded 44 represented functional groups across
265 samples. In strict seven-fold leave-one-study-out validation, 84 of 90
candidate-function predictions replicated in the held-out cohort (93.3%), and
12 FAPROTAX labels met the prevalence criterion in every training definition
and every corresponding held-out study. The mean fold validation rate was
94.4%, exceeding the prevalence-matched null expectation of 91.0% (999 draws;
global empirical P=0.001). These labels included broad or nested categories
related to chemoheterotrophy, cellulolysis, nitrogen fixation, nitrate
reduction, photoheterotrophy/phototrophy, nitrification/ammonia oxidation,
fermentation, iron respiration, and predatory or exoparasitic lifestyles.
However, only a mean 27.2% of original 16S reads received at least one FAPROTAX
assignment, so these results are interpreted as taxonomically inferred
ecological potential rather than direct evidence of metabolic activity.”

## Recommended placement

Use the three-panel figure and detailed function tables as Supplementary
Information. In the main text, retain one short paragraph reporting the strict
LOSO result and the coverage limitation. Avoid the phrases “conserved metabolic
core” and “confirmed functional activity”; prefer “cross-study stable
FAPROTAX-predicted functions” or “inferred functional potential.”
"""
    (output_dir / "FAPROTAX_LOSO_RESULTS_AND_INTERPRETATION.md").write_text(
        report, encoding="utf-8", newline="\n"
    )

    print(global_summary.to_json(orient="records", indent=2))


if __name__ == "__main__":
    main()
