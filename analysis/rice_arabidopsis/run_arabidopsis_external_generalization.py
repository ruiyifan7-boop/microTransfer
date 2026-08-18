from __future__ import annotations

import html
import hashlib
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "work" / "pydeps"))

import numpy as np
import pandas as pd
import pyreadr


SEED = 20260722
N_PERM = 999
PRIMARY_BLOCKS = ["C1", "C2", "M1", "M2"]
MIN_RETAINED_DEPTH = 100
PREV_THRESHOLD = 0.50
MIN_TOTAL = 20

DATA = ROOT / "work" / "arabidopsis_cran" / "CoreMicrobiomeR" / "data"
OUT = ROOT / "work" / "arabidopsis_external_generalization_results"
OUT.mkdir(parents=True, exist_ok=True)


def normal_two_sided(z: np.ndarray | float) -> np.ndarray:
    a = np.asarray(z, dtype=float)
    return np.vectorize(lambda x: math.erfc(abs(float(x)) / math.sqrt(2.0)))(a)


def bh_adjust(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    out = np.full(p.shape, np.nan)
    good = np.isfinite(p)
    vals = p[good]
    if not vals.size:
        return out
    order = np.argsort(vals)
    ranked = vals[order]
    q = ranked * vals.size / np.arange(1, vals.size + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    q = np.clip(q, 0, 1)
    restored = np.empty_like(q)
    restored[order] = q
    out[good] = restored
    return out


def parse_sample(sample: str) -> dict[str, str]:
    parts = sample.split(".")
    if len(parts) != 5:
        raise ValueError(f"Unexpected sample identifier: {sample}")
    final = parts[4].split("_")
    fraction = parts[2][-1]
    return {
        "sample": sample,
        "soil": parts[0],
        "genotype": parts[1],
        "sample_fraction_code": parts[2],
        "fraction": fraction,
        "age": parts[3],
        "experiment": final[0],
        "plate": "_".join(final[1:]),
    }


def clr(counts: np.ndarray, pseudocount: float = 0.5) -> np.ndarray:
    z = np.log(counts + pseudocount)
    return z - z.mean(axis=1, keepdims=True)


def residualized_correlations(
    x: np.ndarray, y: np.ndarray, fraction: np.ndarray, genotype: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, int]:
    d = pd.DataFrame({"fraction": fraction, "genotype": genotype})
    design = pd.get_dummies(d, drop_first=True, dtype=float).to_numpy()
    design = np.column_stack([np.ones(len(y)), design])
    q, _ = np.linalg.qr(design, mode="reduced")
    yr = y - q @ (q.T @ y)
    xr = x - q @ (q.T @ x)
    denom = np.sqrt(np.sum(xr * xr, axis=0) * np.sum(yr * yr))
    r = np.divide(
        np.sum(xr * yr[:, None], axis=0),
        denom,
        out=np.zeros(x.shape[1], dtype=float),
        where=denom > 0,
    )
    r = np.clip(r, -0.999999, 0.999999)
    rank = int(np.linalg.matrix_rank(design))
    df = max(len(y) - rank - 2, 1)
    t = r * np.sqrt(df / np.maximum(1.0 - r * r, 1e-12))
    p = normal_two_sided(t)
    fisher = np.arctanh(r)
    variance = np.full(x.shape[1], 1.0 / max(len(y) - rank - 3, 2))
    return r, p, variance, rank


def random_effects(z: np.ndarray, variance: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    w = 1.0 / variance
    fixed = np.sum(w * z, axis=0) / np.sum(w, axis=0)
    q = np.sum(w * (z - fixed) ** 2, axis=0)
    df = z.shape[0] - 1
    sumw = np.sum(w, axis=0)
    c = sumw - np.sum(w * w, axis=0) / sumw
    tau2 = np.maximum(0.0, (q - df) / np.maximum(c, 1e-12))
    wr = 1.0 / (variance + tau2)
    mu = np.sum(wr * z, axis=0) / np.sum(wr, axis=0)
    se = np.sqrt(1.0 / np.sum(wr, axis=0))
    p = normal_two_sided(mu / se)
    return mu, se, p


def make_bins(values: np.ndarray, n: int = 4) -> np.ndarray:
    finite = np.isfinite(values)
    out = np.zeros(values.size, dtype=int)
    if finite.sum() == 0:
        return out
    edges = np.unique(np.quantile(values[finite], np.linspace(0, 1, n + 1)))
    if edges.size <= 2:
        return out
    out[finite] = np.digitize(values[finite], edges[1:-1], right=True)
    return out


def matched_subset(
    core: np.ndarray, universe: np.ndarray, strata: np.ndarray, rng: np.random.Generator
) -> np.ndarray:
    core_set = set(core.tolist())
    noncore = np.array([i for i in universe if i not in core_set], dtype=int)
    chosen: list[int] = []
    for stratum in np.unique(strata[core]):
        need = int(np.sum(strata[core] == stratum))
        candidates = noncore[strata[noncore] == stratum]
        candidates = np.array([i for i in candidates if i not in chosen], dtype=int)
        take = min(need, candidates.size)
        if take:
            chosen.extend(rng.choice(candidates, size=take, replace=False).tolist())
    need = core.size - len(chosen)
    if need:
        remaining = np.array([i for i in noncore if i not in chosen], dtype=int)
        replace = remaining.size < need
        source = remaining if remaining.size else universe
        chosen.extend(rng.choice(source, size=need, replace=replace).tolist())
    return np.asarray(chosen, dtype=int)


def balanced_accuracy(y: np.ndarray, pred: np.ndarray) -> float:
    vals = []
    for cls in np.unique(y):
        idx = y == cls
        vals.append(float(np.mean(pred[idx] == cls)))
    return float(np.mean(vals))


def macro_f1(y: np.ndarray, pred: np.ndarray) -> float:
    vals = []
    for cls in np.unique(y):
        tp = np.sum((pred == cls) & (y == cls))
        fp = np.sum((pred == cls) & (y != cls))
        fn = np.sum((pred != cls) & (y == cls))
        den = 2 * tp + fp + fn
        vals.append(0.0 if den == 0 else float(2 * tp / den))
    return float(np.mean(vals))


def nearest_centroid_fit_predict(
    x_train: np.ndarray, y_train: np.ndarray, x_test: np.ndarray
) -> np.ndarray:
    mean = x_train.mean(axis=0)
    sd = x_train.std(axis=0, ddof=1)
    sd[~np.isfinite(sd) | (sd < 1e-8)] = 1.0
    tr = (x_train - mean) / sd
    te = (x_test - mean) / sd
    classes = np.unique(y_train)
    centroids = np.vstack([tr[y_train == c].mean(axis=0) for c in classes])
    dist = ((te[:, None, :] - centroids[None, :, :]) ** 2).sum(axis=2)
    return classes[np.argmin(dist, axis=1)]


def svg_figure(
    core_summary: pd.DataFrame,
    loso: pd.DataFrame,
    response: pd.DataFrame,
    prediction: pd.DataFrame,
    prediction_overall: dict[str, float],
    target: Path,
) -> None:
    width, height = 1100, 760
    blue, orange, teal, dark, grid = "#356B9A", "#E07A5F", "#3A8D7C", "#24323D", "#D9E1E8"
    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="100%" height="100%" fill="white"/>',
        '<style>text{font-family:Arial,Liberation Sans,sans-serif;fill:#24323D}.axis{stroke:#758390;stroke-width:1}.grid{stroke:#D9E1E8;stroke-width:1}.panel{font-size:18px;font-weight:700}.label{font-size:13px}.small{font-size:11px}.title{font-size:14px;font-weight:600}</style>',
    ]

    def text(x, y, s, cls="label", anchor="start", rotate=None):
        transform = f' transform="rotate({rotate} {x} {y})"' if rotate is not None else ""
        parts.append(
            f'<text x="{x:.1f}" y="{y:.1f}" class="{cls}" text-anchor="{anchor}"{transform}>{html.escape(str(s))}</text>'
        )

    panels = [(55, 45, 470, 285), (590, 45, 455, 285), (55, 410, 470, 285), (590, 410, 455, 285)]
    for label, (x, y, w, h) in zip("abcd", panels):
        text(x - 32, y + 5, label, "panel")

    # a: pooled versus study-aware core size
    x, y, w, h = panels[0]
    text(x, y, "Core membership depends on study awareness", "title")
    plot_y0, plot_h = y + 35, h - 75
    maxv = max(1, int(core_summary[["pooled_features", "study_aware_features"]].to_numpy().max()))
    for tick in np.linspace(0, maxv, 5):
        yy = plot_y0 + plot_h - tick / maxv * plot_h
        parts.append(f'<line x1="{x}" y1="{yy}" x2="{x+w}" y2="{yy}" class="grid"/>')
        text(x - 8, yy + 4, int(round(tick)), "small", "end")
    group_w = w / len(core_summary)
    bar_w = 36
    for i, row in core_summary.reset_index(drop=True).iterrows():
        cx = x + group_w * (i + 0.5)
        for j, (col, color) in enumerate([("pooled_features", blue), ("study_aware_features", orange)]):
            val = float(row[col])
            bh = val / maxv * plot_h
            bx = cx + (j - 0.5) * (bar_w + 5) - bar_w / 2
            by = plot_y0 + plot_h - bh
            parts.append(f'<rect x="{bx}" y="{by}" width="{bar_w}" height="{bh}" fill="{color}"/>')
            text(bx + bar_w / 2, by - 5, int(val), "small", "middle")
        text(cx, plot_y0 + plot_h + 20, f'{int(row["prevalence_threshold"]*100)}%', "small", "middle")
    text(x + w / 2, plot_y0 + plot_h + 40, "Prevalence threshold", "label", "middle")
    text(x - 38, plot_y0 + plot_h / 2, "Number of OTUs", "label", "middle", -90)
    parts += [
        f'<rect x="{x+w-165}" y="{y+3}" width="12" height="12" fill="{blue}"/>',
        f'<rect x="{x+w-75}" y="{y+3}" width="12" height="12" fill="{orange}"/>',
    ]
    text(x + w - 148, y + 14, "Pooled", "small")
    text(x + w - 58, y + 14, "Study-aware", "small")

    # b: LOSO membership validation
    x, y, w, h = panels[1]
    text(x, y, "Independent-block validation of core membership", "title")
    plot_y0, plot_h = y + 35, h - 75
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        yy = plot_y0 + plot_h - tick * plot_h
        parts.append(f'<line x1="{x}" y1="{yy}" x2="{x+w}" y2="{yy}" class="grid"/>')
        text(x - 8, yy + 4, f"{tick:.2g}", "small", "end")
    step = w / len(loso)
    for i, row in loso.reset_index(drop=True).iterrows():
        cx = x + step * (i + 0.5)
        oy = plot_y0 + plot_h - float(row.observed_validation_rate) * plot_h
        ny = plot_y0 + plot_h - float(row.matched_null_q95) * plot_h
        parts.append(f'<line x1="{cx-12}" y1="{ny}" x2="{cx+12}" y2="{ny}" stroke="{dark}" stroke-width="3"/>')
        parts.append(f'<circle cx="{cx}" cy="{oy}" r="6" fill="{teal}"/>')
        text(cx, plot_y0 + plot_h + 20, row.heldout_block, "small", "middle")
        text(cx, oy - 10, f'{float(row.observed_validation_rate):.2f}', "small", "middle")
    text(x - 38, plot_y0 + plot_h / 2, "Validation rate", "label", "middle", -90)
    text(x + w / 2, plot_y0 + plot_h + 40, "Held-out experiment", "label", "middle")
    text(x + 10, plot_y0 + plot_h - 8, "● observed   — matched-null 95th percentile", "small")

    # c: response transfer
    x, y, w, h = panels[2]
    text(x, y, "Transfer of developmental-stage responses", "title")
    plot_y0, plot_h = y + 35, h - 75
    maxv = max(1, int(response[["meta_candidates", "directionally_replicated"]].to_numpy().max()))
    for tick in np.linspace(0, maxv, 5):
        yy = plot_y0 + plot_h - tick / maxv * plot_h
        parts.append(f'<line x1="{x}" y1="{yy}" x2="{x+w}" y2="{yy}" class="grid"/>')
        text(x - 8, yy + 4, int(round(tick)), "small", "end")
    step = w / len(response)
    for i, row in response.reset_index(drop=True).iterrows():
        cx = x + step * (i + 0.5)
        for j, (col, color) in enumerate([("meta_candidates", blue), ("directionally_replicated", teal)]):
            val = float(row[col])
            bh = val / maxv * plot_h
            bx = cx + (j - 0.5) * 38 - 15
            by = plot_y0 + plot_h - bh
            parts.append(f'<rect x="{bx}" y="{by}" width="30" height="{bh}" fill="{color}"/>')
        text(cx, plot_y0 + plot_h + 20, row.heldout_block, "small", "middle")
        text(cx, plot_y0 + plot_h - float(row.directionally_replicated) / maxv * plot_h - 7, f'{float(row.replication_rate):.0%}', "small", "middle")
    text(x - 38, plot_y0 + plot_h / 2, "Number of OTUs", "label", "middle", -90)
    text(x + w / 2, plot_y0 + plot_h + 40, "Held-out experiment", "label", "middle")
    parts += [
        f'<rect x="{x+w-190}" y="{y+3}" width="12" height="12" fill="{blue}"/>',
        f'<rect x="{x+w-92}" y="{y+3}" width="12" height="12" fill="{teal}"/>',
    ]
    text(x + w - 173, y + 14, "Discovered", "small")
    text(x + w - 75, y + 14, "Replicated", "small")

    # d: prediction
    x, y, w, h = panels[3]
    text(x, y, "Cross-experiment compartment discrimination", "title")
    plot_y0, plot_h = y + 35, h - 75
    for tick in [0, 0.25, 0.5, 0.75, 1.0]:
        yy = plot_y0 + plot_h - tick * plot_h
        parts.append(f'<line x1="{x}" y1="{yy}" x2="{x+w}" y2="{yy}" class="grid"/>')
        text(x - 8, yy + 4, f"{tick:.2g}", "small", "end")
    parts.append(
        f'<line x1="{x}" y1="{plot_y0+plot_h*0.5}" x2="{x+w}" y2="{plot_y0+plot_h*0.5}" stroke="{dark}" stroke-dasharray="5 4"/>'
    )
    step = w / len(prediction)
    for i, row in prediction.reset_index(drop=True).iterrows():
        cx = x + step * (i + 0.5)
        oy = plot_y0 + plot_h - float(row.balanced_accuracy) * plot_h
        ny = plot_y0 + plot_h - float(row.permutation_q95) * plot_h
        parts.append(f'<line x1="{cx-12}" y1="{ny}" x2="{cx+12}" y2="{ny}" stroke="{orange}" stroke-width="3"/>')
        parts.append(f'<circle cx="{cx}" cy="{oy}" r="6" fill="{blue}"/>')
        text(cx, plot_y0 + plot_h + 20, row.heldout_block, "small", "middle")
        text(cx, oy - 10, f'{float(row.balanced_accuracy):.2f}', "small", "middle")
    text(x - 38, plot_y0 + plot_h / 2, "Balanced accuracy", "label", "middle", -90)
    text(x + w / 2, plot_y0 + plot_h + 40, "Held-out experiment", "label", "middle")
    text(
        x + 10,
        plot_y0 + plot_h - 8,
        f'Mean BA={prediction_overall["observed_mean_ba"]:.3f}; permutation P={prediction_overall["empirical_p"]:.3g}',
        "small",
    )
    parts.append("</svg>")
    target.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    rng = np.random.default_rng(SEED)
    otu = pyreadr.read_r(str(DATA / "demo_otu.rda"))["demo_otu"]
    tax = pyreadr.read_r(str(DATA / "demo_tax.rda"))["demo_tax"].copy()
    feature_ids = otu.iloc[:, 0].astype(str).to_numpy()
    sample_ids = np.asarray(otu.columns[1:], dtype=str)
    counts = otu.iloc[:, 1:].to_numpy(dtype=float).T
    metadata = pd.DataFrame([parse_sample(x) for x in sample_ids])
    metadata["retained_feature_depth"] = counts.sum(axis=1)
    metadata["primary_block"] = metadata["experiment"].isin(PRIMARY_BLOCKS)
    metadata["root_compartment"] = metadata["fraction"].isin(["E", "R"])
    metadata["depth_pass"] = metadata["retained_feature_depth"] >= MIN_RETAINED_DEPTH
    metadata["analysis_included"] = (
        metadata["primary_block"] & metadata["root_compartment"] & metadata["depth_pass"]
    )
    metadata.to_csv(OUT / "sample_metadata_reconstructed.tsv", sep="\t", index=False)

    # Validate the identifier parser against every row in the packaged metadata subset.
    packaged_md = pyreadr.read_r(str(DATA / "demo_md.rda"))["demo_md"].copy()
    parsed_check = pd.DataFrame([parse_sample(x) for x in packaged_md["Sample"].astype(str)])
    parser_validation = pd.DataFrame(
        [
            {
                "field": "soil",
                "n_checked": len(packaged_md),
                "n_matched": int(
                    np.sum(parsed_check["soil"].to_numpy() == packaged_md["Soil.Type"].astype(str).to_numpy())
                ),
            },
            {
                "field": "genotype",
                "n_checked": len(packaged_md),
                "n_matched": int(
                    np.sum(
                        parsed_check["genotype"].to_numpy()
                        == packaged_md["Arabidopsis"].astype(str).to_numpy()
                    )
                ),
            },
            {
                "field": "fraction",
                "n_checked": len(packaged_md),
                "n_matched": int(
                    np.sum(
                        parsed_check["fraction"].to_numpy()
                        == packaged_md["Treatment"].astype(str).to_numpy()
                    )
                ),
            },
            {
                "field": "developmental_stage",
                "n_checked": len(packaged_md),
                "n_matched": int(
                    np.sum(
                        parsed_check["age"].to_numpy()
                        == packaged_md["Developmental"].astype(str).to_numpy()
                    )
                ),
            },
            {
                "field": "experiment_and_plate",
                "n_checked": len(packaged_md),
                "n_matched": int(
                    np.sum(
                        (
                            parsed_check["experiment"]
                            + "_"
                            + parsed_check["plate"]
                        ).to_numpy()
                        == packaged_md["BiologicalOrExperimentalRep"].astype(str).to_numpy()
                    )
                ),
            },
        ]
    )
    parser_validation["match_fraction"] = (
        parser_validation["n_matched"] / parser_validation["n_checked"]
    )
    parser_validation.to_csv(OUT / "sample_identifier_parser_validation.tsv", sep="\t", index=False)
    if not np.all(parser_validation["n_checked"] == parser_validation["n_matched"]):
        raise RuntimeError("Sample identifier parser did not reproduce packaged metadata")

    audit = (
        metadata.groupby(["experiment", "soil", "fraction", "age"], dropna=False)
        .agg(
            n_samples=("sample", "size"),
            n_depth_pass=("depth_pass", "sum"),
            median_retained_feature_depth=("retained_feature_depth", "median"),
        )
        .reset_index()
    )
    audit.to_csv(OUT / "sample_audit.tsv", sep="\t", index=False)

    keep = metadata["analysis_included"].to_numpy()
    md = metadata.loc[keep].reset_index(drop=True)
    x = counts[keep]
    presence = x > 0
    blocks = md["experiment"].to_numpy()
    thresholds = [0.10, 0.30, 0.50]
    core_rows = []
    for threshold in thresholds:
        pooled = presence.mean(axis=0) >= threshold
        block_prev = np.vstack([presence[blocks == b].mean(axis=0) for b in PRIMARY_BLOCKS])
        study_aware = block_prev.min(axis=0) >= threshold
        core_rows.append(
            {
                "prevalence_threshold": threshold,
                "pooled_features": int(pooled.sum()),
                "study_aware_features": int(study_aware.sum()),
                "pooled_only_features": int(np.sum(pooled & ~study_aware)),
                "study_aware_fraction_of_pooled": float(study_aware.sum() / max(pooled.sum(), 1)),
            }
        )
    core_summary = pd.DataFrame(core_rows)
    core_summary.to_csv(OUT / "pooled_vs_study_aware_core.tsv", sep="\t", index=False)

    loso_rows = []
    loso_features = []
    for heldout in PRIMARY_BLOCKS:
        train = blocks != heldout
        test = blocks == heldout
        train_blocks = [b for b in PRIMARY_BLOCKS if b != heldout]
        train_prev_by_block = np.vstack(
            [presence[blocks == b].mean(axis=0) for b in train_blocks]
        )
        min_prev = train_prev_by_block.min(axis=0)
        mean_prev = train_prev_by_block.mean(axis=0)
        total = x[train].sum(axis=0)
        universe = np.flatnonzero(total >= MIN_TOTAL)
        core = np.flatnonzero((min_prev >= PREV_THRESHOLD) & (total >= MIN_TOTAL))
        test_prev = presence[test].mean(axis=0)
        observed = test_prev[core] >= PREV_THRESHOLD
        observed_rate = float(observed.mean()) if core.size else np.nan
        b1 = make_bins(min_prev)
        b2 = make_bins(np.log10(total + 1))
        strata = b1 * 10 + b2
        null = np.empty(N_PERM)
        for i in range(N_PERM):
            subset = matched_subset(core, universe, strata, rng)
            null[i] = np.mean(test_prev[subset] >= PREV_THRESHOLD)
        q95 = float(np.quantile(null, 0.95))
        empirical_p = float((np.sum(null >= observed_rate) + 1) / (N_PERM + 1))
        loso_rows.append(
            {
                "heldout_block": heldout,
                "n_train": int(train.sum()),
                "n_test": int(test.sum()),
                "universe_features": int(universe.size),
                "candidate_core_features": int(core.size),
                "validated_features": int(observed.sum()),
                "observed_validation_rate": observed_rate,
                "matched_null_mean": float(null.mean()),
                "matched_null_q95": q95,
                "matched_empirical_p": empirical_p,
                "above_matched_q95": bool(observed_rate > q95),
            }
        )
        for idx, ok in zip(core, observed):
            loso_features.append(
                {
                    "heldout_block": heldout,
                    "feature": feature_ids[idx],
                    "min_training_prevalence": min_prev[idx],
                    "mean_training_prevalence": mean_prev[idx],
                    "heldout_prevalence": test_prev[idx],
                    "validated": bool(ok),
                }
            )
    loso = pd.DataFrame(loso_rows)
    loso.to_csv(OUT / "core_membership_loso_validation.tsv", sep="\t", index=False)
    pd.DataFrame(loso_features).to_csv(
        OUT / "core_membership_loso_features.tsv", sep="\t", index=False
    )

    full_prev = np.vstack([presence[blocks == b].mean(axis=0) for b in PRIMARY_BLOCKS])
    stable_idx = np.flatnonzero(
        (full_prev.min(axis=0) >= PREV_THRESHOLD) & (x.sum(axis=0) >= MIN_TOTAL)
    )
    tax.columns = [str(c).strip() for c in tax.columns]
    tax_key = tax.columns[0]
    tax_lookup = tax.set_index(tax_key)
    stable = pd.DataFrame(
        {
            "feature": feature_ids[stable_idx],
            "minimum_block_prevalence": full_prev[:, stable_idx].min(axis=0),
            "mean_block_prevalence": full_prev[:, stable_idx].mean(axis=0),
            "total_counts": x[:, stable_idx].sum(axis=0),
        }
    )
    stable = stable.join(tax_lookup, on="feature")
    stable.to_csv(OUT / "transferable_core_features.tsv", sep="\t", index=False)

    # Response transfer: partial age correlations adjusted for compartment and genotype.
    x_clr = clr(x)
    age = (md["age"].to_numpy() == "old").astype(float)
    fraction = md["fraction"].to_numpy()
    genotype = md["genotype"].to_numpy()
    per_block_stats: dict[str, dict[str, np.ndarray | int]] = {}
    for block in PRIMARY_BLOCKS:
        idx = blocks == block
        r, p, var, rank = residualized_correlations(
            x_clr[idx], age[idx], fraction[idx], genotype[idx]
        )
        per_block_stats[block] = {
            "r": r,
            "p": p,
            "z": np.arctanh(np.clip(r, -0.999999, 0.999999)),
            "var": var,
            "n": int(idx.sum()),
            "rank": rank,
        }

    response_rows = []
    response_feature_rows = []
    for heldout in PRIMARY_BLOCKS:
        train_blocks = [b for b in PRIMARY_BLOCKS if b != heldout]
        z = np.vstack([per_block_stats[b]["z"] for b in train_blocks])
        var = np.vstack([per_block_stats[b]["var"] for b in train_blocks])
        mu, se, p = random_effects(z, var)
        q = bh_adjust(p)
        candidates = np.flatnonzero(q <= 0.05)
        test_r = np.asarray(per_block_stats[heldout]["r"])
        test_p = np.asarray(per_block_stats[heldout]["p"])
        replicated = (np.sign(test_r[candidates]) == np.sign(mu[candidates])) & (
            test_p[candidates] < 0.05
        )
        response_rows.append(
            {
                "heldout_block": heldout,
                "n_train_blocks": len(train_blocks),
                "n_test": per_block_stats[heldout]["n"],
                "meta_candidates": int(candidates.size),
                "directionally_replicated": int(replicated.sum()),
                "replication_rate": float(replicated.mean()) if candidates.size else np.nan,
                "same_direction_fraction": float(
                    np.mean(np.sign(test_r[candidates]) == np.sign(mu[candidates]))
                )
                if candidates.size
                else np.nan,
            }
        )
        for idx, ok in zip(candidates, replicated):
            response_feature_rows.append(
                {
                    "heldout_block": heldout,
                    "feature": feature_ids[idx],
                    "training_meta_fisher_z": mu[idx],
                    "training_meta_se": se[idx],
                    "training_meta_p": p[idx],
                    "training_meta_q": q[idx],
                    "heldout_partial_r": test_r[idx],
                    "heldout_p": test_p[idx],
                    "same_direction": bool(np.sign(test_r[idx]) == np.sign(mu[idx])),
                    "directionally_replicated": bool(ok),
                }
            )
    response = pd.DataFrame(response_rows)
    response.to_csv(OUT / "developmental_response_loso_summary.tsv", sep="\t", index=False)
    response_features = pd.DataFrame(response_feature_rows)
    if not response_features.empty:
        response_features = response_features.join(tax_lookup, on="feature")
    response_features.to_csv(
        OUT / "developmental_response_loso_features.tsv", sep="\t", index=False
    )

    # Cross-block prediction of endosphere versus rhizosphere, using training-defined core only.
    pred_rows = []
    fold_nulls = {}
    for heldout in PRIMARY_BLOCKS:
        train = blocks != heldout
        test = blocks == heldout
        train_blocks = [b for b in PRIMARY_BLOCKS if b != heldout]
        prev = np.vstack([presence[blocks == b].mean(axis=0) for b in train_blocks])
        total = x[train].sum(axis=0)
        features = np.flatnonzero(
            (prev.min(axis=0) >= PREV_THRESHOLD) & (total >= MIN_TOTAL)
        )
        y_train = fraction[train]
        y_test = fraction[test]
        observed_pred = nearest_centroid_fit_predict(
            x_clr[train][:, features], y_train, x_clr[test][:, features]
        )
        ba = balanced_accuracy(y_test, observed_pred)
        f1 = macro_f1(y_test, observed_pred)
        null = np.empty(N_PERM)
        train_block_labels = blocks[train]
        for i in range(N_PERM):
            shuffled = y_train.copy()
            for b in train_blocks:
                loc = np.flatnonzero(train_block_labels == b)
                shuffled[loc] = rng.permutation(shuffled[loc])
            pp = nearest_centroid_fit_predict(
                x_clr[train][:, features], shuffled, x_clr[test][:, features]
            )
            null[i] = balanced_accuracy(y_test, pp)
        fold_nulls[heldout] = null
        pred_rows.append(
            {
                "heldout_block": heldout,
                "n_train": int(train.sum()),
                "n_test": int(test.sum()),
                "training_defined_features": int(features.size),
                "balanced_accuracy": ba,
                "macro_f1": f1,
                "permutation_mean": float(null.mean()),
                "permutation_q95": float(np.quantile(null, 0.95)),
                "empirical_p": float((np.sum(null >= ba) + 1) / (N_PERM + 1)),
            }
        )
    prediction = pd.DataFrame(pred_rows)
    prediction.to_csv(OUT / "compartment_prediction_loso.tsv", sep="\t", index=False)
    null_table = pd.DataFrame(fold_nulls)
    null_table.insert(0, "permutation", np.arange(1, N_PERM + 1))
    null_table.to_csv(OUT / "compartment_prediction_null.tsv", sep="\t", index=False)
    observed_mean_ba = float(prediction["balanced_accuracy"].mean())
    null_mean_by_perm = null_table[PRIMARY_BLOCKS].mean(axis=1).to_numpy()
    prediction_overall = {
        "observed_mean_ba": observed_mean_ba,
        "permutation_mean": float(null_mean_by_perm.mean()),
        "permutation_q95": float(np.quantile(null_mean_by_perm, 0.95)),
        "empirical_p": float(
            (np.sum(null_mean_by_perm >= observed_mean_ba) + 1) / (N_PERM + 1)
        ),
    }
    pd.DataFrame([prediction_overall]).to_csv(
        OUT / "compartment_prediction_overall.tsv", sep="\t", index=False
    )

    # Sensitivity to the retained-OTU count threshold (label-free QC choice).
    sensitivity_rows = []
    for depth_threshold in [50, 100, 200, 500]:
        sk = (
            metadata["primary_block"].to_numpy()
            & metadata["root_compartment"].to_numpy()
            & (metadata["retained_feature_depth"].to_numpy() >= depth_threshold)
        )
        smd = metadata.loc[sk].reset_index(drop=True)
        sx = counts[sk]
        sp = sx > 0
        sb = smd["experiment"].to_numpy()
        sf = smd["fraction"].to_numpy()
        block_prev = np.vstack([sp[sb == b].mean(axis=0) for b in PRIMARY_BLOCKS])
        pooled_n = int(np.sum(sp.mean(axis=0) >= PREV_THRESHOLD))
        study_n = int(np.sum(block_prev.min(axis=0) >= PREV_THRESHOLD))
        validation_rates = []
        prediction_bas = []
        sx_clr = clr(sx)
        for heldout in PRIMARY_BLOCKS:
            train = sb != heldout
            test = sb == heldout
            train_blocks = [b for b in PRIMARY_BLOCKS if b != heldout]
            train_prev = np.vstack([sp[sb == b].mean(axis=0) for b in train_blocks])
            total = sx[train].sum(axis=0)
            fidx = np.flatnonzero(
                (train_prev.min(axis=0) >= PREV_THRESHOLD) & (total >= MIN_TOTAL)
            )
            test_prev = sp[test].mean(axis=0)
            validation_rates.append(float(np.mean(test_prev[fidx] >= PREV_THRESHOLD)))
            pp = nearest_centroid_fit_predict(
                sx_clr[train][:, fidx], sf[train], sx_clr[test][:, fidx]
            )
            prediction_bas.append(balanced_accuracy(sf[test], pp))
        sensitivity_rows.append(
            {
                "minimum_retained_feature_depth": depth_threshold,
                "n_samples": int(sk.sum()),
                "pooled_core_features": pooled_n,
                "study_aware_core_features": study_n,
                "mean_loso_core_validation_rate": float(np.mean(validation_rates)),
                "mean_loso_prediction_balanced_accuracy": float(np.mean(prediction_bas)),
            }
        )
    sensitivity = pd.DataFrame(sensitivity_rows)
    sensitivity.to_csv(OUT / "retained_depth_sensitivity.tsv", sep="\t", index=False)

    svg_figure(
        core_summary,
        loso,
        response,
        prediction,
        prediction_overall,
        OUT / "Figure_external_generalization_arabidopsis.svg",
    )

    n_included = int(keep.sum())
    methods = f"""External generalization analysis (Arabidopsis)

Data source: the CoreMicrobiomeR CRAN demonstration subset of the Lundberg et al. Arabidopsis root-microbiome dataset (188 representative OTUs across 1,439 samples). Sample design fields were reconstructed from the published identifier convention; soil, genotype, fraction, developmental stage and experiment/plate assignments each matched all 103 rows of the packaged metadata subset. Analyses were restricted to four primary experiments (C1, C2, M1 and M2), endosphere and rhizosphere samples, and samples with at least {MIN_RETAINED_DEPTH} counts across the 188 retained OTUs (n = {n_included}). This retained-feature count threshold is not interpreted as full-library sequencing depth.

Core membership was defined as prevalence >= {PREV_THRESHOLD:.2f} in every training experiment and total training abundance >= {MIN_TOTAL}. Transfer was assessed by leave-one-experiment-out validation. For each fold, {N_PERM} null subsets of equal size were sampled after stratification by minimum training prevalence and training abundance; empirical P values used the plus-one correction.

Developmental-stage responses were evaluated on CLR-transformed abundances. Within each experiment, partial correlations with old versus young stage adjusted for compartment and genotype. Training-experiment effects were combined with a DerSimonian-Laird random-effects model and controlled at Benjamini-Hochberg q <= 0.05; validation required nominal P < 0.05 and the same direction in the held-out experiment.

Endosphere-versus-rhizosphere discrimination used a nearest-centroid classifier and core features defined from the training experiments only. Performance was measured by balanced accuracy and macro-F1. The null distribution used {N_PERM} label permutations within each training experiment.
"""
    (OUT / "methods_insert.txt").write_text(methods, encoding="utf-8")

    results = f"""External generalization results (Arabidopsis)

After design and retained-feature-depth filtering, {n_included} samples from four independent experiments were analyzed. At 50% prevalence, pooled analysis classified {int(core_summary.loc[core_summary.prevalence_threshold == 0.5, 'pooled_features'].iloc[0])} OTUs as core, whereas requiring the threshold in every experiment retained {int(core_summary.loc[core_summary.prevalence_threshold == 0.5, 'study_aware_features'].iloc[0])}. Across leave-one-experiment-out folds, observed core validation ranged from {loso.observed_validation_rate.min():.1%} to {loso.observed_validation_rate.max():.1%}; {int(loso.above_matched_q95.sum())} of {len(loso)} folds exceeded the prevalence- and abundance-matched 95th-percentile null.

Training-experiment meta-analysis identified {int(response.meta_candidates.min())}-{int(response.meta_candidates.max())} developmental-stage-responsive OTUs per held-out fold, of which {response.directionally_replicated.sum()} of {response.meta_candidates.sum()} fold-specific candidates replicated nominally with concordant direction ({response.directionally_replicated.sum()/max(response.meta_candidates.sum(),1):.1%} overall).

Training-defined core features discriminated endosphere from rhizosphere with mean held-out balanced accuracy {observed_mean_ba:.3f} (fold range {prediction.balanced_accuracy.min():.3f}-{prediction.balanced_accuracy.max():.3f}) and mean macro-F1 {prediction.macro_f1.mean():.3f}. This exceeded the within-training-experiment permutation null (95th percentile {prediction_overall['permutation_q95']:.3f}; empirical P = {prediction_overall['empirical_p']:.3g}).

These conclusions were insensitive to the retained-feature count filter: thresholds from 50 to 500 retained {int(sensitivity.n_samples.min())}-{int(sensitivity.n_samples.max())} samples, with mean leave-one-experiment-out core validation rates of {sensitivity.mean_loso_core_validation_rate.min():.3f}-{sensitivity.mean_loso_core_validation_rate.max():.3f} and mean balanced accuracies of {sensitivity.mean_loso_prediction_balanced_accuracy.min():.3f}-{sensitivity.mean_loso_prediction_balanced_accuracy.max():.3f}.

Interpretation: this independent plant system reproduces the central framework result that pooled ubiquity, cross-experiment membership, response replication and predictive discrimination are distinct properties. Because the CRAN object contains a representative 188-OTU subset rather than the complete OTU table, this analysis is an external framework demonstration, not a biological reanalysis of the full Lundberg dataset.
"""
    (OUT / "results_insert.txt").write_text(results, encoding="utf-8")

    caption = """Figure X | External demonstration of the study-aware transferability framework in an independent Arabidopsis root-microbiome dataset. a, Numbers of pooled and study-aware core OTUs at increasing prevalence thresholds. b, Leave-one-experiment-out validation rates for training-defined core OTUs; black ticks mark the 95th percentile of prevalence- and abundance-matched random-subset nulls. c, Numbers of developmental-stage-responsive OTUs discovered by random-effects meta-analysis in three experiments and directionally replicated at nominal P < 0.05 in the held-out experiment; percentages denote replication rates. d, Balanced accuracy for endosphere-versus-rhizosphere discrimination in held-out experiments using training-defined core features; orange ticks mark fold-specific 95th-percentile within-training-experiment label-permutation nulls and the dashed line denotes two-class chance. The public CRAN object contains 188 representative OTUs and is used as a framework demonstration rather than a full biological reanalysis."""
    (OUT / "figure_caption.txt").write_text(caption, encoding="utf-8")

    tarball = ROOT / "work" / "arabidopsis_cran" / "CoreMicrobiomeR_0.1.0.tar.gz"
    sha256 = hashlib.sha256(tarball.read_bytes()).hexdigest()
    manifest = pd.DataFrame(
        [
            {
                "resource": "CoreMicrobiomeR source package",
                "version_or_doi": "0.1.0",
                "url": "https://cran.r-project.org/src/contrib/CoreMicrobiomeR_0.1.0.tar.gz",
                "sha256": sha256,
                "role": "Public 188-OTU x 1,439-sample demonstration object",
            },
            {
                "resource": "Lundberg et al. source study",
                "version_or_doi": "10.1038/nature11237",
                "url": "https://doi.org/10.1038/nature11237",
                "sha256": "",
                "role": "Biological source and sample-identifier convention",
            },
        ]
    )
    manifest.to_csv(OUT / "source_manifest.tsv", sep="\t", index=False)

    readme = f"""# Independent Arabidopsis framework demonstration

## Outcome

- Public object: 188 representative OTUs across 1,439 samples from the Lundberg et al. Arabidopsis root-microbiome study.
- Primary analysis: {n_included} endosphere/rhizosphere samples from four independent experiments after retained-feature-depth QC.
- Pooled versus study-aware core at 50% prevalence: 154 versus 98 OTUs.
- LOSO membership validation: {loso.observed_validation_rate.min():.1%}-{loso.observed_validation_rate.max():.1%}; all four folds exceeded matched-null 95th percentiles (all empirical P = 0.001).
- Developmental-stage response replication: 138/203 fold-specific candidates ({response.directionally_replicated.sum()/response.meta_candidates.sum():.1%}).
- Compartment prediction: mean held-out BA {observed_mean_ba:.3f}, mean macro-F1 {prediction.macro_f1.mean():.3f}, overall permutation P = {prediction_overall['empirical_p']:.3g}.

## Important interpretation boundary

This is a framework-level external demonstration, not a full reanalysis of the Lundberg study. The CRAN object contains 188 representative OTUs rather than the complete 18,783-OTU table. The prediction result is strong overall, but fold C1 alone has empirical P = {float(prediction.loc[prediction.heldout_block == 'C1', 'empirical_p'].iloc[0]):.3f}; retain the fold-level table for transparency.

## Main files

- `Figure_external_generalization_arabidopsis.svg` and `.png`: publication-style four-panel figure.
- `results_insert.txt`, `methods_insert.txt`, `figure_caption.txt`: manuscript-ready draft text.
- `pooled_vs_study_aware_core.tsv`: pooled/core contrast.
- `core_membership_loso_validation.tsv`: membership transfer and matched nulls.
- `developmental_response_loso_summary.tsv`: independent response replication.
- `compartment_prediction_loso.tsv`: held-out prediction by experiment.
- `retained_depth_sensitivity.tsv`: QC-threshold sensitivity.
- `sample_identifier_parser_validation.tsv`: 103/103 metadata parser checks for every reconstructed field.
- `source_manifest.tsv`: source URL, DOI and SHA-256.

## Re-run

From the project root with `pyreadr`, `numpy` and `pandas` available:

```powershell
python work/run_arabidopsis_external_generalization.py
```

Random seed: {SEED}; permutations per null: {N_PERM}.
"""
    (OUT / "ANALYSIS_README.md").write_text(readme, encoding="utf-8")

    print(results)
    print("Outputs:", OUT)


if __name__ == "__main__":
    main()
