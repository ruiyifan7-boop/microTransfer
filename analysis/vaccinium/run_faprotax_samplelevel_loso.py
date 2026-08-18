#!/usr/bin/env python3
"""Run official FAPROTAX and seven-cohort sample-level LOSO validation.

Required inputs
---------------
1. Taxon table: rows are bacterial taxonomic paths, columns are samples, and
   the first column is named ``taxonomy``. Plain TSV or TSV.GZ is accepted.
2. Metadata: tab-separated file containing ``sample_uid`` and ``study``.

The script calls the official FAPROTAX ``collapse_table.py`` without any
manual database extension, derives sample-level function abundances, and then
performs leave-one-study-out prevalence validation. Candidate functions are
defined independently in each training fold as present in at least 50% of
samples in every training cohort. The held-out cohort is never used for
candidate definition.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import random
import shutil
import subprocess
import sys
import tempfile
from typing import Iterable

import numpy as np
import pandas as pd


DEFAULT_SEED = 20260728


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--taxon-table", type=Path)
    parser.add_argument("--metadata", type=Path)
    parser.add_argument("--faprotax-dir", type=Path)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--prevalence-threshold", type=float, default=0.50)
    parser.add_argument("--n-permutations", type=int, default=999)
    parser.add_argument("--seed", type=int, default=DEFAULT_SEED)
    parser.add_argument("--expected-studies", type=int, default=7)
    parser.add_argument("--expected-samples", type=int, default=265)
    parser.add_argument(
        "--self-test",
        action="store_true",
        help="Run a synthetic seven-cohort smoke test instead of real data.",
    )
    args = parser.parse_args()
    if not args.self_test:
        required = ["taxon_table", "metadata", "faprotax_dir", "output_dir"]
        missing = [x for x in required if getattr(args, x) is None]
        if missing:
            parser.error("Missing required arguments: " + ", ".join(missing))
    return args


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def open_text(path: Path):
    return gzip.open(path, "rt", encoding="utf-8-sig", newline="") if path.suffix == ".gz" else path.open(
        "r", encoding="utf-8-sig", newline=""
    )


def materialize_plain_tsv(source: Path, destination: Path) -> Path:
    if source.suffix != ".gz":
        return source
    with gzip.open(source, "rb") as src, destination.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return destination


def write_tsv(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False, na_rep="")


def validate_inputs(
    taxon_table: Path,
    metadata_path: Path,
    expected_studies: int,
    expected_samples: int,
) -> tuple[pd.DataFrame, list[str], pd.Series]:
    metadata = pd.read_csv(metadata_path, sep="\t", dtype=str)
    required = {"sample_uid", "study"}
    if not required.issubset(metadata.columns):
        raise ValueError(f"Metadata must contain {sorted(required)}")
    metadata = metadata.loc[:, ~metadata.columns.duplicated()].copy()
    metadata["sample_uid"] = metadata["sample_uid"].astype(str)
    metadata["study"] = metadata["study"].astype(str)
    if metadata["sample_uid"].duplicated().any():
        raise ValueError("Metadata contains duplicated sample_uid values")
    if metadata[["sample_uid", "study"]].isna().any().any():
        raise ValueError("Metadata contains missing sample_uid or study values")

    with open_text(taxon_table) as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader)
        if not header or header[0].lstrip("#") != "taxonomy":
            raise ValueError("Taxon table first column must be named taxonomy")
        sample_ids = header[1:]
        if len(sample_ids) != len(set(sample_ids)):
            raise ValueError("Taxon table contains duplicated sample columns")
        if not sample_ids:
            raise ValueError("Taxon table contains no samples")
        library_totals = np.zeros(len(sample_ids), dtype=float)
        n_features = 0
        for row_number, row in enumerate(reader, start=2):
            if not row or (len(row) == 1 and not row[0].strip()):
                continue
            if len(row) != len(header):
                raise ValueError(
                    f"Taxon table row {row_number} has {len(row)} fields; expected {len(header)}"
                )
            try:
                values = np.asarray(row[1:], dtype=float)
            except ValueError as exc:
                raise ValueError(f"Non-numeric count at taxon-table row {row_number}") from exc
            if np.any(~np.isfinite(values)) or np.any(values < 0):
                raise ValueError(f"Invalid count at taxon-table row {row_number}")
            library_totals += values
            n_features += 1
    if n_features == 0:
        raise ValueError("Taxon table contains no features")
    if np.any(library_totals <= 0):
        bad = [sample_ids[i] for i in np.flatnonzero(library_totals <= 0)]
        raise ValueError(f"Samples with zero total bacterial counts: {bad[:10]}")

    missing_metadata = sorted(set(sample_ids) - set(metadata["sample_uid"]))
    missing_counts = sorted(set(metadata["sample_uid"]) - set(sample_ids))
    if missing_metadata or missing_counts:
        raise ValueError(
            "Sample mismatch. "
            f"Counts without metadata={missing_metadata[:10]}; "
            f"metadata without counts={missing_counts[:10]}"
        )
    metadata = metadata.set_index("sample_uid").loc[sample_ids].reset_index()
    studies = sorted(metadata["study"].unique())
    if len(studies) != expected_studies:
        raise ValueError(f"Expected {expected_studies} studies; found {len(studies)}: {studies}")
    if expected_samples > 0 and len(sample_ids) != expected_samples:
        raise ValueError(f"Expected {expected_samples} samples; found {len(sample_ids)}")
    return metadata, sample_ids, pd.Series(library_totals, index=sample_ids, name="total_16S_reads")


def run_faprotax(
    taxon_table: Path,
    faprotax_dir: Path,
    output_dir: Path,
) -> dict[str, Path]:
    collapse = faprotax_dir / "collapse_table.py"
    database = faprotax_dir / "FAPROTAX.txt"
    if not collapse.is_file() or not database.is_file():
        raise FileNotFoundError("FAPROTAX directory must contain collapse_table.py and FAPROTAX.txt")

    outputs = {
        "collapsed": output_dir / "sample_function_abundance_raw.tsv",
        "report": output_dir / "faprotax_assignment_report.txt",
        "log": output_dir / "faprotax_run.log",
        "mapping": output_dir / "taxon_function_mapping.tsv",
        "definitions_used": output_dir / "faprotax_definitions_used.txt",
        "definitions_unused": output_dir / "faprotax_definitions_unused.txt",
    }
    command = [
        sys.executable,
        "-X",
        "utf8",
        str(collapse),
        "-i",
        str(taxon_table),
        "-g",
        str(database),
        "-o",
        str(outputs["collapsed"]),
        "-r",
        str(outputs["report"]),
        "-l",
        str(outputs["log"]),
        "--out_groups2records_table",
        str(outputs["mapping"]),
        "--out_group_definitions_used",
        str(outputs["definitions_used"]),
        "--out_group_definitions_unused",
        str(outputs["definitions_unused"]),
        "--group_leftovers_as",
        "unassigned_taxa",
        "--omit_unrepresented_groups",
        "-d",
        "taxonomy",
        "--column_names_are_in",
        "first_data_line",
        "-c",
        "#",
        "-v",
        "-f",
    ]
    environment = os.environ.copy()
    environment["PYTHONUTF8"] = "1"
    completed = subprocess.run(
        command,
        text=True,
        capture_output=True,
        env=environment,
    )
    (output_dir / "faprotax_console.txt").write_text(
        completed.stdout + "\n--- STDERR ---\n" + completed.stderr,
        encoding="utf-8",
    )
    if completed.returncode != 0:
        raise RuntimeError(
            f"FAPROTAX failed with exit code {completed.returncode}; "
            f"see {output_dir / 'faprotax_console.txt'}"
        )
    if not outputs["collapsed"].is_file():
        raise RuntimeError("FAPROTAX completed without producing a collapsed function table")
    return outputs


def read_collapsed(path: Path, expected_samples: list[str]) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", comment=None)
    frame.columns = [str(x).lstrip("#") for x in frame.columns]
    sample_set = set(expected_samples)
    sample_columns = [c for c in frame.columns if c in sample_set]
    label_columns = [c for c in frame.columns if c not in sample_set]
    if len(sample_columns) != len(expected_samples):
        missing = sorted(sample_set - set(sample_columns))
        raise ValueError(f"Collapsed table is missing samples: {missing[:10]}")
    if not label_columns:
        raise ValueError("Collapsed table has no functional-group label column")
    group_col = label_columns[0]
    frame = frame.rename(columns={group_col: "function"})
    frame["function"] = frame["function"].astype(str)
    numeric = frame[sample_columns].apply(pd.to_numeric, errors="raise")
    if (numeric < 0).any().any():
        raise ValueError("Collapsed function table contains negative values")
    numeric.index = frame["function"]
    if numeric.index.duplicated().any():
        numeric = numeric.groupby(level=0).sum()
    return numeric.loc[:, expected_samples].T


def prevalence_matched_null(
    function_presence: pd.DataFrame,
    metadata: pd.DataFrame,
    heldout: str,
    candidates: list[str],
    observed_rate: float,
    n_permutations: int,
    rng: np.random.Generator,
    threshold: float,
) -> tuple[pd.DataFrame, float]:
    if not candidates or n_permutations <= 0:
        return pd.DataFrame(), math.nan
    train_mask = metadata["study"].to_numpy() != heldout
    test_mask = ~train_mask
    train_mean = function_presence.loc[train_mask].mean(axis=0)
    heldout_prev = function_presence.loc[test_mask].mean(axis=0)
    functions = train_mean.index.to_numpy()
    ranks = train_mean.rank(method="average", pct=True)
    bins = np.minimum((ranks * 10).astype(int), 9)
    candidate_bins = bins.loc[candidates].to_numpy()
    candidate_set = set(candidates)

    rows: list[dict[str, object]] = []
    for permutation in range(1, n_permutations + 1):
        selected: list[str] = []
        for bin_id in candidate_bins:
            same_bin = functions[bins.loc[functions].to_numpy() == bin_id]
            preferred = [x for x in same_bin if x not in candidate_set and x not in selected]
            pool = preferred or [x for x in same_bin if x not in selected] or list(functions)
            selected.append(str(rng.choice(pool)))
        rate = float((heldout_prev.loc[selected] >= threshold).mean())
        rows.append(
            {
                "heldout_study": heldout,
                "permutation": permutation,
                "matched_function_count": len(selected),
                "null_validation_rate": rate,
            }
        )
    null = pd.DataFrame(rows)
    p_value = (1 + int((null["null_validation_rate"] >= observed_rate).sum())) / (
        n_permutations + 1
    )
    return null, p_value


def run_loso(
    function_abundance: pd.DataFrame,
    metadata: pd.DataFrame,
    threshold: float,
    n_permutations: int,
    seed: int,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    presence = function_abundance.gt(0)
    studies = sorted(metadata["study"].unique())
    fold_rows: list[dict[str, object]] = []
    detail_frames: list[pd.DataFrame] = []
    null_frames: list[pd.DataFrame] = []
    rng = np.random.default_rng(seed)

    for heldout in studies:
        train_studies = [s for s in studies if s != heldout]
        training_prevalence = pd.DataFrame(
            {
                study: presence.loc[metadata["study"].eq(study).to_numpy()].mean(axis=0)
                for study in train_studies
            }
        ).T
        heldout_prevalence = presence.loc[
            metadata["study"].eq(heldout).to_numpy()
        ].mean(axis=0)
        minimum_training = training_prevalence.min(axis=0)
        mean_training = training_prevalence.mean(axis=0)
        candidates = minimum_training.index[minimum_training >= threshold].tolist()
        validated = [f for f in candidates if heldout_prevalence.loc[f] >= threshold]
        observed_rate = len(validated) / len(candidates) if candidates else math.nan

        details = pd.DataFrame(
            {
                "heldout_study": heldout,
                "function": candidates,
                "minimum_training_prevalence": minimum_training.loc[candidates].to_numpy(),
                "mean_training_prevalence": mean_training.loc[candidates].to_numpy(),
                "heldout_prevalence": heldout_prevalence.loc[candidates].to_numpy(),
            }
        )
        if not details.empty:
            details["validated"] = details["heldout_prevalence"] >= threshold
            details["training_studies"] = ",".join(train_studies)
            detail_frames.append(details)

        null, p_value = prevalence_matched_null(
            presence,
            metadata,
            heldout,
            candidates,
            observed_rate,
            n_permutations,
            rng,
            threshold,
        )
        if not null.empty:
            null_frames.append(null)
        fold_rows.append(
            {
                "heldout_study": heldout,
                "training_studies": ",".join(train_studies),
                "n_train": int(metadata["study"].ne(heldout).sum()),
                "n_test": int(metadata["study"].eq(heldout).sum()),
                "functions_tested": int(function_abundance.shape[1]),
                "discovery_core_functions": len(candidates),
                "validated_core_functions": len(validated),
                "validation_rate": observed_rate,
                "matched_null_mean": float(null["null_validation_rate"].mean())
                if not null.empty
                else math.nan,
                "matched_null_p": p_value,
            }
        )

    fold_summary = pd.DataFrame(fold_rows)
    details_all = pd.concat(detail_frames, ignore_index=True) if detail_frames else pd.DataFrame()
    null_all = pd.concat(null_frames, ignore_index=True) if null_frames else pd.DataFrame()
    if details_all.empty:
        stable = pd.DataFrame(
            columns=[
                "function",
                "core_in_heldout_definitions",
                "validated_in_heldout_definitions",
                "minimum_heldout_prevalence",
                "mean_heldout_prevalence",
                "stable_all_heldouts",
            ]
        )
    else:
        stable = (
            details_all.groupby("function", as_index=False)
            .agg(
                core_in_heldout_definitions=("heldout_study", "nunique"),
                validated_in_heldout_definitions=("validated", "sum"),
                minimum_heldout_prevalence=("heldout_prevalence", "min"),
                mean_heldout_prevalence=("heldout_prevalence", "mean"),
            )
            .sort_values(
                ["validated_in_heldout_definitions", "minimum_heldout_prevalence"],
                ascending=False,
            )
        )
        stable["stable_all_heldouts"] = (
            stable["core_in_heldout_definitions"].eq(len(studies))
            & stable["validated_in_heldout_definitions"].eq(len(studies))
        )
    return fold_summary, details_all, stable, null_all


def build_self_test(base: Path) -> tuple[Path, Path]:
    studies = [f"COHORT_{i}" for i in range(1, 8)]
    sample_ids = [f"{study}__S{j}" for study in studies for j in range(1, 5)]
    metadata = pd.DataFrame(
        {"sample_uid": sample_ids, "study": [s for s in studies for _ in range(4)]}
    )
    taxa = [
        "k__Bacteria;p__Proteobacteria;c__Alphaproteobacteria;o__Rhizobiales;f__Xanthobacteraceae;g__Bradyrhizobium",
        "k__Bacteria;p__Proteobacteria;c__Gammaproteobacteria;o__Burkholderiales;f__Burkholderiaceae;g__Burkholderia",
        "k__Bacteria;p__Actinobacteriota;c__Actinobacteria;o__Frankiales;f__Acidothermaceae;g__Acidothermus",
    ]
    matrix = np.asarray(
        [
            [100 for j in range(len(sample_ids))],
            [50 if j % 3 != 0 else 0 for j in range(len(sample_ids))],
            [30 if (j // 4) % 2 == 0 else 0 for j in range(len(sample_ids))],
        ],
        dtype=int,
    )
    taxon_path = base / "self_test_taxon_table.tsv"
    metadata_path = base / "self_test_metadata.tsv"
    pd.DataFrame(matrix, columns=sample_ids).assign(taxonomy=taxa).loc[
        :, ["taxonomy", *sample_ids]
    ].to_csv(taxon_path, sep="\t", index=False)
    metadata.to_csv(metadata_path, sep="\t", index=False)
    return taxon_path, metadata_path


def main() -> None:
    args = parse_args()
    if args.self_test:
        if args.faprotax_dir is None:
            raise SystemExit("--self-test still requires --faprotax-dir")
        output_dir = args.output_dir or Path("faprotax_loso_self_test")
        output_dir.mkdir(parents=True, exist_ok=True)
        taxon_table, metadata_path = build_self_test(output_dir)
        expected_samples = 28
    else:
        output_dir = args.output_dir
        output_dir.mkdir(parents=True, exist_ok=True)
        taxon_table = args.taxon_table
        metadata_path = args.metadata
        expected_samples = args.expected_samples

    plain_taxon_table = materialize_plain_tsv(
        taxon_table,
        output_dir / "faprotax_input_taxon_table.tsv",
    )
    metadata, sample_ids, library_totals = validate_inputs(
        plain_taxon_table,
        metadata_path,
        args.expected_studies,
        expected_samples,
    )
    faprotax_outputs = run_faprotax(plain_taxon_table, args.faprotax_dir, output_dir)
    functions_by_sample = read_collapsed(faprotax_outputs["collapsed"], sample_ids)

    unassigned_col = "unassigned_taxa"
    if unassigned_col in functions_by_sample.columns:
        unassigned = functions_by_sample[unassigned_col]
        function_abundance = functions_by_sample.drop(columns=[unassigned_col])
    else:
        # FAPROTAX omits the leftovers group when every input record is assigned.
        unassigned = pd.Series(0.0, index=functions_by_sample.index)
        function_abundance = functions_by_sample
    coverage = pd.DataFrame(
        {
            "sample_uid": sample_ids,
            "study": metadata["study"],
            "total_16S_reads": library_totals.loc[sample_ids].to_numpy(),
            "unassigned_16S_reads": unassigned.loc[sample_ids].to_numpy(),
        }
    )
    coverage["assigned_16S_reads"] = coverage["total_16S_reads"] - coverage["unassigned_16S_reads"]
    coverage["assigned_read_fraction"] = (
        coverage["assigned_16S_reads"] / coverage["total_16S_reads"]
    )
    if (coverage["assigned_read_fraction"] < -1e-9).any():
        raise ValueError("FAPROTAX unassigned counts exceed original library totals")

    fractions = function_abundance.div(library_totals.loc[sample_ids], axis=0)
    raw_out = function_abundance.copy()
    raw_out.insert(0, "sample_uid", raw_out.index)
    raw_out.insert(1, "study", metadata.set_index("sample_uid").loc[raw_out.index, "study"])
    fraction_out = fractions.copy()
    fraction_out.insert(0, "sample_uid", fraction_out.index)
    fraction_out.insert(1, "study", metadata.set_index("sample_uid").loc[fraction_out.index, "study"])

    fold_summary, details, stable, null = run_loso(
        function_abundance,
        metadata,
        args.prevalence_threshold,
        args.n_permutations,
        args.seed,
    )

    write_tsv(raw_out.reset_index(drop=True), output_dir / "sample_function_abundance_counts.tsv")
    write_tsv(
        fraction_out.reset_index(drop=True),
        output_dir / "sample_function_abundance_fraction_of_16S.tsv",
    )
    write_tsv(coverage, output_dir / "sample_annotation_coverage.tsv")
    write_tsv(fold_summary, output_dir / "faprotax_loso_fold_summary.tsv")
    write_tsv(details, output_dir / "faprotax_loso_function_details.tsv")
    write_tsv(stable, output_dir / "faprotax_stable_functions.tsv")
    write_tsv(null, output_dir / "faprotax_loso_prevalence_matched_null.tsv")

    database = args.faprotax_dir / "FAPROTAX.txt"
    manifest = {
        "analysis": "sample-level FAPROTAX with seven-cohort LOSO validation",
        "official_database_only": True,
        "faprotax_database_sha256": sha256(database),
        "taxon_table_sha256": sha256(taxon_table),
        "metadata_sha256": sha256(metadata_path),
        "samples": len(sample_ids),
        "studies": metadata["study"].value_counts().sort_index().to_dict(),
        "functions": int(function_abundance.shape[1]),
        "prevalence_threshold": args.prevalence_threshold,
        "n_permutations": args.n_permutations,
        "seed": args.seed,
        "stable_all_heldouts": int(stable["stable_all_heldouts"].sum())
        if not stable.empty
        else 0,
        "mean_annotation_coverage": float(coverage["assigned_read_fraction"].mean()),
        "minimum_annotation_coverage": float(coverage["assigned_read_fraction"].min()),
    }
    (output_dir / "faprotax_run_manifest.json").write_text(
        json.dumps(manifest, indent=2, ensure_ascii=False),
        encoding="utf-8",
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
