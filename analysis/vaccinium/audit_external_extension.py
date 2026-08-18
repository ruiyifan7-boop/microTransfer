#!/usr/bin/env python3
"""Audit and pair high-priority external Vaccinium 16S/ITS cohorts."""

from __future__ import annotations

import csv
import io
import re
import sys
import time
import urllib.parse
import urllib.request
from collections import Counter, defaultdict
from pathlib import Path


OUT = Path("external_extension")
RUN_FIELDS = [
    "run_accession",
    "sample_accession",
    "secondary_sample_accession",
    "experiment_accession",
    "sample_title",
    "experiment_title",
    "library_name",
    "library_layout",
    "instrument_model",
    "read_count",
    "base_count",
    "fastq_bytes",
    "fastq_ftp",
    "fastq_md5",
]

PROJECTS = [
    "PRJNA835240",
    "PRJNA1253324",
    "PRJNA577971",
    "PRJNA578171",
]

EXPECTED_PAIRS = {
    "PRJNA835240": 24,
    "PRJNA1253324": 8,
    "PRJNA577971_PRJNA578171": 34,
}


def ena_rows(project: str) -> list[dict[str, str]]:
    query = urllib.parse.urlencode(
        {
            "result": "read_run",
            "query": f'study_accession="{project}"',
            "fields": ",".join(RUN_FIELDS),
            "format": "tsv",
            "limit": "0",
        }
    )
    url = "https://www.ebi.ac.uk/ena/portal/api/search?" + query
    cache = OUT / "cache" / f"{project}_ENA_runs.tsv"
    cache.parent.mkdir(parents=True, exist_ok=True)
    text = ""
    last_error: Exception | None = None
    for attempt in range(1, 4):
        try:
            request = urllib.request.Request(
                url, headers={"User-Agent": "iMeta-Vaccinium-audit/1.0"}
            )
            with urllib.request.urlopen(request, timeout=120) as response:
                text = response.read().decode("utf-8")
            if text.strip():
                cache.write_text(text, encoding="utf-8")
                break
        except Exception as exc:
            last_error = exc
            time.sleep(2**attempt)
    if not text.strip() and cache.exists():
        text = cache.read_text(encoding="utf-8")
        print(f"WARNING: using cached ENA response for {project}", file=sys.stderr)
    if not text.strip():
        raise RuntimeError(f"ENA query failed for {project}: {last_error}")
    rows = list(csv.DictReader(io.StringIO(text), delimiter="\t"))
    if not rows:
        raise RuntimeError(f"ENA returned no runs for {project}")
    for row in rows:
        row["source_project"] = project
        row["ena_query_url"] = url
    return rows


def marker_and_biosample(row: dict[str, str]) -> tuple[str, str, str]:
    project = row["source_project"]
    library = row["library_name"].strip()

    if project == "PRJNA835240":
        match = re.fullmatch(r"([bf])(\d+)", library, flags=re.IGNORECASE)
        if not match:
            return "", "", "unrecognized_library_name"
        marker = "16S" if match.group(1).lower() == "b" else "ITS"
        return marker, match.group(2), "library_prefix_b_or_f"

    if project == "PRJNA1253324":
        match = re.fullmatch(r"(.+)-(16S|ITS)", library, flags=re.IGNORECASE)
        if not match:
            return "", "", "unrecognized_library_name"
        return match.group(2).upper(), match.group(1), "library_marker_suffix"

    if project in {"PRJNA577971", "PRJNA578171"}:
        marker = "16S" if project == "PRJNA577971" else "ITS"
        biological = re.sub(r"^E-", "", library)
        biological = re.sub(r"_S\d+_L\d+_R1$", "", biological)
        return marker, biological, "paired_bioproject_and_library_stem"

    return "", "", "unsupported_project"


def study_unit(project: str) -> str:
    if project in {"PRJNA577971", "PRJNA578171"}:
        return "PRJNA577971_PRJNA578171"
    return project


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=fields, delimiter="\t", extrasaction="ignore"
        )
        writer.writeheader()
        writer.writerows(rows)


def split_fastq_field(value: str) -> list[str]:
    return [item for item in value.split(";") if item]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    all_rows: list[dict[str, str]] = []
    for project in PROJECTS:
        rows = ena_rows(project)
        for row in rows:
            marker, biological, method = marker_and_biosample(row)
            row["study_unit"] = study_unit(project)
            row["marker"] = marker
            row["biological_sample"] = biological
            row["pairing_method"] = method
        all_rows.extend(rows)

    run_fields = [
        "study_unit",
        "source_project",
        "biological_sample",
        "marker",
        "pairing_method",
    ] + RUN_FIELDS
    write_tsv(OUT / "external_candidate_runs.tsv", all_rows, run_fields)

    grouped: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in all_rows:
        grouped[(row["study_unit"], row["biological_sample"])].append(row)

    pair_rows: list[dict[str, object]] = []
    problems: list[str] = []
    for (study, biological), rows in sorted(grouped.items()):
        by_marker = defaultdict(list)
        for row in rows:
            by_marker[row["marker"]].append(row)
        if len(by_marker["16S"]) != 1 or len(by_marker["ITS"]) != 1:
            problems.append(
                f"{study}/{biological}: "
                f"16S={len(by_marker['16S'])}, ITS={len(by_marker['ITS'])}"
            )
            continue
        b = by_marker["16S"][0]
        f = by_marker["ITS"][0]
        pair_rows.append(
            {
                "study": study,
                "biological_sample": biological,
                "run_16S": b["run_accession"],
                "run_ITS": f["run_accession"],
                "sample_accession_16S": b["sample_accession"],
                "sample_accession_ITS": f["sample_accession"],
                "layout_16S": b["library_layout"],
                "layout_ITS": f["library_layout"],
                "instrument_16S": b["instrument_model"],
                "instrument_ITS": f["instrument_model"],
                "reads_16S_raw": b["read_count"],
                "reads_ITS_raw": f["read_count"],
                "bases_16S_raw": b["base_count"],
                "bases_ITS_raw": f["base_count"],
                "fastq_ftp_16S": b["fastq_ftp"],
                "fastq_ftp_ITS": f["fastq_ftp"],
                "pairing_method": b["pairing_method"],
                "pairing_confidence": "high",
            }
        )

    if problems:
        raise RuntimeError("Pairing failures:\n" + "\n".join(problems))

    pair_fields = list(pair_rows[0])
    write_tsv(OUT / "external_16S_ITS_pairing.tsv", pair_rows, pair_fields)

    counts = Counter(row["study"] for row in pair_rows)
    summary_rows: list[dict[str, object]] = []
    for study in EXPECTED_PAIRS:
        project_rows = [row for row in all_rows if row["study_unit"] == study]
        pairs = counts[study]
        expected = EXPECTED_PAIRS[study]
        total_fastq_bytes = sum(
            int(value)
            for row in project_rows
            for value in split_fastq_field(row["fastq_bytes"])
        )
        summary_rows.append(
            {
                "study": study,
                "source_projects": ";".join(
                    sorted({row["source_project"] for row in project_rows})
                ),
                "runs": len(project_rows),
                "paired_biological_samples": pairs,
                "expected_pairs": expected,
                "pairing_complete": pairs == expected,
                "markers": ";".join(
                    sorted({row["marker"] for row in project_rows})
                ),
                "layout_16S": ";".join(
                    sorted(
                        {
                            row["library_layout"]
                            for row in project_rows
                            if row["marker"] == "16S"
                        }
                    )
                ),
                "layout_ITS": ";".join(
                    sorted(
                        {
                            row["library_layout"]
                            for row in project_rows
                            if row["marker"] == "ITS"
                        }
                    )
                ),
                "fastq_size_GiB": round(total_fastq_bytes / 1024**3, 3),
                "recommended_role": (
                    "pH_meta_analysis_and_external_validation"
                    if study == "PRJNA835240"
                    else "core_and_transportability_external_validation"
                ),
            }
        )
    write_tsv(
        OUT / "external_candidate_summary.tsv",
        summary_rows,
        list(summary_rows[0]),
    )

    urls_dir = OUT / "fastq_urls"
    urls_dir.mkdir(exist_ok=True)
    for study in EXPECTED_PAIRS:
        rows = [row for row in all_rows if row["study_unit"] == study]
        urls: list[str] = []
        md5_rows: list[str] = []
        for row in rows:
            ftp = split_fastq_field(row["fastq_ftp"])
            md5 = split_fastq_field(row["fastq_md5"])
            urls.extend(
                "https://" + item.removeprefix("ftp://").removeprefix("https://")
                for item in ftp
            )
            if len(ftp) == len(md5):
                md5_rows.extend(
                    f"{checksum}  {Path(url).name}"
                    for url, checksum in zip(ftp, md5)
                )
        (urls_dir / f"{study}_fastq_urls.txt").write_text(
            "\n".join(urls) + "\n", encoding="utf-8"
        )
        (urls_dir / f"{study}_fastq_md5.txt").write_text(
            "\n".join(md5_rows) + ("\n" if md5_rows else ""),
            encoding="utf-8",
        )

    total_pairs = len(pair_rows)
    if counts != Counter(EXPECTED_PAIRS):
        raise AssertionError(f"Unexpected pair counts: {dict(counts)}")

    print("External cohort audit complete")
    for row in summary_rows:
        print(
            row["study"],
            "runs=", row["runs"],
            "pairs=", row["paired_biological_samples"],
            "layout=", f"{row['layout_16S']}/{row['layout_ITS']}",
        )
    print("New paired samples before QC:", total_pairs)
    print("Projected total before QC:", 211 + total_pairs)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
