"""Lightweight validation for standardized v0.5.0 final outputs."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from .config import PipelineConfig
from .outputs import expected_output_entries

PAIRTOOLS_STATS_REQUIRED_KEYS = (
    "total",
    "total_unmapped",
    "total_mapped",
    "total_dups",
    "total_nodups",
    "cis",
    "trans",
    "cis_1kb+",
    "cis_10kb+",
)

QC_TSV_REQUIRED_COLUMNS = ("metric", "value", "percent")


class OutputValidationError(RuntimeError):
    """Raised when one or more expected final outputs are missing or invalid."""

    def __init__(self, failures: list[str], results: dict[str, dict[str, Any]]):
        self.failures = failures
        self.results = results
        super().__init__("Output validation failed:\n" + "\n".join(f"- {failure}" for failure in failures))


def _file_status(path: Path) -> tuple[bool, int]:
    if not path.is_file():
        return False, 0
    return True, path.stat().st_size


def _check_nonempty(path: Path, label: str, failures: list[str]) -> dict[str, Any]:
    exists, size = _file_status(path)
    if not exists:
        failures.append(f"{label} is missing: {path}")
    elif size <= 0:
        failures.append(f"{label} is empty: {path}")
    return {"exists": exists, "size_bytes": size}


def _read_stats_keys(path: Path) -> set[str]:
    keys: set[str] = set()
    with open(path, "r", encoding="utf-8") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) >= 2:
                keys.add(fields[0])
    return keys


def validation_checks(config: PipelineConfig) -> list[str]:
    """Describe planned validation checks for dry-run output."""

    entries = expected_output_entries(config)
    checks = [
        "valid_pairs: non-empty .valid.pairs.gz and non-empty .px2 pairix index",
        "cool: non-empty .cool file",
        "pairtools_stats: non-empty stats file with required get_qc.py keys",
        "preseq: non-empty lc_extrap output",
        "qc_tsv: non-empty TSV with metric/value/percent header columns",
    ]
    if entries["mcool"].get("required"):
        checks.append("mcool: non-empty .mcool file")
    if entries["hic"].get("required"):
        checks.append("hic: non-empty Juicer .hic contact-map file")
    if entries["bam"].get("required"):
        checks.append("bam: non-empty final BAM and BAI files")
    return checks


def validate_outputs(config: PipelineConfig) -> dict[str, dict[str, Any]]:
    """Validate expected final outputs for *config*.

    Returns per-output validation metadata. Raises OutputValidationError when a
    required output is missing or invalid.
    """

    entries = expected_output_entries(config)
    failures: list[str] = []
    results: dict[str, dict[str, Any]] = {}

    for name, entry in entries.items():
        if not entry.get("required", False):
            results[name] = {"validated": True, "skipped": True, "reason": "not expected by config"}
            continue

        output_failures: list[str] = []
        path = Path(entry["path"])
        result: dict[str, Any] = {"checks": {}}
        result["checks"]["path"] = _check_nonempty(path, name, output_failures)

        if "index" in entry:
            result["checks"]["index"] = _check_nonempty(Path(entry["index"]), f"{name} index", output_failures)

        if name == "pairtools_stats" and path.is_file() and path.stat().st_size > 0:
            keys = _read_stats_keys(path)
            missing = [key for key in PAIRTOOLS_STATS_REQUIRED_KEYS if key not in keys]
            result["required_keys"] = list(PAIRTOOLS_STATS_REQUIRED_KEYS)
            if missing:
                output_failures.append(
                    "pairtools_stats is missing required key(s): " + ", ".join(missing)
                )

        if name == "qc_tsv" and path.is_file() and path.stat().st_size > 0:
            with open(path, "r", encoding="utf-8") as handle:
                header = handle.readline().rstrip("\n").split("\t")
            result["header"] = header
            missing_columns = [column for column in QC_TSV_REQUIRED_COLUMNS if column not in header]
            if missing_columns:
                output_failures.append("qc_tsv header missing column(s): " + ", ".join(missing_columns))

        result["validated"] = not output_failures
        if output_failures:
            result["errors"] = output_failures
            failures.extend(output_failures)
        results[name] = result

    if failures:
        raise OutputValidationError(failures, results)
    return results
