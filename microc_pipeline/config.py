"""Configuration loading and validation for the v0.5.1 single-sample runner."""

from __future__ import annotations

import os
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any


class ConfigError(ValueError):
    """Raised when a configuration file is invalid."""


@dataclass(frozen=True)
class GenomeConfig:
    name: str
    fasta: Path
    chrom_sizes: Path | None


@dataclass(frozen=True)
class OutputConfig:
    keep_bam: bool = True
    make_hic: bool = True
    make_mcool: bool = True


@dataclass(frozen=True)
class PipelineConfig:
    sample: str
    assay: str
    fastq_r1: Path
    fastq_r2: Path
    genome: GenomeConfig
    threads: int
    bin_sizes: list[int]
    output_root: Path
    sample_output_dir: Path
    outputs: OutputConfig
    config_path: Path
    raw: dict[str, Any]


_MISSING = object()
_SAFE_SAMPLE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]*$")


def _strip_comment(line: str) -> str:
    in_single = False
    in_double = False
    escaped = False
    result = []
    for char in line:
        if char == "'" and not in_double:
            in_single = not in_single
        elif char == '"' and not in_single and not escaped:
            in_double = not in_double
        elif char == "#" and not in_single and not in_double:
            break
        result.append(char)
        escaped = char == "\\" and not escaped
    return "".join(result).rstrip()


def _parse_scalar(value: str) -> Any:
    value = value.strip()
    if value in {"", "null", "Null", "NULL", "~"}:
        return None
    if value in {"true", "True", "TRUE"}:
        return True
    if value in {"false", "False", "FALSE"}:
        return False
    if (value.startswith('"') and value.endswith('"')) or (
        value.startswith("'") and value.endswith("'")
    ):
        return value[1:-1]
    try:
        return int(value)
    except ValueError:
        return value


def _parse_simple_yaml(text: str) -> dict[str, Any]:
    """Parse the small YAML subset used by the documented v0.5.0 config.

    This intentionally supports mappings, nested mappings, scalar values, and
    lists of scalar values so the runner does not need an external Python YAML
    dependency on minimal HPC login nodes. It is not a general YAML parser.
    """

    root: dict[str, Any] = {}
    stack: list[tuple[int, Any]] = [(-1, root)]
    lines = text.splitlines()

    for line_number, original_line in enumerate(lines, start=1):
        line = _strip_comment(original_line)
        if not line.strip():
            continue
        indent = len(line) - len(line.lstrip(" "))
        if indent % 2 != 0:
            raise ConfigError(f"Invalid indentation at line {line_number}: use multiples of two spaces.")
        stripped = line.strip()

        while stack and indent <= stack[-1][0]:
            stack.pop()
        if not stack:
            raise ConfigError(f"Invalid indentation at line {line_number}.")
        parent = stack[-1][1]

        if stripped.startswith("- "):
            if not isinstance(parent, list):
                raise ConfigError(f"Unexpected list item at line {line_number}.")
            parent.append(_parse_scalar(stripped[2:]))
            continue

        if ":" not in stripped:
            raise ConfigError(f"Expected key/value pair at line {line_number}.")
        key, value = stripped.split(":", 1)
        key = key.strip()
        if not key:
            raise ConfigError(f"Empty key at line {line_number}.")
        if not isinstance(parent, dict):
            raise ConfigError(f"Cannot add key under list at line {line_number}.")

        value = value.strip()
        if value == "":
            next_container: Any = {}
            for following in lines[line_number:]:
                following_clean = _strip_comment(following)
                if not following_clean.strip():
                    continue
                following_indent = len(following_clean) - len(following_clean.lstrip(" "))
                following_stripped = following_clean.strip()
                if following_indent <= indent:
                    break
                if following_stripped.startswith("- "):
                    next_container = []
                break
            parent[key] = next_container
            stack.append((indent, next_container))
        else:
            parent[key] = _parse_scalar(value)

    return root


def load_config_file(path: str | os.PathLike[str]) -> dict[str, Any]:
    config_path = Path(path)
    if not config_path.exists():
        raise ConfigError(f"Config file not found: {config_path}")
    try:
        return _parse_simple_yaml(config_path.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ConfigError(f"Could not read config file: {config_path}: {exc}") from exc


def _get_required(mapping: dict[str, Any], dotted: str) -> Any:
    current: Any = mapping
    for part in dotted.split("."):
        if not isinstance(current, dict) or part not in current or current[part] in (None, ""):
            raise ConfigError(f"Missing required config field: {dotted}")
        current = current[part]
    return current


def _get_optional(mapping: dict[str, Any], dotted: str, default: Any = None) -> Any:
    current: Any = mapping
    for part in dotted.split("."):
        if not isinstance(current, dict) or part not in current:
            return default
        current = current[part]
    return current


def _require_string(value: Any, field: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"Invalid config field {field}: expected a non-empty string.")
    return value


def _validate_sample_name(sample: str) -> None:
    if not _SAFE_SAMPLE_RE.fullmatch(sample):
        raise ConfigError(
            "Invalid config field sample: expected a safe sample name matching "
            r"^[A-Za-z0-9][A-Za-z0-9._-]*$ "
            "(start with a letter or digit; use only letters, digits, dot, underscore, or hyphen)."
        )


def _require_bool(value: Any, field: str) -> bool:
    if not isinstance(value, bool):
        raise ConfigError(f"Invalid config field {field}: expected true or false.")
    return value


def _resolve_path(value: Any, field: str) -> Path:
    return Path(_require_string(value, field)).expanduser()


def _validate_existing_file(path: Path, field: str) -> None:
    if not path.is_file():
        raise ConfigError(f"{field} not found: {path}")


def parse_and_validate_config(
    config_path: str | os.PathLike[str], *, check_files: bool = True
) -> PipelineConfig:
    raw = load_config_file(config_path)
    if not isinstance(raw, dict):
        raise ConfigError("Config root must be a mapping.")

    sample = _require_string(_get_required(raw, "sample"), "sample")
    _validate_sample_name(sample)
    assay = _require_string(_get_optional(raw, "assay", "microc"), "assay")
    if assay != "microc":
        raise ConfigError(f"Unsupported assay: {assay}. v0.5.0 supports only assay: microc.")

    fastq_r1 = _resolve_path(_get_required(raw, "fastq.r1"), "fastq.r1")
    fastq_r2 = _resolve_path(_get_required(raw, "fastq.r2"), "fastq.r2")
    genome_name = _require_string(_get_required(raw, "genome.name"), "genome.name")
    genome_fasta = _resolve_path(_get_required(raw, "genome.fasta"), "genome.fasta")
    chrom_sizes_value = _get_optional(raw, "genome.chrom_sizes", None)
    chrom_sizes = None if chrom_sizes_value in (None, "") else _resolve_path(chrom_sizes_value, "genome.chrom_sizes")
    output_root = _resolve_path(_get_required(raw, "output_dir"), "output_dir")

    threads_value = _get_optional(raw, "threads", None)
    if threads_value is None:
        threads_value = os.environ.get("SLURM_CPUS_PER_TASK", 16)
    try:
        threads = int(threads_value)
    except (TypeError, ValueError) as exc:
        raise ConfigError("Invalid config field threads: expected a positive integer.") from exc
    if threads < 1:
        raise ConfigError("Invalid config field threads: expected a positive integer.")

    bin_sizes_value = _get_optional(raw, "bin_sizes", [1000])
    if not isinstance(bin_sizes_value, list) or not bin_sizes_value:
        raise ConfigError("Invalid config field bin_sizes: expected a non-empty list of positive integers.")
    bin_sizes: list[int] = []
    for index, bin_size in enumerate(bin_sizes_value):
        try:
            parsed = int(bin_size)
        except (TypeError, ValueError) as exc:
            raise ConfigError(f"Invalid config field bin_sizes[{index}]: expected a positive integer.") from exc
        if parsed < 1:
            raise ConfigError(f"Invalid config field bin_sizes[{index}]: expected a positive integer.")
        bin_sizes.append(parsed)

    outputs = OutputConfig(
        keep_bam=_require_bool(_get_optional(raw, "outputs.keep_bam", True), "outputs.keep_bam"),
        make_hic=_require_bool(_get_optional(raw, "outputs.make_hic", True), "outputs.make_hic"),
        make_mcool=_require_bool(_get_optional(raw, "outputs.make_mcool", True), "outputs.make_mcool"),
    )

    if check_files:
        _validate_existing_file(fastq_r1, "FASTQ R1")
        _validate_existing_file(fastq_r2, "FASTQ R2")
        _validate_existing_file(genome_fasta, "Genome FASTA")
        if chrom_sizes is not None:
            _validate_existing_file(chrom_sizes, "Genome chromosome sizes")

    return PipelineConfig(
        sample=sample,
        assay=assay,
        fastq_r1=fastq_r1,
        fastq_r2=fastq_r2,
        genome=GenomeConfig(name=genome_name, fasta=genome_fasta, chrom_sizes=chrom_sizes),
        threads=threads,
        bin_sizes=bin_sizes,
        output_root=output_root,
        sample_output_dir=output_root / sample,
        outputs=outputs,
        config_path=Path(config_path),
        raw=raw,
    )
