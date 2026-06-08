"""CLI for the v0.5.0 config-driven single-sample Micro-C runner."""

from __future__ import annotations

import argparse
import sys

from .config import ConfigError, parse_and_validate_config
from .runner import (
    RunnerError,
    run_pipeline,
    validate_and_write_manifest,
    validate_bwa_index,
    validate_chrom_sizes_feasibility,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="microc-pipeline",
        description="Lightweight config-driven single-sample Micro-C preprocessing runner.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    run_parser = subparsers.add_parser("run", help="Run the single-sample Micro-C workflow.")
    run_parser.add_argument("--config", required=True, help="YAML config file for one sample.")
    run_parser.add_argument("--dry-run", action="store_true", help="Print planned commands without running tools.")

    validate_parser = subparsers.add_parser("validate-config", help="Validate a single-sample YAML config.")
    validate_parser.add_argument("--config", required=True, help="YAML config file to validate.")

    validate_outputs_parser = subparsers.add_parser(
        "validate-outputs",
        help="Validate standardized final outputs for a completed single-sample run.",
    )
    validate_outputs_parser.add_argument("--config", required=True, help="YAML config file for the completed run.")
    return parser


def command_validate_config(args: argparse.Namespace) -> int:
    config = parse_and_validate_config(args.config, check_files=True)
    validate_bwa_index(config)
    validate_chrom_sizes_feasibility(config)
    print(f"Config validation succeeded for sample: {config.sample}")
    print(f"Assay: {config.assay}")
    print(f"Output directory: {config.sample_output_dir}")
    return 0


def command_validate_outputs(args: argparse.Namespace) -> int:
    config = parse_and_validate_config(args.config, check_files=False)
    validate_and_write_manifest(config)
    print(f"Output validation succeeded for sample: {config.sample}")
    print(f"Manifest: {config.sample_output_dir / 'output_manifest.json'}")
    return 0


def command_run(args: argparse.Namespace) -> int:
    config = parse_and_validate_config(args.config, check_files=True)
    validate_bwa_index(config)
    validate_chrom_sizes_feasibility(config)
    run_pipeline(config, dry_run=args.dry_run)
    if not args.dry_run:
        print(f"Run completed and outputs validated for sample: {config.sample}")
    return 0


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        if args.command == "validate-config":
            return command_validate_config(args)
        if args.command == "validate-outputs":
            return command_validate_outputs(args)
        if args.command == "run":
            return command_run(args)
        parser.error(f"Unknown command: {args.command}")
    except (ConfigError, RunnerError) as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0
