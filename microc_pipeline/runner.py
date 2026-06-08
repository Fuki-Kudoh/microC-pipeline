"""Command construction and execution for the v0.5.0 single-sample workflow."""

from __future__ import annotations

import json
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from shlex import quote

from . import __version__
from .config import ConfigError, PipelineConfig
from .outputs import build_manifest, final_output_paths, final_outputs
from .validation import OutputValidationError, validate_outputs, validation_checks

OUTPUT_SUBDIRS = ("bam", "cool", "fastqc", "hic", "genome", "logs", "pairs", "qc", "stats", "temp")
REQUIRED_COMMANDS = (
    "fastqc",
    "trim_galore",
    "samtools",
    "bwa",
    "pairtools",
    "preseq",
    "bgzip",
    "pairix",
    "cooler",
    "python3",
)


class RunnerError(RuntimeError):
    """Raised when execution cannot continue."""


def q(value: object) -> str:
    return quote(str(value))


def trim_galore_output_name(fastq: Path, read_number: int) -> str:
    name = fastq.name
    for suffix in (".fastq.gz", ".fq.gz", ".fastq", ".fq"):
        if name.endswith(suffix):
            name = name[: -len(suffix)]
            break
    return f"{name}_val_{read_number}.fq.gz"


def check_dependencies(config: PipelineConfig) -> None:
    required = list(REQUIRED_COMMANDS)
    if config.outputs.make_hic:
        required.append("juicer_tools")
    missing = [command for command in required if shutil.which(command) is None]
    if missing:
        raise RunnerError("Missing required command(s): " + ", ".join(missing))


def ensure_output_dirs(config: PipelineConfig) -> dict[str, Path]:
    sample_dir = config.sample_output_dir
    dirs = {name: sample_dir / name for name in OUTPUT_SUBDIRS}
    sample_dir.mkdir(parents=True, exist_ok=True)
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def derive_chrom_sizes(config: PipelineConfig, dirs: dict[str, Path], *, dry_run: bool) -> tuple[Path, list[str]]:
    target = dirs["genome"] / f"{config.genome.name}.chrom.sizes"
    commands: list[str] = []
    if config.genome.chrom_sizes is not None:
        commands.append(f"cp {q(config.genome.chrom_sizes)} {q(target)}")
        if not dry_run:
            shutil.copyfile(config.genome.chrom_sizes, target)
        return target, commands

    fai = Path(str(config.genome.fasta) + ".fai")
    if not fai.is_file() or fai.stat().st_size == 0:
        fasta_dir = config.genome.fasta.parent
        if not dry_run and not fasta_dir.exists():
            raise RunnerError(f"FASTA directory does not exist: {fasta_dir}")
        if not dry_run and not fasta_dir.is_dir():
            raise RunnerError(f"FASTA directory is not a directory: {fasta_dir}")
        if not dry_run and not fasta_dir.stat().st_mode:
            raise RunnerError(f"Cannot inspect FASTA directory: {fasta_dir}")
        if dry_run or os_access_writable(fasta_dir):
            commands.append(f"samtools faidx {q(config.genome.fasta)}")
            if not dry_run:
                run_shell(commands[-1])
        else:
            raise RunnerError(
                f"FASTA index not found or empty: {fai}. Create it first with: samtools faidx {q(config.genome.fasta)}"
            )
    commands.append(f"cut -f1,2 {q(fai)} > {q(target)}")
    if not dry_run:
        with open(fai, "r", encoding="utf-8") as source, open(target, "w", encoding="utf-8") as sink:
            for line in source:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 2:
                    sink.write(f"{fields[0]}\t{fields[1]}\n")
    return target, commands


def os_access_writable(path: Path) -> bool:
    import os

    return os.access(path, os.W_OK)


def build_metadata(config: PipelineConfig, chrom_sizes: Path, commands: list[str], *, dry_run: bool) -> dict[str, object]:
    return {
        "sample": config.sample,
        "assay": config.assay,
        "config_path": str(config.config_path),
        "output_dir": str(config.sample_output_dir),
        "genome_name": config.genome.name,
        "genome_fasta": str(config.genome.fasta),
        "chrom_sizes": str(chrom_sizes),
        "fastq_r1": str(config.fastq_r1),
        "fastq_r2": str(config.fastq_r2),
        "threads": config.threads,
        "bin_sizes": config.bin_sizes,
        "outputs": {
            "keep_bam": config.outputs.keep_bam,
            "make_hic": config.outputs.make_hic,
            "make_mcool": config.outputs.make_mcool,
        },
        "final_outputs": final_output_paths(config),
        "pipeline_version": __version__,
        "dry_run": dry_run,
        "created_at": datetime.now(timezone.utc).isoformat(),
        "commands": commands,
    }


def workflow_commands(config: PipelineConfig, dirs: dict[str, Path], chrom_sizes: Path) -> list[str]:
    sample = config.sample
    threads = config.threads
    first_bin_size = config.bin_sizes[0]
    repo_root = Path(__file__).resolve().parents[1]
    get_qc = repo_root / "get_qc.py"

    trim_dir = dirs["temp"] / "trimmed_fastq"
    sample_temp = dirs["temp"] / sample
    trimmed_r1 = trim_dir / trim_galore_output_name(config.fastq_r1, 1)
    trimmed_r2 = trim_dir / trim_galore_output_name(config.fastq_r2, 2)
    sam = sample_temp / f"{sample}.sam"
    pairsam = sample_temp / f"{sample}.pairsam"
    sorted_pairsam = sample_temp / f"{sample}.sorted.pairsam"
    dedup_pairsam = sample_temp / f"{sample}.dedup.pairsam"
    undedup_bam = dirs["bam"] / f"{sample}.undedup.bam"
    undedup_pt_bam = dirs["bam"] / f"{sample}.undedup.PT.bam"
    raw_bam = dirs["bam"] / f"{sample}.bam"
    final_bam = dirs["bam"] / f"{sample}.PT.bam"
    final_bai = Path(str(final_bam) + ".bai")
    outputs = final_outputs(config)
    pairs = outputs.valid_pairs_uncompressed
    pairs_gz = outputs.valid_pairs
    cool = outputs.cool

    commands = [
        f"mkdir -p {q(trim_dir)} {q(sample_temp)}",
        f"fastqc -t {threads} -o {q(dirs['fastqc'])} {q(config.fastq_r1)} {q(config.fastq_r2)}",
        f"trim_galore -j {threads} -o {q(trim_dir)} --paired {q(config.fastq_r1)} {q(config.fastq_r2)}",
        f"bwa mem -5SP -T0 -t {threads} {q(config.genome.fasta)} {q(trimmed_r1)} {q(trimmed_r2)} -o {q(sam)}",
        f"pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 --nproc-in {threads} --nproc-out {threads} --chroms-path {q(chrom_sizes)} {q(sam)} > {q(pairsam)}",
        f"pairtools sort --nproc {threads} --tmpdir={q(sample_temp)} {q(pairsam)} > {q(sorted_pairsam)}",
        f"rm {q(pairsam)}",
        f"pairtools split --nproc-in {threads} --nproc-out {threads} --output-sam {q(undedup_bam)} {q(sorted_pairsam)}",
        f"samtools sort -@ {threads} -T {q(undedup_bam)} -o {q(undedup_pt_bam)} {q(undedup_bam)}",
        f"rm {q(undedup_bam)}",
        f"samtools index {q(undedup_pt_bam)}",
        f"preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 -output {q(outputs.preseq)} {q(undedup_pt_bam)}",
        f"rm {q(undedup_pt_bam)}",
        f"pairtools dedup --nproc-in {threads} --nproc-out {threads} --mark-dups --output-stats {q(outputs.pairtools_stats)} --output {q(dedup_pairsam)} {q(sorted_pairsam)}",
        f"rm {q(sorted_pairsam)}",
        f"python3 {q(get_qc)} -p {q(outputs.pairtools_stats)} --format tsv > {q(outputs.qc_tsv)}",
        f"pairtools split --nproc-in {threads} --nproc-out {threads} --output-pairs {q(pairs)} --output-sam {q(raw_bam)} {q(dedup_pairsam)}",
        f"rm {q(dedup_pairsam)}",
        f"samtools sort -@ {threads} -T {q(raw_bam)} -o {q(final_bam)} {q(raw_bam)}",
        f"rm {q(raw_bam)}",
        f"samtools index {q(final_bam)}",
    ]
    if config.outputs.make_hic:
        commands.append(f"juicer_tools pre -j {threads} {q(pairs)} {q(outputs.hic)} {q(config.genome.name)}")
    commands.extend(
        [
            f"bgzip {q(pairs)}",
            f"pairix {q(pairs_gz)}",
            f"cooler cload pairix -p {threads} {q(str(chrom_sizes) + ':' + str(first_bin_size))} {q(pairs_gz)} {q(cool)}",
        ]
    )
    if config.outputs.make_mcool:
        commands.append(f"cooler zoomify -p {threads} -o {q(outputs.mcool)} {q(cool)}")
    if not config.outputs.keep_bam:
        commands.append(f"rm -f {q(final_bam)} {q(final_bai)}")
    return commands


def run_shell(command: str, log_handle=None) -> None:
    if log_handle is not None:
        print(f"$ {command}", file=log_handle, flush=True)
    completed = subprocess.run(command, shell=True, stdout=log_handle, stderr=subprocess.STDOUT)
    if completed.returncode != 0:
        raise RunnerError(f"Command failed with exit code {completed.returncode}: {command}")


def print_plan(config: PipelineConfig, chrom_commands: list[str], commands: list[str]) -> None:
    print(f"Sample: {config.sample}")
    print(f"Assay: {config.assay}")
    print(f"Output directory: {config.sample_output_dir}")
    print(f"Threads: {config.threads}")
    if len(config.bin_sizes) > 1:
        print(
            "Warning: v0.5.0 uses only the first configured bin size for .cool generation: "
            f"{config.bin_sizes[0]}"
        )
    print("Planned commands:")
    for command in [*chrom_commands, *commands]:
        print(command)
    print("Expected final outputs:")
    for name, path in final_output_paths(config).items():
        print(f"{name}: {path}")
    print("Planned validation checks:")
    for check in validation_checks(config):
        print(f"- {check}")


def write_output_manifest(config: PipelineConfig, validation_results: dict[str, dict[str, object]] | None = None) -> None:
    manifest_path = final_outputs(config).output_manifest
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(build_manifest(config, validation_results), indent=2) + "\n",
        encoding="utf-8",
    )


def validate_and_write_manifest(config: PipelineConfig) -> None:
    try:
        results = validate_outputs(config)
    except OutputValidationError as exc:
        write_output_manifest(config, exc.results)
        raise RunnerError(str(exc)) from exc
    write_output_manifest(config, results)


def run_pipeline(config: PipelineConfig, *, dry_run: bool = False) -> None:
    dry_dirs = {name: config.sample_output_dir / name for name in OUTPUT_SUBDIRS}
    if dry_run:
        dirs = dry_dirs
    else:
        check_dependencies(config)
        dirs = ensure_output_dirs(config)
    chrom_sizes, chrom_commands = derive_chrom_sizes(config, dirs, dry_run=dry_run)
    commands = workflow_commands(config, dirs, chrom_sizes)
    all_commands = [*chrom_commands, *commands]

    if dry_run:
        print_plan(config, chrom_commands, commands)
        return

    metadata = build_metadata(config, chrom_sizes, all_commands, dry_run=False)
    metadata_path = final_outputs(config).run_metadata
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")

    log_path = dirs["logs"] / f"{config.sample}.log"
    with open(log_path, "a", encoding="utf-8") as log_handle:
        print(f"Sample: {config.sample}", file=log_handle)
        print(f"Assay: {config.assay}", file=log_handle)
        print(f"Config path: {config.config_path}", file=log_handle)
        print(f"Output directory: {config.sample_output_dir}", file=log_handle)
        print(f"Threads: {config.threads}", file=log_handle)
        if len(config.bin_sizes) > 1:
            print(
                "Warning: v0.5.0 uses only the first configured bin size for .cool generation: "
                f"{config.bin_sizes[0]}",
                file=log_handle,
            )
        for command in chrom_commands:
            print(f"$ {command}", file=log_handle, flush=True)
        for command in commands:
            run_shell(command, log_handle=log_handle)

    validate_and_write_manifest(config)


def validate_chrom_sizes_feasibility(config: PipelineConfig) -> None:
    if config.genome.chrom_sizes is not None:
        return
    fai = Path(str(config.genome.fasta) + ".fai")
    if fai.is_file() and fai.stat().st_size > 0:
        return
    if not os_access_writable(config.genome.fasta.parent):
        raise ConfigError(
            f"FASTA index not found or empty: {fai}; FASTA directory is not writable. "
            f"Create the index first with: samtools faidx {q(config.genome.fasta)}"
        )
