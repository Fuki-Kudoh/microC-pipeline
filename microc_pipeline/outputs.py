"""Standardized v0.5.0 output paths and manifest helpers."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

from . import __version__
from .config import PipelineConfig

HIC_FORMAT_NOTE = "Juicer .hic contact-map file format; not enzyme-aware Hi-C assay mode."


@dataclass(frozen=True)
class FinalOutputs:
    """Standardized final output paths for a single configured sample."""

    sample_dir: Path
    valid_pairs: Path
    valid_pairs_index: Path
    cool: Path
    mcool: Path
    hic: Path
    pairtools_stats: Path
    preseq: Path
    qc_tsv: Path
    bam: Path
    bam_index: Path
    output_manifest: Path
    run_metadata: Path

    @property
    def valid_pairs_uncompressed(self) -> Path:
        return self.valid_pairs.with_suffix("")


def final_outputs(config: PipelineConfig) -> FinalOutputs:
    """Return standardized final output paths for *config*."""

    sample = config.sample
    sample_dir = config.sample_output_dir
    return FinalOutputs(
        sample_dir=sample_dir,
        valid_pairs=sample_dir / "pairs" / f"{sample}.valid.pairs.gz",
        valid_pairs_index=sample_dir / "pairs" / f"{sample}.valid.pairs.gz.px2",
        cool=sample_dir / "cool" / f"{sample}.cool",
        mcool=sample_dir / "cool" / f"{sample}.mcool",
        hic=sample_dir / "hic" / f"{sample}.hic",
        pairtools_stats=sample_dir / "stats" / f"{sample}.pairtools.stats.txt",
        preseq=sample_dir / "stats" / f"{sample}.preseq.lc_extrap.txt",
        qc_tsv=sample_dir / "qc" / f"{sample}.qc.tsv",
        bam=sample_dir / "bam" / f"{sample}.PT.bam",
        bam_index=sample_dir / "bam" / f"{sample}.PT.bam.bai",
        output_manifest=sample_dir / "output_manifest.json",
        run_metadata=sample_dir / "run_metadata.json",
    )


def final_output_paths(config: PipelineConfig) -> dict[str, str]:
    """Return the run_metadata final_outputs mapping for expected outputs."""

    outputs = final_outputs(config)
    paths = {
        "valid_pairs": str(outputs.valid_pairs),
        "valid_pairs_index": str(outputs.valid_pairs_index),
        "cool": str(outputs.cool),
        "pairtools_stats": str(outputs.pairtools_stats),
        "preseq": str(outputs.preseq),
        "qc_tsv": str(outputs.qc_tsv),
        "output_manifest": str(outputs.output_manifest),
    }
    if config.outputs.make_mcool:
        paths["mcool"] = str(outputs.mcool)
    if config.outputs.make_hic:
        paths["hic"] = str(outputs.hic)
    if config.outputs.keep_bam:
        paths["bam"] = str(outputs.bam)
        paths["bam_index"] = str(outputs.bam_index)
    return paths


def expected_output_entries(config: PipelineConfig) -> dict[str, dict[str, Any]]:
    """Return manifest output entries expected for *config*."""

    outputs = final_outputs(config)
    entries: dict[str, dict[str, Any]] = {
        "valid_pairs": {
            "path": str(outputs.valid_pairs),
            "index": str(outputs.valid_pairs_index),
            "required": True,
        },
        "cool": {"path": str(outputs.cool), "required": True},
        "pairtools_stats": {"path": str(outputs.pairtools_stats), "required": True},
        "preseq": {"path": str(outputs.preseq), "required": True},
        "qc_tsv": {"path": str(outputs.qc_tsv), "required": True},
    }
    if config.outputs.make_mcool:
        entries["mcool"] = {"path": str(outputs.mcool), "required": True}
    else:
        entries["mcool"] = {"path": str(outputs.mcool), "required": False, "expected": False}
    if config.outputs.make_hic:
        entries["hic"] = {
            "path": str(outputs.hic),
            "required": True,
            "note": HIC_FORMAT_NOTE,
        }
    else:
        entries["hic"] = {
            "path": str(outputs.hic),
            "required": False,
            "expected": False,
            "note": HIC_FORMAT_NOTE,
        }
    if config.outputs.keep_bam:
        entries["bam"] = {
            "path": str(outputs.bam),
            "index": str(outputs.bam_index),
            "required": True,
        }
    else:
        entries["bam"] = {
            "path": str(outputs.bam),
            "index": str(outputs.bam_index),
            "required": False,
            "expected": False,
        }
    return entries


def build_manifest(config: PipelineConfig, validation: dict[str, dict[str, Any]] | None = None) -> dict[str, Any]:
    """Build an output manifest, optionally merging validation results."""

    entries = expected_output_entries(config)
    validation = validation or {}
    for name, result in validation.items():
        if name in entries:
            entries[name].update(result)
    for entry in entries.values():
        entry.setdefault("validated", False)
    return {
        "sample": config.sample,
        "assay": config.assay,
        "output_dir": str(config.sample_output_dir),
        "pipeline_version": __version__,
        "outputs": entries,
    }
