# Design notes

This document records current design assumptions for `microC-pipeline` as the repository moves through the v0.4.0 config-driven single-sample milestone.

## Current implementation status

- `bin/microc-pipeline` provides a lightweight config-driven single-sample Micro-C entrypoint.
- `mdp.sh` remains the retained legacy/minimal Slurm script for backward compatibility.
- The config-driven runner removes the legacy need to edit hard-coded FASTQ names by accepting explicit `fastq.r1` and `fastq.r2` paths in YAML.
- The current runnable workflows remain single-sample oriented.
- Restartable chunk-based execution is not available in the current repository.
- Multi-sample sample-sheet orchestration is not implemented.
- Snakemake and Nextflow engines are not implemented.

## v0.4.0 design boundary

The v0.4.0 implementation is intentionally small. It validates one YAML config, creates one per-sample output directory, derives or copies chromosome sizes, writes structured metadata, and executes the same core Micro-C preprocessing path as the legacy script.

The default workflow is Micro-C-oriented and does not require restriction enzyme information. Restriction-fragment-aware Hi-C behavior may be considered in a future milestone, but v0.4.0 does not implement restriction fragment generation, `pairtools restrict`, or enzyme-aware processing.

## Engine direction planning note

No final engine choice has been made for a future production-facing workflow.

Candidate future directions include:

- **Bash-first**: keep the implementation close to the existing script while adding stronger structure, validation, and restart behavior.
- **Snakemake**: use a rule-based workflow manager with explicit inputs, outputs, and cluster execution support.
- **Nextflow**: use a process-based workflow manager with strong portability patterns for local, HPC, and containerized execution.

HiC-Nap restartable chunk-processing concepts are expected to inform future development, especially FASTQ chunking, per-chunk status tracking, conservative restart behavior, global merge/deduplication, and progress logging. Those concepts are planning inputs only in this repository today; they are not implemented here yet.

## Near-term design boundaries

Near-term work should avoid turning v0.4.0 into a large framework. Later milestones may standardize valid pairs and matrix outputs, add restartable chunk execution, choose a workflow engine, add sample-sheet support, improve QC, and introduce packaging or containers.
