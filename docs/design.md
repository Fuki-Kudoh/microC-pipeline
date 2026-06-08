# Design notes

This document records current design assumptions for `microC-pipeline` as the repository is prepared for the `v0.3.0` project-structure milestone.

## Current implementation status

- `mdp.sh` remains the retained legacy/minimal Slurm script.
- The current runnable workflow is single-sample oriented and uses historical FASTQ naming assumptions.
- The modern production-facing pipeline engine is not implemented yet.
- Restartable chunk-based execution is not available in the current repository.

## Engine direction planning note

No final engine choice has been made for the future production-facing workflow.

Candidate future directions include:

- **Bash-first**: keep the implementation close to the existing script while adding stronger structure, validation, and restart behavior.
- **Snakemake**: use a rule-based workflow manager with explicit inputs, outputs, and cluster execution support.
- **Nextflow**: use a process-based workflow manager with strong portability patterns for local, HPC, and containerized execution.

HiC-Nap restartable chunk-processing concepts are expected to inform future development, especially FASTQ chunking, per-chunk status tracking, conservative restart behavior, global merge/deduplication, and progress logging. Those concepts are planning inputs only in this repository today; they are not implemented here yet.

## Near-term design boundaries

The `v0.3.0` milestone is limited to project identity and lightweight structure. It should not add the modern engine, sample-sheet execution, workflow-manager files, containers, CI, or full QC reports.
