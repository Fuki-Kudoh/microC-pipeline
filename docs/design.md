# Design notes

## Current design

`microC-pipeline` currently has two interfaces:

1. A retained legacy Slurm script, `mdp.sh`, kept for backward compatibility.
2. A lightweight config-driven single-sample runner, `bin/microc-pipeline`, for modernized Micro-C preprocessing.

The v0.5.0 milestone kept the v0.4.0 single-sample runner model and made final valid pairs, matrix, stats, QC, and manifest files first-class products. The v0.5.1 hardening pass adds stricter preflight checks for sample names and BWA indexes plus best-effort tool version capture in `run_metadata.json`, without changing the core processing commands.

## Micro-C-first boundary

The runner is Micro-C-first:

- no restriction enzyme field is required or accepted as a workflow concept;
- no restriction fragments are generated;
- `pairtools restrict` is not run;
- enzyme-aware Hi-C processing is not implemented.

Optional `.hic` output is only the Juicer contact-map file format. It does not indicate Hi-C assay mode.

## What v0.5.x deliberately does not add

The v0.5.x output and preflight milestones do **not** add:

- multi-sample sample-sheet execution;
- restartable chunk-based execution;
- Snakemake or Nextflow workflows;
- containers;
- CI;
- full HTML QC reports;
- project-level QC summaries;
- enzyme-aware Hi-C behavior;
- restriction fragment generation;
- `pairtools restrict`;
- loop, compartment, TAD, or differential contact analysis;
- downstream biological interpretation.

These remain future milestones.

## Metadata and manifest split

`run_metadata.json` records how a run was configured and executed: config path, inputs, genome resources, threads, commands, output toggles, standardized final output paths, and best-effort tool version information.

`output_manifest.json` records what final outputs were expected and whether lightweight validation passed. Keeping these concerns separate makes output inspection possible without treating execution provenance as validation state.

## Future direction

Future milestones may add restartability, multi-sample orchestration, workflow-manager integrations, richer QC, packaging, tests, and deeper validators. Those changes should preserve the v0.5.0 final output names where possible.
