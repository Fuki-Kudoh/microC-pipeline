# Changelog

All notable changes to this project will be documented in this file.

## v0.5.0 - 2026-06-08

### Added
- Added standardized config-driven final output names for valid pairs, matrices, stats, QC, optional `.hic`, and optional BAM outputs.
- Added `SAMPLE.valid.pairs.gz` as the first-class valid pairs output and `SAMPLE.valid.pairs.gz.px2` as its Pairix index.
- Added explicit `SAMPLE.cool`, optional `SAMPLE.mcool`, and optional `SAMPLE.hic` final output paths.
- Added `SAMPLE.pairtools.stats.txt`, `SAMPLE.preseq.lc_extrap.txt`, and `SAMPLE.qc.tsv` final stats/QC paths.
- Added lightweight output validation for required final products.
- Added `output_manifest.json` to record expected outputs and validation status.
- Added `bin/microc-pipeline validate-outputs` for checking completed runs without rerunning preprocessing.
- Added `get_qc.py --format tsv` while preserving the default text output.

### Changed
- Updated `run_metadata.json` content to include standardized final output paths and the output manifest path.
- Updated README and documentation for v0.5.0 output semantics, validation, manifest behavior, and `.hic` file-format wording.

### Notes
- `.hic` is documented as the Juicer contact-map file format, not enzyme-aware Hi-C assay support.
- Multi-sample execution, restartable chunk processing, workflow-manager files, containers, CI, full HTML QC reports, and enzyme-aware Hi-C behavior remain future work.

## v0.4.0 - 2026-06-08

### Added
- Added a config-driven single-sample Micro-C runner at `bin/microc-pipeline`.
- Added `config/example.single-sample.yaml` as an example single-sample YAML config without restriction enzyme settings.
- Added config validation for required fields, value types, unsupported assay values, and missing input files.
- Added dry-run mode to print planned commands without running external tools.
- Added per-sample output directory handling under `output_dir/sample`.
- Added structured `run_metadata.json` output for real runs.
- Added `docs/configuration.md` for the v0.4.0 config format and runner behavior.

### Changed
- Updated README, design notes, output documentation, roadmap, and configuration directory notes for the v0.4.0 single-sample interface.
- Clarified that the default Micro-C workflow does not require restriction enzyme information.

### Notes
- Restriction-fragment-aware Hi-C behavior, restriction fragment generation, and `pairtools restrict` are not implemented in v0.4.0.
- Multi-sample execution, restartable chunk processing, workflow-manager files, containers, CI, and full QC reports remain future work.

## v0.3.0 - Unreleased

### Added
- Added planning documentation for the future project design, expected outputs, configuration area, and examples area.
- Documented the current engine decision status: `mdp.sh` remains the retained legacy script, while Bash-first, Snakemake, and Nextflow remain candidate future directions.

### Changed
- Updated README repository contents and project identity language for the next public-facing milestone.
- Updated the roadmap to reflect completed v0.3.0 project structure and identity work without marking future engine implementation tasks complete.

### Notes
- The modern production-facing pipeline engine is planned but not implemented in this milestone.
- Restartable chunk-based execution, sample-sheet execution, containers, CI, and full QC reports remain future work.

## v0.2.0 milestone - not tagged

### Added
- Added MIT license.
- Added `.gitignore` entries for large genomics input, intermediate, and output files.
- Added `docs/roadmap.md` describing the path toward a production-grade pipeline.

### Changed
- Clarified that this repository currently provides a legacy/minimal Slurm-based Micro-C preprocessing script rather than a modern restartable workflow manager pipeline.
- Updated README documentation to match the current repository contents, inputs, outputs, requirements, and limitations.
- Made `mdp.sh` safer for public use by using stricter shell behavior, checking required commands, using configurable thread count from `SLURM_CPUS_PER_TASK`, and removing in-script software installation.
- Improved FASTA index and chromosome sizes handling in `mdp.sh`.
- Improved argument parsing and error handling in `get_qc.py`.

### Notes
- This was an internal cleanup milestone before active `v0.3.0` development.
- No Git tag or GitHub Release was created for this milestone.
- This milestone preserves the legacy single-sample script model.
- Restartable chunk-based processing is planned for a future milestone and should not be treated as part of v0.2.0.
