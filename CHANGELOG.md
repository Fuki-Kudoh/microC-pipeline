# Changelog

All notable changes to this project will be documented in this file.

## v0.2.0 - 2026-06-08

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
- This release preserves the legacy single-sample script model.
- Restartable chunk-based processing is planned for a future milestone and should not be treated as part of v0.2.0.
