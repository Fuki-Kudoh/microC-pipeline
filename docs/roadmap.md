# microC-pipeline roadmap

This roadmap describes the intended path from the current legacy/minimal Slurm script toward a production-grade Micro-C / Hi-C preprocessing pipeline.

The long-term goal is to make `microC-pipeline` a restartable, reproducible, inspectable workflow that takes paired-end FASTQ files to validated contact files, contact matrices, and QC reports.

> Working vision: **A restartable, reproducible Micro-C/Hi-C preprocessing pipeline from FASTQ to contact maps and QC reports.**

## Current state: v0.1.x

The current repository contains a compact legacy Slurm workflow centered on:

- `mdp.sh`: a single-sample Slurm-oriented Micro-C preprocessing script
- `get_qc.py`: a lightweight Pairtools-stats summarizer
- `README.md`: public-facing usage and limitations

The current workflow can run a conventional path:

```text
FASTQ
→ FastQC
→ Trim Galore
→ BWA-MEM
→ pairtools parse/sort/dedup/split
→ BAM / pairs.gz / hic / cool / mcool
→ lightweight QC summary
```

This version is intentionally simple. It is useful as a readable legacy script, but it is not yet a modern restartable workflow manager pipeline.

## Design principles for v1.0.0

The v1.0.0 release should be built around the following principles.

1. **Reproducibility**
   - Tool versions should be documented or pinned.
   - Outputs should be deterministic when inputs and settings are unchanged.
   - Commands and configuration should be recorded with each run.

2. **Restartability**
   - Failed or interrupted runs should resume safely.
   - Completed steps should be skipped only when both status and outputs validate.
   - Large long-running jobs should avoid starting from zero after interruption.

3. **Config-driven execution**
   - Hard-coded paths and sample naming should be replaced by config files and sample sheets.
   - The same pipeline should support local workstations and Slurm-based HPC systems.

4. **Transparent outputs**
   - Every major output should have a documented path and meaning.
   - Intermediate files should be either retained by option or safely cleaned.
   - The pipeline should make it clear what is final, temporary, or diagnostic.

5. **Strong QC**
   - The pipeline should not only produce contact files, but also help decide whether a library is usable.
   - QC should summarize mapping, duplication, valid pairs, cis/trans balance, distance decay, and matrix-level properties.

6. **Micro-C first, with future Hi-C compatibility considered carefully**
   - The current workflow should remain Micro-C-oriented by default.
   - Any future Hi-C-specific assumptions should be documented explicitly and added only in a later milestone.

## Relationship to HiC-Nap

HiC-Nap is currently a separate restartable chunk-based preprocessing prototype. Its design is directly relevant to the future of `microC-pipeline`.

Planned direction:

```text
HiC-Nap = restartable chunk-processing prototype / engine
microC-pipeline = final integrated production-facing project
```

Key HiC-Nap concepts to bring into `microC-pipeline`:

- FASTQ chunking for large datasets
- Per-chunk status tracking
- Conservative restart behavior
- Chunk-level trimming, alignment, and Pairtools processing
- Global merge and global deduplication
- BGZF-compressed, pairix-indexed valid pairs output
- Line-based progress logging suitable for SSH, tmux, and Slurm logs

## Milestone v0.2.0: legacy repository polish

Goal: make the current public repository internally consistent and safe to inspect.

Planned or completed tasks:

- [x] Clarify that the current repository is a legacy/minimal Slurm workflow.
- [x] Remove README references to missing files or unrealized features.
- [x] Add `.gitignore` for large genomics inputs and outputs.
- [x] Harden `mdp.sh` with stricter shell behavior and dependency checks.
- [x] Make thread count configurable from `SLURM_CPUS_PER_TASK`.
- [x] Remove in-script software installation.
- [x] Improve FASTA index and chromosome sizes handling.
- [x] Improve `get_qc.py` argument parsing and error handling.
- [x] Add a license file.
- [x] Update README after license addition.
- [x] Add a changelog.

Expected status after this milestone:

```text
microC-pipeline is public-safe and useful as a documented legacy Slurm script.
```

## Milestone v0.3.0: project structure and identity

Goal: prepare the repository for active development toward v1.0.0.

Completed project structure and identity tasks:

- [x] Add a `v0.3.0 - Unreleased` changelog section.
- [x] Add `docs/design.md` with design notes and engine-decision status.
- [x] Add `docs/outputs.md` with current legacy output notes and planned output-documentation direction.
- [x] Explicitly label `mdp.sh` as the retained legacy/minimal Slurm script in public-facing documentation.
- [x] Add `config/README.md` as a placeholder for future configuration examples and schemas.
- [x] Add `examples/README.md` as a placeholder for future tiny synthetic examples.
- [x] Document that Bash-first, Snakemake, and Nextflow remain candidate future engine directions.
- [x] Document that HiC-Nap restartable chunk-processing concepts are planning inputs for future development.

Remaining or future tasks, not completed in this milestone:

- [ ] Add a modern executable entrypoint, for example `bin/microc-pipeline`.
- [ ] Define GitHub repository topics and description.
- [ ] Make a final production engine choice if project ownership gives that direction.
- [ ] Implement the modern pipeline engine.
- [ ] Implement restartable chunk-based execution.
- [ ] Implement sample-sheet execution, workflow-manager files, containers, CI, or full QC reports.

Expected status after this milestone:

```text
The project has a clear identity: legacy script retained, modern pipeline planned.
```

## Milestone v0.4.0: config-driven single-sample workflow

Goal: remove hard-coded sample naming and path assumptions for one Micro-C sample.

Completed tasks:

- [x] Add the `bin/microc-pipeline` executable entrypoint.
- [x] Add a single-sample YAML config format.
- [x] Add `config/example.single-sample.yaml` without restriction enzyme settings.
- [x] Add support for explicit `fastq.r1` and `fastq.r2` paths.
- [x] Add configurable genome name, FASTA, optional chromosome sizes, bin sizes, threads, output directory, and output toggles.
- [x] Default `assay` to `microc` and reject unsupported assay values.
- [x] Add config validation with clear missing-field, type, assay, and missing-file errors.
- [x] Add dry-run command printing.
- [x] Add per-sample output directory handling under `output_dir/sample`.
- [x] Add structured run metadata for real runs.
- [x] Keep the retained legacy `mdp.sh` script available and documented.

Explicit non-goals for v0.4.0:

- [ ] Multi-sample sample-sheet execution.
- [ ] Restartable chunk-based execution.
- [ ] Snakemake or Nextflow workflow-manager implementation.
- [ ] Containers or CI.
- [ ] Full HTML QC reports or project-level QC summaries.
- [ ] Restriction-fragment-aware Hi-C behavior, restriction fragment generation, or `pairtools restrict`.

Example current config:

```yaml
sample: SAMPLE_ID
assay: microc
fastq:
  r1: fastq/SAMPLE_R1.fastq.gz
  r2: fastq/SAMPLE_R2.fastq.gz
genome:
  name: mm10
  fasta: /path/to/mm10.fa
  chrom_sizes: null
threads: 16
bin_sizes:
  - 1000
output_dir: results
outputs:
  keep_bam: true
  make_hic: true
  make_mcool: true
```

The v0.4.0 workflow is Micro-C-oriented and does not require restriction enzyme information. Restriction-fragment-aware Hi-C behavior may be considered in a future milestone, but it is not implemented now.

Expected status after this milestone:

```text
Users can run a single Micro-C sample without editing the legacy script itself.
```

## Milestone v0.5.0: valid pairs and matrix outputs as first-class products

Goal: define and validate the core final outputs.

Planned final outputs:

```text
results/SAMPLE/
  pairs/SAMPLE.valid.pairs.gz
  pairs/SAMPLE.valid.pairs.gz.px2
  cool/SAMPLE.cool
  cool/SAMPLE.mcool
  hic/SAMPLE.hic
  stats/SAMPLE.pairtools.stats.txt
  stats/SAMPLE.preseq.lc_extrap.txt
  qc/SAMPLE.qc.tsv
```

Optional outputs:

```text
results/SAMPLE/
  bam/SAMPLE.PT.bam
  bam/SAMPLE.PT.bam.bai
```

Completed v0.5.0 tasks:

- [x] Standardized valid pairs output naming as `pairs/SAMPLE.valid.pairs.gz`.
- [x] Ensured valid pairs are BGZF-compressed and pairix-indexed as `pairs/SAMPLE.valid.pairs.gz.px2`.
- [x] Added output validation for pairs, pairix index, cool, mcool, hic, stats, QC TSV, and optional BAM where feasible.
- [x] Added optional BAM retention via `outputs.keep_bam`.
- [x] Documented all standardized outputs, validation checks, and the manifest schema in `docs/outputs.md`.

Expected status after this milestone:

```text
The pipeline has clear, validated, documented final outputs.
```

### Milestone v0.5.1: preflight and reproducibility hardening

Completed v0.5.1 tasks:

- [x] Added safe sample name validation.
- [x] Added BWA index preflight validation.
- [x] Added best-effort tool version collection in `run_metadata.json`.
- [x] Removed minor unreachable/overly-local implementation details without changing core processing commands.

This patch does not add restartable chunking, multi-sample support, workflow-manager support, or full QC reports.


## Milestone v0.6.0: restartable chunk engine

Goal: integrate the restartability concepts developed in HiC-Nap.

Planned tasks:

- Add FASTQ chunking.
- Add per-chunk processing:
  - trimming
  - alignment
  - pairtools parse
  - pairtools sort
  - pairtools parse/sort-compatible selection steps when appropriate
- Add per-chunk status files.
- Add conservative restart behavior.
- Add global merge and global deduplication.
- Produce sample-level valid pairs.
- Support partial/nightly runs with a maximum number of chunks.
- Add progress output suitable for long SSH/tmux/Slurm sessions.

Expected status after this milestone:

```text
Large datasets can be processed incrementally and resumed safely.
```

## Milestone v0.7.0: QC report and decision support

Goal: turn QC into a strong feature rather than an afterthought.

Planned QC metrics:

- total read pairs
- mapped read pairs and mapping rate
- unmapped read pairs
- duplicate rate
- no-duplicate read pairs
- valid pairs
- cis/trans ratio
- short-range cis fraction
- long-range cis fraction
- library complexity estimate
- contact distance decay, P(s)
- chromosome-level contact balance
- matrix summary statistics

Planned outputs:

```text
results/SAMPLE/qc/
  SAMPLE.qc.tsv
  SAMPLE.qc.json
  SAMPLE.qc.html
  figures/
    mapping_summary.png
    cis_trans_summary.png
    distance_decay.png
    contact_map_preview.png
```

Expected status after this milestone:

```text
Users can judge library quality from the pipeline output without manual log digging.
```

## Milestone v0.8.0: multi-sample support

Goal: support project-scale processing.

Planned tasks:

- Add sample sheet support.
- Allow per-sample FASTQ, genome, assay, and condition metadata.
- Process failed samples independently on rerun.
- Generate project-level QC summary.
- Add optional grouping by condition or batch.

Example sample sheet:

```csv
sample,fastq1,fastq2,genome,assay,condition
WT_rep1,fastq/WT1_R1.fastq.gz,fastq/WT1_R2.fastq.gz,mm10,microc,WT
WT_rep2,fastq/WT2_R1.fastq.gz,fastq/WT2_R2.fastq.gz,mm10,microc,WT
DS_rep1,fastq/DS1_R1.fastq.gz,fastq/DS1_R2.fastq.gz,mm10,microc,DS
```

Expected status after this milestone:

```text
A complete project can be run from a sample sheet.
```

## Milestone v0.9.0: reproducibility, packaging, and tests

Goal: prepare for a stable v1.0.0 release.

Planned tasks:

- Add conda/mamba environment file.
- Add container support if feasible.
- Add small smoke-test dataset or synthetic fixture.
- Add GitHub Actions for syntax and smoke tests.
- Add test commands for shell and Python components.
- Add `CHANGELOG.md` entries for all releases.
- Add `CITATION.cff` if a formal citation or DOI is available.
- Add documentation for local and Slurm execution.

Expected status after this milestone:

```text
The pipeline is testable, installable, and release-ready.
```

## v1.0.0 release criteria

The v1.0.0 release should satisfy all of the following criteria.

- [ ] FASTQ to valid pairs, cool/mcool, hic, and QC can be run from documented inputs.
- [ ] Single-sample and multi-sample modes are supported.
- [ ] Restart/resume behavior is safe and documented.
- [ ] Local and Slurm execution modes are supported or clearly documented.
- [ ] Configuration is file-driven; users do not need to edit pipeline code.
- [ ] Core outputs are validated before steps are marked complete.
- [ ] QC report is sufficient for first-pass library assessment.
- [ ] A small test or smoke-test dataset is available.
- [ ] CI runs syntax checks and at least a minimal smoke test.
- [ ] Documentation covers installation, quickstart, outputs, troubleshooting, and design assumptions.
- [ ] License is present.
- [ ] Citation guidance is present, even if only repository citation is recommended.
- [ ] The pipeline has been validated on at least one real Micro-C or Hi-C dataset.

## Non-goals before v1.0.0

The following are useful but not required for v1.0.0:

- cloud execution
- full nf-core compliance
- graphical user interface
- automatic biological interpretation
- downstream differential contact analysis
- loop calling
- compartment calling
- TAD/domain calling

These may be added after the core preprocessing pipeline is stable.

## Strategic direction

The project should grow in two layers:

1. **Stable core engine**
   - reliable FASTQ to contact outputs
   - restartable execution
   - validated files
   - clear logs

2. **Research-facing usability**
   - QC reports
   - project summaries
   - sensible examples
   - documentation that makes the pipeline useful to other researchers

The immediate next step after v0.1.1 is to preserve the cleaned legacy script while designing the modern engine. The most valuable future work is to bring the restartable chunk-processing model from HiC-Nap into `microC-pipeline` and make it the foundation of v1.0.0.
