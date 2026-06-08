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

6. **Micro-C / Hi-C compatibility**
   - The workflow should support Micro-C and Hi-C-like paired-end proximity ligation data.
   - Enzyme/restriction-fragment assumptions should be explicit and configurable.

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

Planned tasks:

- Add `CHANGELOG.md`.
- Add `docs/` with design notes and output specifications.
- Move the current `mdp.sh` into `legacy/` or explicitly label it as legacy in place.
- Add a modern entrypoint plan, for example `bin/microc-pipeline`.
- Add `config/` examples.
- Add `examples/` with tiny synthetic or documented toy inputs.
- Define repository topics and description on GitHub.
- Decide whether the main engine will be Bash-first, Snakemake, or Nextflow.

Expected status after this milestone:

```text
The project has a clear identity: legacy script retained, modern pipeline planned.
```

## Milestone v0.4.0: config-driven single-sample workflow

Goal: remove hard-coded sample naming and path assumptions.

Planned tasks:

- Add a config file format, for example `config.yaml`.
- Add support for explicit FASTQ paths.
- Add configurable genome name, FASTA, chrom sizes, enzyme, bin sizes, and output directory.
- Add a dry-run or command-printing mode.
- Add structured run metadata output.
- Keep a single-sample mode as the simplest supported workflow.

Example future config:

```yaml
sample: SAMPLE_ID
fastq1: fastq/SAMPLE_R1.fastq.gz
fastq2: fastq/SAMPLE_R2.fastq.gz
genome_name: mm10
genome_fasta: /path/to/mm10.fa
enzyme: MboI
threads: 16
bin_sizes:
  - 1000
  - 5000
  - 10000
  - 25000
  - 100000
outputs:
  keep_bam: false
  make_hic: true
  make_mcool: true
```

Expected status after this milestone:

```text
Users can run a single sample without editing the script itself.
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
  qc/SAMPLE.qc.tsv
```

Optional outputs:

```text
results/SAMPLE/
  bam/SAMPLE.PT.bam
  bam/SAMPLE.PT.bam.bai
```

Planned tasks:

- Standardize valid pairs output naming.
- Ensure valid pairs are BGZF-compressed and pairix-indexed.
- Add output validation for pairs, pairix index, cool, mcool, and hic where feasible.
- Add optional BAM retention via config.
- Document all outputs in `docs/outputs.md`.

Expected status after this milestone:

```text
The pipeline has clear, validated, documented final outputs.
```

## Milestone v0.6.0: restartable chunk engine

Goal: integrate the restartability concepts developed in HiC-Nap.

Planned tasks:

- Add FASTQ chunking.
- Add per-chunk processing:
  - trimming
  - alignment
  - pairtools parse
  - pairtools sort
  - pairtools restrict/select when appropriate
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
- Allow per-sample FASTQ, genome, enzyme, and condition metadata.
- Process failed samples independently on rerun.
- Generate project-level QC summary.
- Add optional grouping by condition or batch.

Example sample sheet:

```csv
sample,fastq1,fastq2,genome,enzyme,condition
WT_rep1,fastq/WT1_R1.fastq.gz,fastq/WT1_R2.fastq.gz,mm10,MboI,WT
WT_rep2,fastq/WT2_R1.fastq.gz,fastq/WT2_R2.fastq.gz,mm10,MboI,WT
DS_rep1,fastq/DS1_R1.fastq.gz,fastq/DS1_R2.fastq.gz,mm10,MboI,DS
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
