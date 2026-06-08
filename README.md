# microC-pipeline

`microC-pipeline` is a Micro-C preprocessing repository. The v0.5.0 development milestone keeps the lightweight config-driven single-sample runner introduced in v0.4.0 and makes its valid pairs, matrix, stats, QC, and manifest files first-class validated products.

The preferred interface is:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml
```

The retained legacy interface remains available for backward compatibility:

```bash
sbatch --time=24:00:00 --cpus-per-task=16 --mem=64g \
  mdp.sh SAMPLE mm10 /path/to/mm10.fa
```

This repository is still **not** a modern restartable workflow-manager pipeline. It does not implement restartable chunk-based execution, multi-sample sample sheets, Snakemake or Nextflow workflows, containers, CI, full HTML QC reports, enzyme-aware Hi-C processing, restriction fragment generation, loop calling, compartment calling, TAD calling, differential contact analysis, or downstream biological interpretation.

## What this repository does

For one Micro-C sample, the config-driven runner follows the same core preprocessing path as the retained legacy script:

1. FASTQ quality checks with FastQC.
2. Adapter/quality trimming with Trim Galore.
3. Alignment with BWA-MEM using `bwa mem -5SP -T0`.
4. Pair parsing and sorting with Pairtools.
5. Undeduplicated BAM split, sorting, and indexing.
6. Library complexity estimation with Preseq.
7. Pairtools deduplication and stats.
8. Machine-readable QC TSV generation with `get_qc.py --format tsv`.
9. Final deduplicated valid pairs and BAM generation.
10. Optional Juicer `.hic` contact-map file generation.
11. BGZF compression, Pairix indexing, `.cool`, and optional `.mcool` generation.
12. Lightweight validation and `output_manifest.json` creation.

The default workflow is Micro-C-oriented. No restriction enzyme is required, the config does not include an enzyme field, and restriction-fragment-aware Hi-C processing is not implemented. The optional `.hic` output is the Juicer contact-map **file format**; it does not imply enzyme-aware Hi-C assay mode.

## Repository contents

```text
README.md                         Project overview and usage notes
bin/microc-pipeline               v0.5.0 config-driven single-sample CLI
microc_pipeline/                  Python package for config, output paths, validation, and command execution
config/example.single-sample.yaml Example single-sample Micro-C config
config/README.md                  Configuration directory notes
docs/configuration.md             Config reference
docs/design.md                    Design notes and current boundaries
docs/outputs.md                   Standardized output layout and manifest schema
docs/roadmap.md                   Development roadmap
mdp.sh                            Retained legacy/minimal Slurm Micro-C preprocessing script
get_qc.py                         Pairtools-stats QC summarizer with text and TSV output
CHANGELOG.md                      Repository-level change history
LICENSE                           MIT license
examples/README.md                Placeholder for future tiny synthetic examples
.gitignore                        Ignore rules for large genomics outputs and local files
```

## Config-driven usage

Start from `config/example.single-sample.yaml`:

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

Validate a config:

```bash
bin/microc-pipeline validate-config --config config/example.single-sample.yaml
```

Print planned commands, expected final outputs, and planned validation checks without running external tools:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml --dry-run
```

Run one sample:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml
```

Validate outputs from an existing completed run without rerunning preprocessing:

```bash
bin/microc-pipeline validate-outputs --config config/example.single-sample.yaml
```

See `docs/configuration.md` for required fields, defaults, dry-run behavior, genome FASTA index handling, output toggles, and current limitations.

## Standardized final output layout

The config-driven runner writes under a per-sample directory and validates the final products before reporting success:

```text
results/SAMPLE/
  pairs/SAMPLE.valid.pairs.gz
  pairs/SAMPLE.valid.pairs.gz.px2
  cool/SAMPLE.cool
  cool/SAMPLE.mcool              # when outputs.make_mcool: true
  hic/SAMPLE.hic                 # when outputs.make_hic: true
  stats/SAMPLE.pairtools.stats.txt
  stats/SAMPLE.preseq.lc_extrap.txt
  qc/SAMPLE.qc.tsv
  bam/SAMPLE.PT.bam              # when outputs.keep_bam: true
  bam/SAMPLE.PT.bam.bai          # when outputs.keep_bam: true
  run_metadata.json
  output_manifest.json
```

`run_metadata.json` records how the run was configured and executed, including standardized final output paths. `output_manifest.json` records which final outputs were expected and whether validation passed. Detailed output semantics, validation checks, and manifest fields are documented in `docs/outputs.md`.

## Requirements

Install and configure dependencies before launching a real run. The config-driven runner checks for required commands and exits early if any are missing:

- `fastqc`
- `trim_galore`
- `samtools`
- `bwa`
- `pairtools`
- `preseq`
- `bgzip`
- `pairix`
- `cooler`
- `python3`

The runner checks `juicer_tools` only when `outputs.make_hic: true`.

The required genome FASTA should already have BWA indexes available for alignment. If `genome.chrom_sizes` is not provided, `genome.fasta` should have an existing non-empty `.fai` index, or the FASTA directory must be writable so the runner can create the index with `samtools faidx`.

### HPC modules are site-specific

Many HPC clusters expose tools through environment modules, but module names vary by site. Neither the config-driven runner nor `mdp.sh` installs software during a run; missing commands are reported as dependency errors.

## Retained legacy `mdp.sh` usage

`mdp.sh` remains available as the retained legacy/minimal direct Slurm script. It uses historical root-level output directories and historical FASTQ naming assumptions.

Run the script from the repository root. It expects paired FASTQs in `fastq/` using this naming convention:

```text
fastq/{SAMPLE_ID}_R1_001.fastq.gz
fastq/{SAMPLE_ID}_R2_001.fastq.gz
```

Submit a legacy Slurm job:

```bash
sbatch --time=24:00:00 --cpus-per-task=16 --mem=64g \
  mdp.sh <SAMPLE_ID> <GENOME_NAME> <GENOME_FASTA>
```

The v0.5.0 standardized output naming and validation layer apply to the config-driven runner, not to the retained legacy root-level layout.
