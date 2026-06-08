# microC-pipeline

`microC-pipeline` is a Micro-C preprocessing repository. The v0.4.0 development milestone adds a lightweight config-driven single-sample entrypoint while retaining the historical `mdp.sh` Slurm script for backward compatibility.

The preferred v0.4.0 interface is:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml
```

The retained legacy interface remains:

```bash
sbatch --time=24:00:00 --cpus-per-task=16 --mem=64g \
  mdp.sh SAMPLE mm10 /path/to/mm10.fa
```

This repository is still **not** a modern restartable workflow-manager pipeline. It does not implement restartable chunk-based execution, multi-sample sample sheets, Snakemake or Nextflow workflows, containers, CI, full HTML QC reports, or downstream biological interpretation.

## What this repository does

For one Micro-C sample, the config-driven runner follows the same core preprocessing path as the retained legacy script:

1. FASTQ quality checks with FastQC.
2. Adapter/quality trimming with Trim Galore.
3. Alignment with BWA-MEM using `bwa mem -5SP -T0`.
4. Pair parsing and sorting with Pairtools.
5. Undeduplicated BAM split, sorting, and indexing.
6. Library complexity estimation with Preseq.
7. Pairtools deduplication and stats.
8. A small QC text summary from `get_qc.py`.
9. Final deduplicated pairs and BAM generation.
10. Optional `.hic` generation with Juicer tools.
11. `.pairs.gz`, Pairix index, `.cool`, and optional `.mcool` generation.

The default workflow is Micro-C-oriented. No restriction enzyme is required, the v0.4.0 config does not include an enzyme field, and restriction-fragment-aware Hi-C processing is not implemented in this milestone.

## Repository contents

```text
README.md                         Project overview and usage notes
bin/microc-pipeline               v0.4.0 config-driven single-sample CLI
microc_pipeline/                  Python package for config validation and command execution
config/example.single-sample.yaml Example single-sample Micro-C config
config/README.md                  Configuration directory notes
docs/configuration.md             v0.4.0 config reference
docs/design.md                    Design notes and current boundaries
docs/outputs.md                   Current config-driven and legacy output notes
docs/roadmap.md                   Development roadmap
mdp.sh                            Retained legacy/minimal Slurm Micro-C preprocessing script
get_qc.py                         Small Pairtools-stats QC text summarizer
CHANGELOG.md                      Repository-level change history
LICENSE                           MIT license
examples/README.md                Placeholder for future tiny synthetic examples
.gitignore                        Ignore rules for large genomics outputs and local files
```

## Config-driven v0.4.0 usage

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

Print planned commands without running external tools:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml --dry-run
```

Run one sample:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml
```

The runner requires explicit FASTQ paths in the config rather than the legacy hard-coded `fastq/{SAMPLE_ID}_R1_001.fastq.gz` and `fastq/{SAMPLE_ID}_R2_001.fastq.gz` names. See `docs/configuration.md` for required fields, defaults, dry-run behavior, genome FASTA index handling, and current limitations.

## Config-driven output layout

The v0.4.0 runner writes under a per-sample directory:

```text
results/SAMPLE/
  bam/
  cool/
  fastqc/
  hic/
  genome/
  logs/
  pairs/
  qc/
  stats/
  temp/
  run_metadata.json
```

For real runs, `run_metadata.json` records the resolved sample, assay, config path, output directory, genome resources, FASTQ paths, threads, bin sizes, output toggles, pipeline version, dry-run status, timestamp, and command list. The metadata does not include restriction enzyme information.

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

Example:

```bash
sbatch --time=24:00:00 --cpus-per-task=16 --mem=64g \
  mdp.sh SAMPLE mm10 /path/to/mm10.fa
```

Arguments:

- `SAMPLE_ID`: sample prefix matching FASTQs under `fastq/`.
- `GENOME_NAME`: short genome label used for generated files and Juicer tools, for example `mm10` or `hg38`.
- `GENOME_FASTA`: path to the reference FASTA.

Threading defaults to `SLURM_CPUS_PER_TASK` when set, otherwise `16`.

The legacy script creates and uses root-level directories such as `BAM/`, `cool/`, `fastqc/`, `hic/`, `genome/`, `pairs/`, `stats/`, `temp/`, and `logs/`.

## QC summary

After Pairtools deduplication, both workflows can create a compact text summary from Pairtools stats with `get_qc.py`:

```bash
python3 get_qc.py -p stats/SAMPLE.txt > stats/qc_SAMPLE.txt
```

The config-driven runner writes its QC summary under `results/SAMPLE/qc/`.

## Notes and limitations

- v0.4.0 supports only one sample per config/run command.
- Multi-sample support is not implemented yet.
- Restartable chunk-based processing is not implemented yet.
- The repository does not include workflow-manager files, containers, CI, full QC reports, or project-level QC summaries.
- The default workflow is Micro-C-oriented and does not require restriction enzyme information.
- Enzyme-aware Hi-C processing, restriction fragment generation, and `pairtools restrict` are not implemented.
- Downstream loop calling, compartment calling, TAD calling, differential contact analysis, and biological interpretation are outside the current scope.

Large genomics outputs and inputs are ignored by Git so the public repository stays small and safe.
