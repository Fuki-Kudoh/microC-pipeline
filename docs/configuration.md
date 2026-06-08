# Configuration guide

The v0.4.0 workflow adds a small config-driven interface for one Micro-C sample:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml
```

Use `validate-config` before running a full preprocessing job:

```bash
bin/microc-pipeline validate-config --config config/example.single-sample.yaml
```

Use dry-run mode to print the planned command sequence without running external tools or creating large outputs:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml --dry-run
```

The v0.4.0 workflow is Micro-C-oriented and does not require restriction enzyme information. Restriction-fragment-aware Hi-C behavior may be considered in a future milestone, but it is not implemented now.

## Example config

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

Path values are interpreted by the command shell from the directory where you run `bin/microc-pipeline`. Absolute paths are recommended for shared genome resources.

## Required fields

| Field | Meaning |
| --- | --- |
| `sample` | Single sample identifier used for the per-sample output directory and output filenames. |
| `fastq.r1` | Explicit path to read 1 FASTQ. |
| `fastq.r2` | Explicit path to read 2 FASTQ. |
| `genome.name` | Short genome label, for example `mm10` or `hg38`. |
| `genome.fasta` | Path to the reference FASTA. BWA indexes must already be available for alignment. |
| `output_dir` | Root output directory. The workflow writes to `output_dir/sample`. |

## Optional fields and defaults

| Field | Default | Meaning |
| --- | --- | --- |
| `assay` | `microc` | v0.4.0 supports only `microc`. Unsupported values fail clearly. |
| `genome.chrom_sizes` | `null` | Optional chromosome sizes file. If absent or null, the runner derives chromosome sizes from `${genome.fasta}.fai`. |
| `threads` | `SLURM_CPUS_PER_TASK`, otherwise `16` | Thread count passed to supported tools. |
| `bin_sizes` | `[1000]` | Configured contact matrix bin sizes. v0.4.0 uses only the first value for `.cool` generation and warns when more are provided. |
| `outputs.keep_bam` | `true` | Keep or remove the final `bam/SAMPLE.PT.bam` and index after downstream outputs are created. |
| `outputs.make_hic` | `true` | Run `juicer_tools pre` to create `hic/SAMPLE.hic`. |
| `outputs.make_mcool` | `true` | Run `cooler zoomify` to create a multiresolution Cooler output. |

## Assay behavior

If `assay` is omitted, the runner treats the sample as Micro-C:

```yaml
assay: microc
```

No `hic` assay mode is implemented in v0.4.0. For example, `assay: hic` fails with:

```text
Unsupported assay: hic. v0.4.0 supports only assay: microc.
```

## Genome FASTA index and chromosome sizes

The runner uses `genome.fasta` as the reference FASTA. If `genome.chrom_sizes` is provided, that file is copied into the sample output directory. If `genome.chrom_sizes` is absent or null, the runner requires a non-empty `${genome.fasta}.fai`, or creates one with `samtools faidx` when the FASTA directory is writable.

For real runs, the workflow writes the chromosome sizes file to:

```text
results/SAMPLE/genome/{genome.name}.chrom.sizes
```

That file is used by Pairtools and Cooler.

## Output directory behavior

The config-driven workflow writes under a per-sample directory:

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

If `output_dir: .` is set explicitly, the per-sample directory is created in the current directory as `./SAMPLE/`. The retained legacy `mdp.sh` script continues to use its historical repository-root output directories.

## Current limitations

- v0.4.0 supports one sample per command.
- Multi-sample sample sheets are not implemented.
- Restartable chunk-based execution is not implemented.
- Snakemake and Nextflow workflow-manager implementations are not included.
- Containers and CI are not included.
- Full HTML QC reports and project-level QC summaries are not included.
- Restriction-fragment-aware Hi-C behavior, restriction fragment generation, and `pairtools restrict` are not implemented.
- Loop calling, compartment calling, TAD calling, differential contact analysis, and biological interpretation are outside the preprocessing scope.
