# Configuration

The config-driven runner accepts one YAML file for one Micro-C sample:

```bash
bin/wulf3c run --config config/example.single-sample.yaml
```

Validate configuration without running preprocessing:

```bash
bin/wulf3c validate-config --config config/example.single-sample.yaml
```

Validate outputs from an existing completed run. This command uses the config only to infer expected output paths and toggles, so it does not require the original FASTQs, genome FASTA, or BWA index sidecars to still be present:

```bash
bin/wulf3c validate-outputs --config config/example.single-sample.yaml
```

## Example

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

## Fields

| Field | Required | Default | Description |
| --- | --- | --- | --- |
| `sample` | Yes | None | Sample identifier used in the per-sample output directory and standardized output names. Must match `^[A-Za-z0-9][A-Za-z0-9._-]*$`. |
| `assay` | No | `microc` | Only `microc` is supported. No enzyme-aware Hi-C mode is implemented. |
| `fastq.r1` | Yes | None | R1 FASTQ path. |
| `fastq.r2` | Yes | None | R2 FASTQ path. |
| `genome.name` | Yes | None | Genome label passed to tools that need a genome name. |
| `genome.fasta` | Yes | None | Reference FASTA path. Non-empty BWA index sidecars (`.amb`, `.ann`, `.bwt`, `.pac`, `.sa`) must already exist before `validate-config` or `run`. |
| `genome.chrom_sizes` | No | `null` | Optional chromosome sizes file. If omitted, the runner uses `genome.fasta.fai` or creates it with `samtools faidx` when the FASTA directory is writable. |
| `threads` | No | `$SLURM_CPUS_PER_TASK` or `16` | Positive integer thread count. |
| `bin_sizes` | No | `[1000]` | Positive integer bin sizes. v0.5.0 uses the first value for `.cool`; optional `.mcool` is produced by `cooler zoomify`. |
| `output_dir` | Yes | None | Root output directory. Final outputs are under `output_dir/sample`. |
| `outputs.keep_bam` | No | `true` | Retain and validate `bam/SAMPLE.PT.bam` and `bam/SAMPLE.PT.bam.bai`. If false, final BAM files are removed and are not expected by validation. |
| `outputs.make_hic` | No | `true` | Create and validate `hic/SAMPLE.hic`. This is the Juicer contact-map file format, not enzyme-aware Hi-C assay support. |
| `outputs.make_mcool` | No | `true` | Create and validate `cool/SAMPLE.mcool`. |

## Output toggles and validation

The output toggles define both the commands that are planned and the final files expected by validation:

- `outputs.keep_bam: true` makes `bam/SAMPLE.PT.bam` and `bam/SAMPLE.PT.bam.bai` required final outputs.
- `outputs.keep_bam: false` removes final BAM products after downstream outputs are produced; BAM files are skipped by output validation.
- `outputs.make_hic: true` makes `hic/SAMPLE.hic` a required validated output.
- `outputs.make_hic: false` skips `.hic` creation and validation.
- `outputs.make_mcool: true` makes `cool/SAMPLE.mcool` a required validated output.
- `outputs.make_mcool: false` skips `.mcool` creation and validation.

Always-required final outputs are the BGZF-compressed valid pairs file, Pairix index, `.cool`, Pairtools stats, Preseq output, QC TSV, `run_metadata.json`, and `output_manifest.json`. For real runs, `run_metadata.json` includes best-effort `tool_versions` records for required command-line tools.

## Dry run

Dry-run mode validates the config and BWA index preflight, then prints planned commands, expected final outputs, and planned validation checks without running external tools:

```bash
bin/wulf3c run --config config/example.single-sample.yaml --dry-run
```

Dry-run mode is intended for inspection. It may create no files at all.

## Boundaries

The config format is intentionally single-sample and Micro-C-first. It does not include a sample sheet, restart/resume settings, workflow-manager settings, container settings, restriction enzyme fields, restriction fragment files, or downstream biological interpretation options.
