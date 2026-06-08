# Output notes

This document distinguishes the current v0.4.0 config-driven output layout from retained legacy `mdp.sh` outputs and future planned outputs.

## Current config-driven outputs

The preferred v0.4.0 single-sample runner writes outputs under `output_dir/sample`. With the example `output_dir: results` and `sample: SAMPLE`, the layout is:

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

Directory roles:

- `bam/`: undeduplicated and final sorted/indexed BAM files. If `outputs.keep_bam: false`, the final `SAMPLE.PT.bam` and `.bai` are removed after downstream files are created.
- `cool/`: first-bin `.cool` output and optional zoomified `.mcool` output.
- `fastqc/`: FastQC outputs.
- `hic/`: optional `SAMPLE.hic` when `outputs.make_hic: true`.
- `genome/`: copied or derived `{genome.name}.chrom.sizes` file used by Pairtools and Cooler.
- `logs/`: per-sample log file with resolved settings and command lines.
- `pairs/`: deduplicated pairs, BGZF-compressed pairs, and Pairix index.
- `qc/`: compact QC text summary produced by `get_qc.py`.
- `stats/`: Pairtools stats and Preseq complexity output.
- `temp/`: intermediate SAM, pairsam, and trimmed FASTQ files.
- `run_metadata.json`: structured run metadata for real runs.

The v0.4.0 runner always creates `.cool` for the first configured bin size. If multiple `bin_sizes` are configured, it warns that only the first value is used. It creates `.hic` only when `outputs.make_hic: true` and runs `cooler zoomify` only when `outputs.make_mcool: true`.

`run_metadata.json` includes sample, assay, config path, output directory, genome FASTA, chromosome sizes path, FASTQ paths, thread count, bin sizes, output toggles, pipeline version, dry-run status, timestamp, and command records. It does not include restriction enzyme information.

## Retained legacy `mdp.sh` outputs

The retained legacy/minimal Slurm script writes outputs and intermediates into repository-root directories such as:

```text
BAM/
cool/
fastqc/
hic/
genome/
pairs/
stats/
temp/
logs/
```

These paths reflect the historical direct script workflow. They are kept for backward compatibility and are separate from the v0.4.0 per-sample `results/SAMPLE/` layout.

## Future planned outputs

Later milestones may standardize final valid pairs and matrix products as first-class outputs, validate output files, and add richer QC files. A candidate future layout remains:

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

Optional future outputs may include retained final BAM files under `results/SAMPLE/bam/`. Restartable chunk outputs, project-level QC summaries, and full reports remain future work.
