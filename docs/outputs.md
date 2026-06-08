# Outputs

This document describes the v0.5.0 standardized output layout for the config-driven single-sample runner and the retained legacy `mdp.sh` layout.

## Scope

The v0.5.0 output standardization applies to:

```bash
bin/microc-pipeline run --config config/example.single-sample.yaml
```

It does not change the retained legacy `mdp.sh` output layout.

The config-driven runner is Micro-C-first. It does not require restriction enzyme information, does not generate restriction fragments, and does not run `pairtools restrict`. Optional `.hic` output means the Juicer contact-map file format only; it is not a separate enzyme-aware Hi-C assay mode.

## Final output layout

For a sample named `SAMPLE` and `output_dir: results`, final products are written under `results/SAMPLE/`:

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
  run_metadata.json
  output_manifest.json
```

Optional final BAM products are retained only when `outputs.keep_bam: true`:

```text
results/SAMPLE/
  bam/SAMPLE.PT.bam
  bam/SAMPLE.PT.bam.bai
```

When `outputs.keep_bam: false`, the final BAM and BAI are removed after downstream commands complete. They are not required by output validation in that mode.

## Final products

| Product | Path | Controlled by | Notes |
| --- | --- | --- | --- |
| Valid pairs | `pairs/SAMPLE.valid.pairs.gz` | Always expected | BGZF-compressed deduplicated valid pairs. |
| Pairix index | `pairs/SAMPLE.valid.pairs.gz.px2` | Always expected | Pairix index for the valid pairs file. |
| Cooler matrix | `cool/SAMPLE.cool` | Always expected | Created with the first configured bin size. |
| Multi-resolution cooler | `cool/SAMPLE.mcool` | `outputs.make_mcool` | Created with explicit `cooler zoomify -o cool/SAMPLE.mcool`. |
| Juicer contact map | `hic/SAMPLE.hic` | `outputs.make_hic` | Juicer `.hic` file format, not enzyme-aware Hi-C mode. |
| Pairtools stats | `stats/SAMPLE.pairtools.stats.txt` | Always expected | Deduplication stats used by `get_qc.py`. |
| Preseq complexity | `stats/SAMPLE.preseq.lc_extrap.txt` | Always expected | `preseq lc_extrap` output. |
| QC TSV | `qc/SAMPLE.qc.tsv` | Always expected | Machine-readable QC summary with `metric`, `value`, and `percent` columns. |
| Final BAM | `bam/SAMPLE.PT.bam` | `outputs.keep_bam` | Deduplicated final BAM. |
| Final BAM index | `bam/SAMPLE.PT.bam.bai` | `outputs.keep_bam` | BAI for final BAM. |

## Temporary and intermediate files

The runner may create temporary files under these directories while commands are executing:

```text
results/SAMPLE/temp/
results/SAMPLE/bam/
results/SAMPLE/pairs/
```

Examples include trimmed FASTQs, SAM, Pairtools `.pairsam` intermediates, undeduplicated BAM files, an uncompressed `pairs/SAMPLE.valid.pairs`, and raw BAM files. These are implementation details and should not be treated as stable final products.

## Output validation

A real `run` validates expected final outputs before printing success. The same checks can be run against an existing output directory without rerunning preprocessing:

```bash
bin/microc-pipeline validate-outputs --config config/example.single-sample.yaml
```

Validation is intentionally lightweight and conservative. It checks file presence, non-empty sizes, and small format indicators needed by this repository. It does not prove biological correctness.

Required checks are:

- `pairs/SAMPLE.valid.pairs.gz` exists and is greater than zero bytes.
- `pairs/SAMPLE.valid.pairs.gz.px2` exists and is greater than zero bytes.
- `cool/SAMPLE.cool` exists and is greater than zero bytes.
- `cool/SAMPLE.mcool` exists and is greater than zero bytes when `outputs.make_mcool: true`.
- `hic/SAMPLE.hic` exists and is greater than zero bytes when `outputs.make_hic: true`.
- `stats/SAMPLE.pairtools.stats.txt` exists, is greater than zero bytes, and includes the Pairtools stats keys needed by `get_qc.py`.
- `stats/SAMPLE.preseq.lc_extrap.txt` exists and is greater than zero bytes.
- `qc/SAMPLE.qc.tsv` exists, is greater than zero bytes, and has header columns `metric`, `value`, and `percent`.
- `bam/SAMPLE.PT.bam` and `bam/SAMPLE.PT.bam.bai` exist and are greater than zero bytes when `outputs.keep_bam: true`.

The current validation layer does not require deep optional checks such as `cooler info`, `cooler ls`, `pairix -l`, or `samtools quickcheck`; those may be added later as best-effort checks.

## `output_manifest.json`

`output_manifest.json` records expected and observed final outputs. It is written after a real run validates outputs and is also updated by `validate-outputs`.

Representative schema:

```json
{
  "sample": "SAMPLE",
  "assay": "microc",
  "output_dir": "results/SAMPLE",
  "pipeline_version": "v0.5.0-dev",
  "outputs": {
    "valid_pairs": {
      "path": "results/SAMPLE/pairs/SAMPLE.valid.pairs.gz",
      "index": "results/SAMPLE/pairs/SAMPLE.valid.pairs.gz.px2",
      "required": true,
      "validated": true,
      "checks": {
        "path": {"exists": true, "size_bytes": 12345},
        "index": {"exists": true, "size_bytes": 678}
      }
    },
    "cool": {
      "path": "results/SAMPLE/cool/SAMPLE.cool",
      "required": true,
      "validated": true
    },
    "mcool": {
      "path": "results/SAMPLE/cool/SAMPLE.mcool",
      "required": true,
      "validated": true
    },
    "hic": {
      "path": "results/SAMPLE/hic/SAMPLE.hic",
      "required": true,
      "validated": true,
      "note": "Juicer .hic contact-map file format; not enzyme-aware Hi-C assay mode."
    },
    "pairtools_stats": {
      "path": "results/SAMPLE/stats/SAMPLE.pairtools.stats.txt",
      "required": true,
      "validated": true
    },
    "preseq": {
      "path": "results/SAMPLE/stats/SAMPLE.preseq.lc_extrap.txt",
      "required": true,
      "validated": true
    },
    "qc_tsv": {
      "path": "results/SAMPLE/qc/SAMPLE.qc.tsv",
      "required": true,
      "validated": true
    },
    "bam": {
      "path": "results/SAMPLE/bam/SAMPLE.PT.bam",
      "index": "results/SAMPLE/bam/SAMPLE.PT.bam.bai",
      "required": true,
      "validated": true
    }
  }
}
```

If a toggle disables an optional output, the manifest keeps an entry with `required: false`, `expected: false`, and a skipped validation result.

## `run_metadata.json`

`run_metadata.json` records how the run was configured and executed. It includes sample, assay, config path, genome resources, FASTQ paths, thread count, bin sizes, output toggles, pipeline version, timestamp, commands, and a `final_outputs` mapping. It does not include restriction enzyme information.

`run_metadata.json` is about execution provenance. `output_manifest.json` is about expected and validated products.

## QC summary formats

`get_qc.py` remains backward compatible. Its default output is the historical text table:

```bash
python3 get_qc.py -p stats/SAMPLE.pairtools.stats.txt --format text
```

The config-driven runner uses the new TSV mode:

```bash
python3 get_qc.py -p stats/SAMPLE.pairtools.stats.txt --format tsv
```

The TSV header is:

```text
metric	value	percent
```

## Retained legacy `mdp.sh` outputs

The retained legacy script remains available:

```bash
sbatch --time=24:00:00 --cpus-per-task=16 --mem=64g \
  mdp.sh SAMPLE mm10 /path/to/mm10.fa
```

It uses historical root-level directories such as `fastqc/`, `BAM/`, `pairs/`, `cool/`, `hic/`, `stats/`, and `qc/`, plus historical naming. The v0.5.0 `output_manifest.json`, standardized names, and `validate-outputs` command are for the config-driven runner.

## Future planned output improvements

Later milestones may add deeper external validators, full QC reports, project-level summaries, restartable execution, multi-sample manifests, and workflow-manager integrations. Those are intentionally outside v0.5.0.
