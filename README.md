# microC-pipeline

`microC-pipeline` is a legacy/minimal Slurm-based Micro-C preprocessing workflow. The repository is intentionally small and centers on one shell script, `mdp.sh`, plus a lightweight Pairtools-stats summarizer, `get_qc.py`.

Active, restartable development has moved to **HiC-Nap**. This repository is retained as a compact legacy workflow for users who want a simple Slurm script that can be inspected and adapted for an HPC site.

## What this repository does

For a single sample, `mdp.sh` runs a classic Micro-C preprocessing path:

1. FASTQ quality checks with FastQC.
2. Adapter/quality trimming with Trim Galore.
3. Alignment with BWA-MEM.
4. Pair parsing, sorting, deduplication, and stats with Pairtools.
5. BAM sorting/indexing with Samtools.
6. Library complexity estimation with Preseq.
7. `.hic` generation with Juicer tools.
8. `.pairs.gz`, `.cool`, and `.mcool` generation with Pairix and Cooler.
9. A small text QC summary from Pairtools stats using `get_qc.py`.

This is **not** a modern, restartable workflow manager pipeline. It does not include sample-sheet orchestration, containers, continuous integration, or a packaged software environment.

## Repository contents

```text
README.md   Project overview and usage notes
mdp.sh      Legacy/minimal Slurm Micro-C preprocessing script
get_qc.py   Small Pairtools-stats QC text summarizer
```

## Expected input and output layout

Run the script from the repository root. It expects paired FASTQs in `fastq/` using the historical naming convention:

```text
fastq/{SAMPLE_ID}_R1_001.fastq.gz
fastq/{SAMPLE_ID}_R2_001.fastq.gz
```

The script creates and uses these output/work directories:

```text
BAM/      sorted/indexed BAM outputs
cool/     .cool and zoomified .mcool outputs
fastqc/   FastQC reports
fastq/    input FASTQs and Trim Galore outputs
hic/      .hic outputs
genome/   generated chromosome sizes files
pairs/    .pairs.gz and pairix index outputs
stats/    pairtools, preseq, and QC summary text outputs
temp/     intermediate SAM/pairsam files
logs/     per-sample pipeline logs
```

Large genomics outputs and inputs are ignored by Git so the public repository stays small and safe.

## Requirements

Install and configure dependencies before submitting a job. The script checks for these commands and exits early if any are missing:

- `fastqc`
- `trim_galore`
- `samtools`
- `bwa`
- `pairtools`
- `preseq`
- `juicer_tools`
- `bgzip`
- `pairix`
- `cooler`
- `python3`

The required genome FASTA should already have any aligner indexes needed by BWA. `GENOME_FASTA` should also have an existing non-empty `.fai` index, or the directory containing `GENOME_FASTA` must be writable so `mdp.sh` can create the index with `samtools faidx`. The script derives `genome/{GENOME_NAME}.chrom.sizes` from that FASTA index for Pairtools and Cooler.

### HPC modules are site-specific

Many HPC clusters expose tools through environment modules, but module names vary by site. `mdp.sh` includes a clearly marked optional module-loading section that you may edit for your cluster. The script does not install software during a run; if a command such as `cooler` is unavailable, it fails during dependency checks.

## Usage

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

Threading defaults to `SLURM_CPUS_PER_TASK` when set, otherwise `16`:

```bash
THREADS="${SLURM_CPUS_PER_TASK:-16}"
```

## QC summary

After Pairtools deduplication, `mdp.sh` writes Pairtools stats to:

```text
stats/{SAMPLE_ID}.txt
```

It then creates a compact text summary with:

```bash
python3 get_qc.py -p stats/SAMPLE.txt > stats/qc_SAMPLE.txt
```

You can rerun that command manually after a completed or partially completed run if `stats/SAMPLE.txt` exists.

## Notes and limitations

- This is a legacy/minimal Slurm script, not a fully restartable production workflow.
- Paths and FASTQ naming are intentionally simple and historical.
- Intermediate files may be removed after downstream files are created, matching the previous script behavior.
- Module names and dependency installation are intentionally left to each HPC site.
- Full biological validation requires running on real data and inspecting the outputs; syntax checks alone do not validate biological correctness.
- No license file is included in this repository at this time.

## Version note

Legacy script version history is retained in comments at the bottom of `mdp.sh`. This README does not assign a separate release version.
