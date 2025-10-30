````markdown
# microC-pipeline

Minimal Micro-C preprocessing and QC summarization on Slurm.  
**v0.1.0 – Minimal runnable release.**

> Runs alignment → pairing → cooler generation, plus a tiny QC report.  
> Opinionated defaults; designed for HPC queues (Slurm). Local mode planned.

---

## ✨ Features (v0.1.x)
- Single-sample run via one script: `mdp.sh`
- Produce BAM / pairs / `.cool` with sensible paths
- Tiny QC export (`QC.tsv` + `QC.html`) via `get_qc.py`
- Conda environment (`env.yml`) for quick setup
- Minimal CI (shellcheck + pyflakes)

Planned for v0.2.0: sample sheet orchestration, Nextflow preview, richer QC, container image.

---

## 🧭 Pipeline flow

```mermaid
flowchart LR
    A[FASTQ (R1/R2)] --> B[Align (bwa|bowtie2)]
    B --> C[Sort/Index (samtools)]
    C --> D[Pairing (pairtools)]
    D --> E[Contact matrix (.cool, cooler)]
    E --> F[QC collect (get_qc.py)]
    subgraph Outputs
      C
      D
      E
      F
    end
````

**I/O layout (per sample)**

```
input/{SAMPLE}_R1.fastq.gz
input/{SAMPLE}_R2.fastq.gz
results/{SAMPLE}/
  ├─ align/     # BAM + index
  ├─ pairs/     # .pairs / stats
  ├─ cool/      # .cool
  └─ logs/      # step-wise logs
qc/{SAMPLE}/
  ├─ QC.tsv
  └─ QC.html
```

---

## ⚙️ Requirements

* **System**: Linux on HPC (Slurm). Bash, Python ≥3.11.
* **Reference**: Genome FASTA (`mm10`/`hg38`, etc.) + index (`.fai`, bwa/bowtie2 index).
* **Tools** (installed via Conda/Mamba):

  * `bwa` *or* `bowtie2`, `samtools`, `pairtools`, `cooler`, `bedtools`, `pigz`
* **Python libs**: `pandas`, `matplotlib` (from `requirements.txt`)

> Tip: On Biowulf or similar, point to institute-provided FASTA when possible.

---

## 🚀 Quickstart

```bash
# 1) Clone
git clone git@github.com:Fuki-Kudoh/microC-pipeline.git
cd microC-pipeline

# 2) Environment (conda/mamba)
mamba env create -f env.yml
mamba activate microc

# 3) Put input FASTQs
ls input/
# SAMPLE_R1.fastq.gz  SAMPLE_R2.fastq.gz

# 4) Submit (example: mm10)
sbatch --time=24:00:00 --cpus-per-task=16 --mem=32g \
  mdp.sh SAMPLE mm10 /path/to/mm10.fa

# 5) Collect QC (after run finishes)
python get_qc.py -i results/SAMPLE -o qc/SAMPLE
```

Minimal Slurm resources that work for smoke runs:

* `--cpus-per-task 16` / `--mem 32g` / `--time 24:00:00`
  Adjust for your read depth and queue policy.

---

## 🧪 Sample-sheet (optional, early)

For multiple samples, you can fan out jobs with a tiny launcher:

```
samples.csv
sample,fastq1,fastq2,genome,fa
S1,input/S1_R1.fastq.gz,input/S1_R2.fastq.gz,mm10,/path/to/mm10.fa
```

```bash
bash run_samples.sh    # submits one Slurm job per row
```

---

## 🔧 Script interface

### `mdp.sh`

```text
Usage:
  mdp.sh <SAMPLE_ID> <GENOME> <GENOME_FASTA>

Env:
  INPUT_DIR   default: input
  OUT_DIR     default: results
  LOG_DIR     default: logs
```

**Behavior**

* Validates dependencies (bwa/bowtie2, samtools, pairtools, cooler, bedtools, pigz)
* Checks FASTA + `.fai`, and input FASTQs
* Writes step logs under `logs/` and `results/<SAMPLE>/logs/`
* Idempotent-ish: will not overwrite obvious final outputs unless `--force` (stub)

> `--dry-run` (stub) prints planned commands without running them.

### `get_qc.py`

```text
Usage:
  python get_qc.py -i results/<SAMPLE> -o qc/<SAMPLE>
```

Emits:

* `QC.tsv` – compact metrics (reads_total, mapped_rate, dup_rate, etc.)
* `QC.html` – quick text report + small plots (`reads_total.png`, …)

---

## 📦 Environments

### Conda/Mamba (recommended)

```bash
mamba env create -f env.yml
mamba activate microc
python -m pip install -r requirements.txt
```

### Manual (fallback)

```bash
mamba create -n microc python=3.11 bwa bowtie2 samtools pairtools cooler bedtools pigz -c conda-forge -c bioconda
mamba activate microc
pip install -r requirements.txt
```

> Pin additional versions in `env.yml` if you need stricter reproducibility.

---

## 🧰 Troubleshooting

* **`Permission denied (publickey)` when cloning**
  Configure SSH (e.g., 1Password SSH Agent) and add your public key on GitHub.

* **`FASTA not found` or missing `.fai`**
  Ensure `GENOME_FASTA` points to the actual `.fa/.fasta`; run `samtools faidx`.

* **Memory/timeouts**
  Start with `16 CPU / 32 GB / 24 h` and scale based on read depth.
  Very deep libraries may need `--cpus-per-task 32` and `--mem 64g+`.

* **Cooler errors**
  Verify chromosome naming consistency between FASTA and outputs (`chr1` vs `1`).

If a job fails, check `results/<SAMPLE>/logs/` and Slurm output files.
Please open an issue with the failing step, command, and log snippet.

---

## 📄 Citation & License

* **License**: MIT (see `LICENSE`)
* **How to cite**: see `CITATION.cff` (v0.1.0)

> If you use this pipeline in academic work, please reference the release tag and commit hash.

---

## 🗺️ Roadmap

* **v0.1.x (hardening)**: `--dry-run/--force`, better logs, sample-sheet UX, parse real QC logs, FAQs
* **v0.2.0 (preview)**: Nextflow wrapper (`slurm`/`local`), container image, benchmark notes
* **v0.3+**: richer QC HTML, automated reports, downstream modules

---

## 🧪 Repro note

To fully reproduce results, record:

* Commit: `git rev-parse --short HEAD`
* Release tag: `v0.1.0`
* `conda list --explicit` export
* FASTA source and index build command

---

## 🙌 Contributions

PRs and issues are welcome. Please keep changes small and focused.
For feature ideas: open a discussion with example inputs/outputs and expected runtime.

```
::contentReference[oaicite:0]{index=0}
```
