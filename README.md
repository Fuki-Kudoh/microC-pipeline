# microC-pipeline

Preprocessing, QC, and visualization(.hic and .cool) of Micro-C

## Quickstart
```bash
# dependency（conda）
mamba create -n microc \
  python=3.11 bwa samtools pairtools cooler bedtools pigz -c conda-forge -c bioconda
mamba activate microc

# Run
sbatch --time=24:00:00 --cpus-per-task=16 --mem=32g \
  mdp.sh SAMPLE_ID mm10 /path/to/mm10.fa

```
## Requirements

bash, python3, Slurm

tools: bwa, samtools, pairtools, cooler, bedtools, pigz

reference: FASTA + index（Ex: mm10 / hg38）

## Inputs / Outputs

Input: {SAMPLE}_R1.fastq.gz, {SAMPLE}_R2.fastq.gz  under input/ directory

Output: results/{SAMPLE}/

QC: python get_qc.py -i results/{SAMPLE} -o qc/{SAMPLE} → QC.tsv, QC.html
