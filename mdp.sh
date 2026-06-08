#! /bin/bash
set -euo pipefail

usage() {
    cat <<'USAGE'
Usage: mdp.sh <SAMPLE_ID> <GENOME_NAME> <GENOME_FASTA>

Arguments:
  SAMPLE_ID      Sample prefix for fastq/{SAMPLE_ID}_R1_001.fastq.gz and fastq/{SAMPLE_ID}_R2_001.fastq.gz
  GENOME_NAME    Short genome label, for example mm10 or hg38
  GENOME_FASTA   Path to the reference FASTA
USAGE
}

if [[ $# -ne 3 ]]; then
    usage >&2
    exit 2
fi

date

sample_ID="$1"
genome_name="$2"
genome="$3"
THREADS="${SLURM_CPUS_PER_TASK:-16}"

##############################
# This is for alignment.
# Example:
# sbatch --time=24:00:00 --cpus-per-task=16 --mem=64g mdp.sh <sample_ID> <genome_name> <genome_fasta>
#
# directory-+-fastq--+--{sample_ID}_R1_001.fastq.gz
#           |        +--{sample_ID}_R2_001.fastq.gz
#           +-get_qc.py
#           +-mdp.sh
##############################
error_handling() {
    echo "Error occurred in command: '${BASH_COMMAND}' at line: $1" >&2
    exit 1
}

require_command() {
    local command_name="$1"
    if ! command -v "${command_name}" >/dev/null 2>&1; then
        echo "Required command not found: ${command_name}" >&2
        exit 127
    fi
}
##############################

trap 'error_handling $LINENO' ERR

# Optional/site-specific module section. Edit these names for your HPC site if needed.
# The script does not assume module loading is available and does not install software.
if command -v module >/dev/null 2>&1; then
    module load fastqc trimgalore samtools bwa pairtools preseq juicer || true
fi

for command_name in fastqc trim_galore samtools bwa pairtools preseq juicer_tools bgzip pairix cooler python3; do
    require_command "${command_name}"
done

if [[ ! -f "${genome}" ]]; then
    echo "Genome FASTA not found: ${genome}" >&2
    exit 1
fi

fastq_r1="fastq/${sample_ID}_R1_001.fastq.gz"
fastq_r2="fastq/${sample_ID}_R2_001.fastq.gz"
trimmed_r1="fastq/${sample_ID}_R1_001_val_1.fq.gz"
trimmed_r2="fastq/${sample_ID}_R2_001_val_2.fq.gz"
chrom_sizes="genome/${genome_name}.chrom.sizes"
genome_fai="${genome}.fai"

if [[ ! -f "${fastq_r1}" ]]; then
    echo "R1 FASTQ not found: ${fastq_r1}" >&2
    exit 1
fi

if [[ ! -f "${fastq_r2}" ]]; then
    echo "R2 FASTQ not found: ${fastq_r2}" >&2
    exit 1
fi

mkdir -p "temp/${sample_ID}" BAM fastqc hic pairs stats genome cool logs
exec > >(tee -a "logs/${sample_ID}.log") 2>&1

printf 'Sample: %s\nGenome name: %s\nGenome FASTA: %s\nThreads: %s\n' "${sample_ID}" "${genome_name}" "${genome}" "${THREADS}"

if [[ ! -s "${genome_fai}" ]]; then
    genome_dir="$(dirname "${genome}")"
    if [[ -w "${genome_dir}" ]]; then
        samtools faidx "${genome}"
    else
        echo "FASTA index not found or empty: ${genome_fai}" >&2
        echo "The FASTA directory is not writable: ${genome_dir}" >&2
        echo "Please create the FASTA index first with: samtools faidx '${genome}'" >&2
        exit 1
    fi
fi
cut -f1,2 "${genome_fai}" > "${chrom_sizes}"

fastqc -t "${THREADS}" -o fastqc/ "${fastq_r1}" "${fastq_r2}"
trim_galore -j "${THREADS}" -o fastq --paired "${fastq_r1}" "${fastq_r2}"

# alignment -t is thread number, 30 min
bwa mem -5SP -T0 -t "${THREADS}" "${genome}" \
    "${trimmed_r1}" "${trimmed_r2}" \
    -o "temp/${sample_ID}/${sample_ID}.sam"

# parse --nproc-in and --nproc-out are thread numbers, 1.5 h
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
    --nproc-in "${THREADS}" --nproc-out "${THREADS}" --chroms-path "${chrom_sizes}" \
    "temp/${sample_ID}/${sample_ID}.sam" > "temp/${sample_ID}/${sample_ID}.pairsam"

# sort --nproc is thread number, 40 min
pairtools sort --nproc "${THREADS}" --tmpdir="temp/${sample_ID}" "temp/${sample_ID}/${sample_ID}.pairsam" \
    > "temp/${sample_ID}/${sample_ID}.sorted.pairsam"
rm "temp/${sample_ID}/${sample_ID}.pairsam"

# complexity session
# split 10 min
pairtools split --nproc-in "${THREADS}" --nproc-out "${THREADS}" \
    --output-sam "BAM/${sample_ID}.undedup.bam" "temp/${sample_ID}/${sample_ID}.sorted.pairsam"

# sort and indexing
samtools sort -@ "${THREADS}" -T "BAM/${sample_ID}.undedup.bam" \
    -o "BAM/${sample_ID}.undedup.PT.bam" "BAM/${sample_ID}.undedup.bam"
rm "BAM/${sample_ID}.undedup.bam"
samtools index "BAM/${sample_ID}.undedup.PT.bam"

# complexity
preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 \
    -output "stats/comp_${sample_ID}.txt" "BAM/${sample_ID}.undedup.PT.bam"
rm "BAM/${sample_ID}.undedup.PT.bam"

# Go back process
# dedup 50 min
pairtools dedup --nproc-in "${THREADS}" --nproc-out "${THREADS}" --mark-dups \
    --output-stats "stats/${sample_ID}.txt" \
    --output "temp/${sample_ID}/${sample_ID}.dedup.pairsam" "temp/${sample_ID}/${sample_ID}.sorted.pairsam"
rm "temp/${sample_ID}/${sample_ID}.sorted.pairsam"

# check stats
python3 get_qc.py -p "stats/${sample_ID}.txt" > "stats/qc_${sample_ID}.txt"

# split 10 min
pairtools split --nproc-in "${THREADS}" --nproc-out "${THREADS}" \
    --output-pairs "pairs/${sample_ID}.pairs" \
    --output-sam "BAM/${sample_ID}.bam" "temp/${sample_ID}/${sample_ID}.dedup.pairsam"
rm "temp/${sample_ID}/${sample_ID}.dedup.pairsam"

# sort and indexing
samtools sort -@ "${THREADS}" -T "BAM/${sample_ID}.bam" -o "BAM/${sample_ID}.PT.bam" "BAM/${sample_ID}.bam"
rm "BAM/${sample_ID}.bam"
samtools index "BAM/${sample_ID}.PT.bam"

# .hic contact map
juicer_tools pre -j "${THREADS}" "pairs/${sample_ID}.pairs" "hic/${sample_ID}.hic" "${genome_name}"
bgzip "pairs/${sample_ID}.pairs"

# .cool contact matrix
pairix "pairs/${sample_ID}.pairs.gz"
cooler cload pairix -p "${THREADS}" "${chrom_sizes}:1000" "pairs/${sample_ID}.pairs.gz" "cool/${sample_ID}.cool"
cooler zoomify -p "${THREADS}" "cool/${sample_ID}.cool"

#######################
# Version history and updating schedule
# v1.0 2023/11/01 From fastq to hic without qc
# v1.1 2023/11/02 add fastqc, galore, quality check
# v1.2 2024/01/02 add error handling
# v1.3 2024/01/04 Integrate cooler format
# v1.4 2024/07/19 bug fix
# v1.5 2024/11/01 change preseq timing
