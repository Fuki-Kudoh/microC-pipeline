#! /bin/bash

echo ${date}
sample_ID=$1
genome_name=$2
genome=$3

##############################
#This is for alignment.
#sbatch --time=24:00:00 --cpus-per-task=64 --mem=64g mdp.sh <sample_ID> <genome> <genome location>
#directory-+-fastq--+--{sample_ID}_R1_001.fastq.gz
#          |        +--{sample_ID}_R2_001.fastq.gz
#          +-get_qc.py
#          +-mdp.sh
##############################
error_handling() {
    echo "Error occurred in command: '$BASH_COMMAND' at line: $1"
    exit 1
}
##############################

module load fastqc trimgalore samtools bwa pairtools preseq juicer
mkdir temp temp/${sample_ID} BAM fastqc hic pairs stats genome cool

trap 'error_handling $LINENO' ERR

fastqc -t 64 -o fastqc/ fastq/${sample_ID}_R1_001.fastq.gz fastq/${sample_ID}_R2_001.fastq.gz
trim_galore -j 64 -o fastq --paired fastq/${sample_ID}_R1_001.fastq.gz fastq/${sample_ID}_R2_001.fastq.gz

#alignment -t is thread number, 30min
bwa mem -5SP -T0 -t64 ${genome} \
fastq/${sample_ID}_R1_001_val_1.fq.gz fastq/${sample_ID}_R2_001_val_2.fq.gz \
-o temp/${sample_ID}/${sample_ID}.sam

#parse --nproc-in and --nproc-out is thread number, 1.5h
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
--nproc-in 64 --nproc-out 64 --chroms-path ${genome} \
temp/${sample_ID}/${sample_ID}.sam > temp/${sample_ID}/${sample_ID}.pairsam
#sort --nproc is thread number. 40 min

pairtools sort --nproc 64 --tmpdir=temp/${sample_ID} temp/${sample_ID}/${sample_ID}.pairsam \
> temp/${sample_ID}/${sample_ID}.sorted.pairsam
rm temp/${sample_ID}/${sample_ID}.pairsam
#dedup 50 min
pairtools dedup --nproc-in 64 --nproc-out 64 --mark-dups --output-stats stats/${sample_ID}.txt \
--output temp/${sample_ID}/${sample_ID}.dedup.pairsam temp/${sample_ID}/${sample_ID}.sorted.pairsam
rm temp/${sample_ID}/${sample_ID}.sorted.pairsam
#check stats
touch stats/qc_${sample_ID}.txt
python3 get_qc.py -p stats/${sample_ID}.txt > stats/qc_${sample_ID}.txt

#split 10 min
pairtools split --nproc-in 64 --nproc-out 64 --output-pairs pairs/${sample_ID}.pairs \
--output-sam BAM/${sample_ID}.bam temp/${sample_ID}/${sample_ID}.dedup.pairsam
rm temp/${sample_ID}/${sample_ID}.dedup.pairsam

#sort and indexing
samtools sort -@64 -T BAM/${sample_ID}.bam -o BAM/${sample_ID}.PT.bam BAM/${sample_ID}.bam
rm BAM/${sample_ID}.bam
samtools index BAM/${sample_ID}.PT.bam
#complexity
preseq lc_extrap -bam -pe -extrap 2.1e9 -step 1e8 -seg_len 1000000000 \
-output stats/comp_${sample_ID}.txt BAM/${sample_ID}.PT.bam
#.hic contact map
juicer_tools pre -j 64 pairs/${sample_ID}.pairs hic/${sample_ID}.hic ${genome_name}
bgzip pairs/${sample_ID}.pairs

#.cool contact matrix
samtools faidx ${genome} --output genome/${genome_name}.fa.faidx
cut -f1,2 genome/${genome_name}.fa.fai > genome/${genome_name}.genome

module load mamba_install
mamba_install
mamba_install update --install
source myconda
conda install -c conda-forge -c bioconda cooler
pairix pairs/${sample_ID}.pairs.gz
cooler cload pairix -p 64 genome/${genome_name}.genome:1000 pairs/${sample_ID}.pairs.gz cool/${sample_ID}.cool
cooler zoomify -p 64 cool/${sample_ID}.cool

#######################
# Version history and updating schedule
# v1.0 2023/11/01 From fastq to hic without qc
# v1.1 2023/11/02 add fastqc, galore, quality check
# v1.2 2024/01/02 add error handling
# v1.3 2024/01/04 Integrate cooler format
# v1.4 2024/07/19 bug fix
