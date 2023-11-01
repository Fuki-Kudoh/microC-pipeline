#! /bin/bash

echo ${date}
sample_ID=$1
genome_name=$2
genome=$3
option=$4

##############################
#This is for alignment.
#sbatch --time=24:00:00 --cpus-per-task=64 --mem=64g microC.sh <sample_ID> <genome> <genome location> [option]
#directory-+-fastq--+--{sample_ID}_R1_001.fastq.gz
#          |        +--{sample_ID}_R2_001.fastq.gz
#          +-BAM
#          +-fastqc
#          +-hic
#          +-pairs
#          +-stats
#          +-temp
##############################

module load samtools bwa pairtools preseq juicer
mkdir temp/${sample_ID}

trim_galore -j 4 --fastqc_args "--outdir fastqc" --paired fastq/${sample_ID}_R1.001.fastq.gz fastq/${sample_ID}_R2.001.fastq.gz

#alignment -t is thread number, 30min
bwa mem -5SP -T0 -t64 ${genome}.fa \
fastq/${sample_ID}_R1_001.fastq.gz fastq/${sample_ID}_R2_001.fastq.gz \
-o temp/${sample_ID}/${sample_ID}.sam

#parse --nproc-in and --nproc-out is thread number, 1.5h
pairtools parse --min-mapq 40 --walks-policy 5unique --max-inter-align-gap 30 \
--nproc-in 64 --nproc-out 64 --chroms-path ${genome}.fa \
temp/${sample_ID}/${sample_ID}.sam > temp/${sample_ID}/${sample_ID}.pairsam
#sort --nproc is thread number. 40 min
rm temp/${sample_ID}/${sample_ID}.sam

pairtools sort --nproc 64 --tmpdir=temp/${sample_ID} temp/${sample_ID}/${sample_ID}.pairsam \
> temp/${sample_ID}/${sample_ID}.sorted.pairsam
rm temp/${sample_ID}/${sample_ID}.pairsam
#dedup 50 min
pairtools dedup --nproc-in 64 --nproc-out 64 --mark-dups --output-stats stats/${sample_ID}.txt \
--output temp/${sample_ID}/${sample_ID}.dedup.pairsam temp/${sample_ID}/${sample_ID}.sorted.pairsam
rm temp/${sample_ID}/${sample_ID}.sorted.pairsam
#check stats
#wget https://github.com/dovetail-genomics/Micro-C/blob/973c67ae2be329ded96007cddeec39667f60698e/get_qc.py
#python3 get_qc.py -p stats/${sample_ID}.txt > stats/qc_${sample_ID}.txt

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

#cooltools
#module load cooltools
#bgzip pairs/${sample_ID}.pairs
#pairtools sort --nproc 64 -o pairs/${sample_ID}.sorted.pairs.gz pairs/${sample_ID}.pairs.gz
#rm pairs/${sample_ID}.pairs
#cooler cload pairix -p 64 ${genome}:1000 pairs/${sample_ID}.sorted.pairs.gz cooler/${sample_ID}.cool

#######################
# Version history and updating schedule
# v1.0 2023/11/01 From fastq to hic without qc
# v1.1 TBD        add Galore and fastqc
# v1.2 TBD        add quality check
# v1.3
