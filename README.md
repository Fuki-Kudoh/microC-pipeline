# microC-dovetail-pipeline
This is just a pipeline that processes micro-C data by following the documentation by dovetail. \
https://micro-c.readthedocs.io/en/latest/ 

   **Command:**
   .. code-block:: console
   
    sbatch --time=24:00:00 --cpus-per-task=64 --mem=64g microC.sh <sample_ID>

Set paired fastq files as below. \
directory---fastq--+--{sample_ID}_R1.fastq.gz \
                   +--{sample_ID}_R2.fastq.gz


directory-+-fastq--+--{sample_ID}_R1.fastq.gz \
          |        +--{sample_ID}_R2.fastq.gz \                
          +-BAM    \
          +-hic    \
          +-pairs  \
          +-stats  \
          +-temp

