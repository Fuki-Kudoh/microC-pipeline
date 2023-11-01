# microC-dovetail-pipeline
This is just a pipeline that processes micro-C data by following the documentation by dovetail. \
https://micro-c.readthedocs.io/en/latest/

   **Command:**
   .. code-block:: console

    sbatch --time=24:00:00 --cpus-per-task=64 --mem=64g microC.sh <sample_ID> <genome> <genome location> [option]

    **Example:**
    .. code-block:: console

    sbatch --time=24:00:00 --cpus-per-task=64 --mem=64g microC.sh CTRL-1 mm10 /fdb/bwa/indexes/mm10.fa qc
