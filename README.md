# microC-pipeline
This is just a pipeline that processes micro-C data by following the documentation by dovetail. \
https://micro-c.readthedocs.io/en/latest/



**Command:**
   .. code-block::

    sbatch --time=24:00:00 --cpus-per-task=64 --mem=64g mdp.sh <sample_ID> <genome> <genome location>

**Example:**
    .. code-block::

    sbatch --time=24:00:00 --cpus-per-task=64 --mem=64g mdp.sh CTRL-1 mm10 /fdb/bwa/indexes/mm10.fa
