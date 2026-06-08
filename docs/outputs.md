# Output notes

This document describes current legacy outputs and planned future output organization. It is a planning document, not a guarantee that the modern output layout has been implemented.

## Current legacy outputs

The retained `mdp.sh` workflow writes outputs and intermediates into repository-root directories such as:

```text
BAM/
cool/
fastqc/
hic/
genome/
pairs/
stats/
temp/
logs/
```

The uppercase `BAM/` path is intentional: it matches the historical legacy script layout. These paths reflect the legacy/minimal Slurm workflow and may include intermediates that are removed after downstream products are created.

## Planned future output direction

A future production-facing pipeline should make final, temporary, and diagnostic outputs explicit. A candidate future layout is documented in `docs/roadmap.md`, but it is not implemented yet.

Expected future documentation should define:

- final contact-pair outputs and indexes
- contact matrix outputs
- optional BAM retention
- run metadata
- QC tables, JSON, and reports
- temporary files and cleanup behavior
