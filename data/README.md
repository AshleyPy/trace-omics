# Data inputs

TRACE-Omics begins from harmonised gene-level RNA-seq count matrices. Public source accessions are listed in `study_registry.csv`; raw biological data are intentionally excluded from Git.

For each analytical disease task, provide a genes-by-samples integer count matrix together with metadata containing at least:

```text
sample_id
participant_id
study
disease_group
outcome
repository_accession
```

`outcome` uses the exact values `control` and `case`, with `case` treated as the positive class. Sample identifiers must align exactly with count-matrix columns. Study and sample identifiers remain provenance variables and are excluded from the predictive feature set.

The final dissertation used distinct frozen partition manifests for the approximately 70:30 original/pooling analyses and approximately 80:20 study/matched-architecture analyses. Numerical reproduction requires those exact manifests rather than regeneration from a seed alone.

See `docs/DATA_CONTRACT.md` and the templates in `data/templates/`.
