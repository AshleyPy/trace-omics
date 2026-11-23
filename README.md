# TRACE-Omics: reproducibility and CI patch

This patch contains `R/core.R`, `R/supplied_functions.R`, `R/utils_io.R`,
`R/study_models.R`, `R/supplied_architecture_functions.R`,
`.github/workflows/r-parse.yml`, `tests/parse_all.R`, `tests/test_contract_base.R`,
`config/example.yml`, `scripts/run_trace.R` and `docs/PUBLICATION_CHECKLIST.md`.

## What is generic?

The common R core uses count matrix and metadata paths supplied through a YAML/RDS
configuration, canonical `sample_id`, `participant_id`, `study` and `outcome` fields
(with explicit column mappings), optional class mappings, and a separate sample
partition manifest. Paths are checked and outputs must be written outside the
repository. Nothing in the core requires the disease to be CD, TB, UC, RA or T1D.

The frozen dissertation manifests must be supplied for numerical reproduction.
A new stratified split can be generated only when explicitly requested with
`partition_mode: generate`, and is **not** equivalent to the dissertation split.

## How to check in RStudio or from the command line

From the repository root (after merging):

```sh
Rscript tests/parse_all.R
Rscript tests/test_contract_base.R
```

To load a task from a YAML configuration (requires `yaml`):

```r
source("R/core.R")
d <- trace_load("config/my_task.yml")
# inspect: dim(d$counts), table(d$metadata$study, d$metadata$outcome)
# next: use trace_predict_pipeline() with an explicit pipeline specification
```

The attached `config/example.yml` is illustrative: its placeholder paths will
not exist on another computer. Set your own private absolute/relative paths,
with an output directory outside the Git repository.

## Migration advice

Commit these changes on a **new branch**, check the diff, then merge after
`Rscript tests/parse_all.R` and the tests pass. Do not overwrite the private
full originals of `study_models.R` and `supplied_architecture_functions.R`
without confirming that no additional source has been lost.

## Limits

Static syntax and smoke checks do not prove computational reproducibility.
Machine-learning results require packages, a locked environment, authorised
input datasets and comparison against the dissertation's frozen outputs.
See `docs/PUBLICATION_CHECKLIST.md` for the scientific validation tasks.

## Run one pre-specified pipeline

```sh
Rscript scripts/run_trace.R --config config/my_task.yml
```

A YAML file must provide a **pre-specified** `pipeline` and a validated frozen
partition manifest. The runner trains a single pipeline on development samples
and evaluates the untouched holdout once. It does **not** implement nested
optimisation or every TRACE-Omics dissertation component. If needed, set
`one_sample_per_participant: true` explicitly only when verified. For a new
non-reference split, `partition_strata` is configurable.
