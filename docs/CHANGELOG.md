# TRACE-Omics publication patch: change log

## Source provenance

- `R/core.R`: replaced the supplied 186-line orchestration source with a configurable task loader, participant-level split safeguards, explicit project-root resolution, cohort-aware CV, protected output paths, optional label mappings, recorded input checksums and a generic prediction interface. This is a functional refactor requiring numerical equivalence tests.
- `R/supplied_functions.R`: copied from the supplied complete 2,045-line file, retained the historical core transformations and FS procedures, and inserted input checks/clarifying comments. Do not claim statistical equivalence to external packages without numerical tests.
- `R/utils_io.R`: copied from the supplied complete 836-line file, corrected invalid R string escapes in regex patterns, replaced `dplyr`/`stringr` label mapping with base-R equivalents, tightened count and prediction checks, expanded forbidden predictor metadata, and respected positive class settings.
- `R/study_models.R`: reconstructed/refactored from the 30-line excerpt, with parameterised random forest, deterministic screening, outcome validation and development-only threshold selection. Compare against the **complete private original** before overwriting.
- `R/supplied_architecture_functions.R`: reconstructed/refactored from the excerpt, with explicit RF settings, lazy Python imports, checked positive-probability columns, configurable compute settings, and the reference 10-row LightGBM grid. Compare against the **complete private original** before overwriting.
- `.github/workflows/r-parse.yml`: use GitHub's currently supported checkout and R setup actions and run syntax checks plus base-R synthetic smoke tests. No patient data or model libraries required for this stage.
- `tests/parse_all.R`: parse every R file and show exactly which filename failed.
- `tests/test_contract_base.R`: synthetic partition, fold, label, input and audit tests. Requires R but no BioConductor/ML packages.
- `scripts/run_trace.R`: minimal configuration-driven pipeline demonstration with secure holdout outputs; deliberately not the full nested-CV dissertation framework.

## Validation record

- Static lexical scanner (unrecognised escapes, balanced delimiters, terminated string literals): **PASS** for each included R file.
- YAML parsed in Python: **PASS** for the GitHub Actions workflow and example configuration.
- Real `Rscript parse()` and R execution: **NOT RUN** in this editing environment because R is not installed.
- Statistical, predictive and cross-disease numerical equivalence: **NOT TESTED**.

A green GitHub Actions parse check is necessary but not sufficient for publishing a research framework with validated predictive performance.
