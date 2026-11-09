# TRACE-Omics publication readiness checklist

These five `R/` files and the workflow form a **candidate patch**, not an independently validated release.

## Before copying

1. Keep a backup branch and compare each file against the original. Three files (`core.R`, `supplied_functions.R`, `utils_io.R`) were provided in full and edited. The other two (`study_models.R`, `supplied_architecture_functions.R`) were reconstructed/refactored from excerpts pasted in the chat. **Do not overwrite additional functions that may exist in the private originals**; merge changes after diff review.
2. Keep `legacy/` dissertation scripts unchanged. They remain the historical numerical reference.
3. Verify that `R/` names and the `tests/` scripts in this patch do not replace a test that your repository specifically relies on.
4. Ensure research-data sharing agreements permit every file in Git history. The generated count matrices, metadata, manifest files and prediction files must remain private unless explicitly approved.

## Automated checks

- `Rscript tests/parse_all.R`: parses all `.R` files; requires R but not analysis packages or datasets.
- `Rscript tests/test_contract_base.R`: basic seeded partitioning and labelling checks; uses base R only. This does **not** validate model predictions.
- `.github/workflows/r-parse.yml`: executes both scripts on every push/PR using the current R release.

These checks may expose pre-existing syntax errors elsewhere in `legacy/`; the check deliberately reports those errors rather than concealing them.

## Reproducibility validation

- Execute dependency-locked R environment (`renv.lock`) and separate Python environment for optional TabPFN / LightGBM; record Bioconductor and system library versions.
- Confirm exact frozen dissertation partitions and CV fold manifests, not just the seed.
- Exercise preprocessing under folds with deliberately permuted evaluation labels to verify evaluation labels cannot influence fitted transforms.
- Check feature selection, resampling, prediction and class-probability ordering at the sample level.
- Compare selected-gene sets, AUROC/BACC/AUPRC and other outputs against frozen dissertation references, with documented tolerances.
- Test on synthetic counts; verify that unseen cohorts fail explicitly for methods that cannot project them.
- Confirm the output root is outside the Git repository and no patient-level information exists in public Git history.

## Scientific scope

Current runners target binary gene-level bulk RNA-seq classification. T1D, TB, UC, RA and other diseases work only when input counts, phenotype definitions and cohort structure meet the documented contract. The generality of the code alone does not demonstrate successful model transportability.

The custom frozen ComBat-Seq procedure and some existing legacy transformations have **not** been audited for numerical equivalence to released package routines. Do not assert full reproduction or statistical correctness until this has been independently checked.
