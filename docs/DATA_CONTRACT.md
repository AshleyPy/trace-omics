# Input data contract

The computational starting point is a harmonised gene-level count matrix. Raw-read alignment is outside this package. Public accession identifiers appear in Chapter 3; public availability alone does not identify the exact retained participant subset or gene mapping used in the dissertation.

## Matrix

Save a numeric matrix as RDS with genes as rows, samples as columns, unique gene identifiers and unique sample identifiers. Counts must be finite non-negative integers. Supply the same gene universe and ordering used by the intended experiment. The main pooled comparison used 4,450 genes. Each column must have a positive library total.

## Metadata

Supply CSV columns:

| Column | Meaning |
|---|---|
| sample_id | Unique analysis-record ID, matching one count-matrix column |
| participant_id | Globally harmonised participant identity across all contributing studies/tasks |
| study | Study label, such as TB-1 |
| disease_group | TB, RA, UC, CD or T1D |
| disease | Exact binary label: control or case |

The input adapter creates the literal `outcome` and `cohort_group` fields expected by surviving functions. Scientific outputs retain disease status and study terminology. Source labels must be harmonised explicitly; the driver does not guess whether progression labels, treatment timepoints or ambiguous phenotypes qualify as binary cases.

The portable driver requires one retained record per participant. Resolve longitudinal samples and shared controls using a documented rule before running. Keep the crosswalk outside the public repository when its sharing conditions require this.

## Partition manifest

Supply a CSV containing `sample_id,partition`, with one row per input sample and partition set to `development` or `holdout`. Participants must stay within a single partition. Restore original manifests for historical reproduction. A newly generated split with the same seed can differ when row order, preprocessing or software changes.

For pooled Model A, remove the TB development records while retaining the identical shared evaluation population. Model B contains all five development populations. Record any deduplication before making this comparison. For LODO, verify participant exclusion across disease-specific records as well as exclusion of the disease label.

## Frozen prediction file

Stage 04 accepts one system and one evaluation population per CSV with `sample_id,disease,probability,threshold`. Probabilities must be case probabilities between zero and one. The threshold must be a single value chosen during development. The script verifies sample IDs and labels against the holdout manifest. Prepare a separate file per classifier when using the architecture driver’s long-format predictions.

## Computational state

For complete replay retain input hashes, retained genes, mappings to syntactic feature names, source/fold manifests, fitted transformation objects, selected features, fitted classifiers, tuning archives, thresholds, per-observation predictions, package versions and checkpoint identifiers. Keep source counts available when a supplied function’s saved state is only a summary of its fitted transformation.
