
# TRACE-Omics: A Framework for Robust and Generalisable Transcriptomic Machine Learning

**A research framework for evaluating the robustness, reproducibility and cross-cohort generalisability of machine-learning models developed using high-dimensional transcriptomic data.**

**Author:** Ashley Ehlers  
**Institution:** Stellenbosch University, South Africa  
**Research area:** Bioinformatics | Computational Biology | Machine Learning | Transcriptomics  
**Implementation:** R

---

## Overview

TRACE-Omics is a research framework developed to investigate how analytical decisions, dataset characteristics and validation strategies influence the performance and generalisability of machine-learning models trained on transcriptomic data.

Although machine-learning algorithms can achieve high predictive accuracy during internal validation, their performance may deteriorate when applied to independent cohorts, populations or experimental conditions. This discrepancy raises important questions about the reliability and transportability of transcriptomic prediction models.

TRACE-Omics addresses these challenges by examining how preprocessing, feature selection, classifier specification, class imbalance management and validation design influence predictive performance.

The framework was developed as part of doctoral research in Bioinformatics and Computational Biology at Stellenbosch University.

### Central research question

**Under what conditions do machine-learning models trained on transcriptomic data produce robust and generalisable predictions across independent datasets and study populations?**

The research emphasises that internal predictive accuracy alone is insufficient evidence of model reliability. Instead, model evaluation should consider performance stability, cohort heterogeneity, analytical sensitivity and generalisability beyond the development dataset.

---

## Research objectives

TRACE-Omics was developed to:

1. Evaluate the influence of analytical pipeline choices on transcriptomic machine-learning performance.
2. Investigate how preprocessing and feature-selection strategies affect predictive performance and feature stability.
3. Compare machine-learning classifiers under different analytical and sampling conditions.
4. Assess the influence of class imbalance and cohort composition on model evaluation.
5. Examine cross-cohort generalisability using independent study populations.
6. Distinguish internal predictive performance from evidence of transportability.
7. Support transparent, reproducible and scientifically defensible machine-learning workflows.

---

## Research scope

The broader doctoral research investigated transcriptomic machine learning across five infectious and inflammatory disease contexts:

- Tuberculosis (TB)
- Ulcerative colitis (UC)
- Rheumatoid arthritis (RA)
- Crohn's disease (CD)
- Type 1 diabetes (T1D)

The analytical framework incorporated multiple datasets, preprocessing approaches, feature-selection methods, classification algorithms and validation strategies.

### Analytical components

| Component | Approaches investigated |
|---|---|
| Data preprocessing | Variance-stabilising transformation and transcriptomic harmonisation strategies |
| Batch-effect management | SVA, ComBat and related approaches |
| Feature selection | Statistical, information-theoretic, model-based and hybrid methods |
| Classification | Random forests, gradient boosting, support vector machines and other supervised learners |
| Class imbalance | SMOTE, upsampling, downsampling and alternative balancing strategies |
| Hyperparameter optimisation | Nested cross-validation and model-specific optimisation |
| Generalisability | Cross-cohort validation and leave-one-cohort-out evaluation |
| Performance assessment | AUROC, AUPRC, balanced accuracy, calibration and related measures |
| Feature stability | Feature-set agreement and stability assessment |

**Important:** These components describe the broader TRACE-Omics doctoral research. Not all components are currently exposed through the public repository's executable pipeline.

---

## Public repository

This repository contains a public-facing implementation of selected TRACE-Omics components, including R source files, configuration-driven data handling, automated checks and an executable workflow for evaluating a pre-specified model pipeline.

The public implementation is intended to support methodological transparency and further software development.

It should not be interpreted as a complete executable reproduction of every experiment conducted during the dissertation.

### Repository structure

```text
trace-omics/
|
|-- R/
|   |-- core.R
|   |-- supplied_functions.R
|   |-- utils_io.R
|   |-- study_models.R
|   `-- supplied_architecture_functions.R
|
|-- config/
|   `-- example.yml
|
|-- scripts/
|   `-- run_trace.R
|
|-- tests/
|   |-- parse_all.R
|   `-- test_contract_base.R
|
|-- docs/
|   `-- PUBLICATION_CHECKLIST.md
|
|-- .github/
|   `-- workflows/
|       `-- r-parse.yml
|
|-- README.md
|-- CONTRIBUTING.md
`-- LICENSE
```

### Main components

**R/core.R**

Provides the common data-loading interface and configuration-driven functionality.

**R/supplied_functions.R**

Contains supporting R functions provided as part of the public implementation.

**R/utils_io.R**

Provides input/output and related utility functionality.

**R/study_models.R**

Contains study-related modelling code included in the public release.

**R/supplied_architecture_functions.R**

Contains additional architecture-related functions included in the release.

**config/example.yml**

Illustrates the configuration structure required to specify input data and analysis settings.

**scripts/run_trace.R**

Provides a command-line entry point for executing one pre-specified pipeline.

**tests/**

Contains syntax and contract-level checks.

**docs/PUBLICATION_CHECKLIST.md**

Documents additional scientific validation requirements for reproducibility and publication.

---

## Data requirements

TRACE-Omics uses a configuration-driven interface designed to accommodate transcriptomic datasets without requiring disease-specific identifiers in the common R core.

The implementation expects input information specifying:

- A gene-expression count matrix.
- Sample-level metadata.
- Canonical sample and study identifiers.
- Outcome labels and, where necessary, explicit class mappings.
- A sample partition manifest.
- An appropriate output directory.

The common metadata interface uses the following canonical fields:

| Field | Description |
|---|---|
| `sample_id` | Unique sample identifier |
| `participant_id` | Participant identifier |
| `study` | Dataset or cohort identifier |
| `outcome` | Outcome or classification label |

Explicit column mappings can be used when the original dataset follows a different naming convention.

### Data independence

The common R core does not require a particular disease label or dataset name.

Although the doctoral research investigated TB, UC, RA, CD and T1D, the configuration interface is not restricted to those disease contexts.

Users are responsible for ensuring that their data satisfy the required structure, outcome definitions and analytical assumptions.

### Data privacy

Input datasets and analysis outputs should be stored outside the Git repository.

Do not commit:

- Restricted patient-level information.
- Personally identifiable information.
- Authentication credentials.
- Confidential institutional datasets.
- Private local configuration files containing sensitive information.

All datasets must be accessed and used in accordance with their applicable licences, ethical approvals and data-sharing conditions.

---

## Installation and requirements

### 1. Clone the repository

```bash
git clone https://github.com/AshleyPy/trace-omics.git
cd trace-omics
```

### 2. Prepare the R environment

An installation of R is required.

The configuration loader requires the `yaml` package when using YAML configuration files.

```r
install.packages("yaml")
```

Additional R packages may be required depending on the selected modelling pipeline and associated functions.

For numerical reproducibility, users should record and, where possible, restore the relevant R version and package versions.

The public repository does not, by itself, guarantee a fully locked computational environment.

### 3. Configure the input data

Use `config/example.yml` as a reference for constructing a project-specific configuration file.

For example:

```text
config/my_task.yml
```

The example configuration contains placeholder paths and is not expected to execute without modification.

Specify valid input paths and ensure that the output directory is outside the repository.

---

## Running TRACE-Omics

### Option 1: Load a dataset in R

From the repository root:

```r
source("R/core.R")

d <- trace_load("config/my_task.yml")

dim(d$counts)

table(
  d$metadata$study,
  d$metadata$outcome
)
```

This example loads configured data and allows users to inspect the input structure.

Subsequent modelling requires an explicit and valid pipeline specification.

### Option 2: Execute a pre-specified pipeline

From a terminal or command prompt:

```bash
Rscript scripts/run_trace.R --config config/my_task.yml
```

The configuration must provide:

1. Valid input data paths and metadata mappings.
2. A pre-specified modelling pipeline.
3. A validated frozen partition manifest.
4. Valid output settings.

The runner trains the specified pipeline using development samples and evaluates the untouched holdout once.

### Important implementation limitations

The public runner:

- Executes one pre-specified pipeline.
- Uses an explicitly defined development/holdout partition.
- Does not implement the complete nested hyperparameter optimisation workflow.
- Does not reproduce every TRACE-Omics dissertation experiment.
- Does not establish external transportability merely by reporting holdout performance.

A successful execution should therefore not be interpreted as complete reproduction of the doctoral research findings.

---

## Partition management and data leakage prevention

Reliable machine-learning evaluation requires strict separation between model development and final evaluation.

TRACE-Omics uses a separate sample partition manifest to define the relevant analysis partitions.

### Frozen partitions

The original dissertation analyses used specific, frozen sample partitions.

Numerical reproduction requires the corresponding validated dissertation partition manifests and authorised input data.

### Generating new partitions

A new stratified partition can be generated only when explicitly requested using:

```yaml
partition_mode: generate
```

For newly generated, non-reference partitions, stratification can be configured using `partition_strata`.

**A newly generated partition is not equivalent to a frozen dissertation partition and must not be presented as reproducing the original results.**

### Participant-level considerations

Where appropriate and scientifically justified, the configuration can specify:

```yaml
one_sample_per_participant: true
```

This setting should be used only after verifying that it is appropriate for the input data and intended analysis.

Participant and study identifiers must not be inadvertently introduced into the predictive feature space.

All data-adaptive modelling operations should respect the development/evaluation boundary.

---

## Testing and continuous integration

The repository includes automated checks intended to identify syntax problems and selected implementation-level issues.

### Run syntax checks

```bash
Rscript tests/parse_all.R
```

### Run contract tests

```bash
Rscript tests/test_contract_base.R
```

### GitHub Actions

The repository includes:

```text
.github/workflows/r-parse.yml
```

This workflow supports automated checking of the R source.

### Interpretation of test results

Passing syntax checks and contract tests does not independently demonstrate:

- Numerical reproducibility.
- Correct implementation of every statistical method.
- Absence of all forms of data leakage.
- Equivalence with dissertation results.
- Clinical validity.
- Generalisability to unseen study populations.

These properties require additional scientific evaluation.

---

## Reproducibility and scientific validation

TRACE-Omics distinguishes between software execution, computational reproducibility and evidence of predictive generalisability.

### Software execution

The public code can be assessed using its available syntax, contract and workflow checks.

### Computational reproducibility

Numerical reproduction requires the appropriate combination of:

- Authorised input datasets.
- Frozen sample partition manifests.
- Matching analysis specifications.
- Reproducible software dependencies.
- Consistent random seeds where applicable.
- Prediction-level and metric-level comparisons against reference outputs.

### Generalisability

Performance on a development dataset or an internal holdout does not independently establish transportability.

Evidence of transportability requires evaluation under clearly defined distributional changes, such as independent study cohorts or distinct target populations.

These evaluations must preserve the intended prediction task and avoid using information from the evaluation cohort during model development unless explicitly permitted by the study design.

For additional validation requirements, see:

`docs/PUBLICATION_CHECKLIST.md`

---

## Scientific interpretation

The central motivation for TRACE-Omics is that high internal predictive accuracy does not necessarily imply robust or transportable performance.

Analytical decisions can influence the apparent predictive signal, selected features and estimated model performance.

Consequently, the interpretation of transcriptomic machine-learning results should consider:

- The relationship between the training and target populations.
- Differences in cohort composition.
- The influence of preprocessing and harmonisation.
- Feature-selection stability.
- Sensitivity to model specification.
- Class imbalance and outcome prevalence.
- Validation design and potential information leakage.
- Discrimination and calibration under external evaluation.

The broader research seeks to strengthen the methodological foundations of transcriptomic prediction rather than promote a particular algorithm as universally superior.

---

## Current limitations

The following limitations should be considered when using this repository:

1. The public implementation does not contain every component of the complete TRACE-Omics dissertation framework.
2. Original numerical results cannot be reproduced without the required frozen manifests, data and software specifications.
3. Example configurations require user-supplied paths and appropriate datasets.
4. The public runner executes a single pre-specified pipeline rather than the complete nested optimisation procedure.
5. Automated software checks are not substitutes for statistical or scientific validation.
6. The public implementation is a research resource and is not validated for clinical deployment.

These limitations are documented to support transparent interpretation and responsible reuse.

---
