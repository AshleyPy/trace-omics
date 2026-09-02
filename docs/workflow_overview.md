TRACE-Omics workflow overview

TRACE-Omics stands for Transcriptomic Robustness Assessment and Cross-cohort Evaluation for Omics learning.

TRACE-Omics is the reproducible analytical framework developed for this dissertation to evaluate when transcriptomic machine-learning relationships remain reliable beyond the data used to estimate them. The framework treats generalisability as an empirical property of the complete learning system and evaluates it through controlled changes in sample support, transcriptomic representation, disease composition, cohort composition, and model architecture.

The framework was developed using multicohort RNA-seq data from tuberculosis (TB), rheumatoid arthritis (RA), ulcerative colitis (UC), Crohn's disease (CD), and type 1 diabetes (T1D).

Contents

Conceptual framework

Disease datasets

Input data requirements

Transcriptomic representations

Feature selection and feature-space characterisation

Experiment 1: Sample-support perturbation

Experiment 2: Disease-composition perturbation

Experiment 3: Cohort-composition perturbation

Experiment 4: Model-architecture perturbation

Performance metrics

Reproducibility controls

Repository structure

Recommended execution logic

Computational environment

Public repository guidance

Citation

Conceptual framework

TRACE-Omics is organised into four sequential stages:

Representation specification
        |
        v
Feature-set characterisation
        |
        v
Model development and optimisation
        |
        v
Generalisation and reliability assessment

TRACE-Omics stage

Core question

Main evidence generated

Representation specification

What statistical form of the transcriptome is presented to the learner?

PCA, cohort/outcome silhouette structure, residual sequencing-depth dependence

Feature-set characterisation

Is the selected molecular information reproducible and statistically well conditioned?

Jaccard agreement, VIF, effective rank, bootstrap uncertainty

Model development and optimisation

How strongly and consistently can the development data support prediction?

Repeated CV, nested CV, class-handling comparisons, optimisation summaries

Generalisation and reliability assessment

Does the fitted relationship persist when the target distribution changes?

Frozen held-out evaluation, LODO, LOCO, calibration, permutation importance, ALE, architecture comparison

The framework separates internal learnability from transportability and distinguishes predictive reliability from molecular reliability.

Disease datasets

TRACE-Omics uses publicly available RNA-seq datasets from NCBI GEO and EMBL-EBI resources.

Disease

Cohorts

Final analytical population

TB

GC6 (GSE94438), ACS (GSE79362), Catalysis (GSE89403)

465

RA

Sales (E-MTAB-13361), Shchetynsky (GSE90081), Guo (GSE120178)

230

UC

Wacker (GSE177044), Halfvarson (E-MTAB-11349), Mo (GSE112057)

1,253

CD

Halfvarson (E-MTAB-11349), Mo (GSE112057)

495

T1D

Leal (GSE123658), Speake (GSE60424)

55

TB required participant-level curation because GC6 and ACS contained longitudinal progression data and Catalysis contained repeated treatment measurements. One eligible observation per participant was retained.

UC and CD were constructed as separate binary disease-versus-control tasks even where they originated from the same source studies.

Most cohorts used whole-blood RNA-seq. The Sales RA cohort contained monocyte-derived transcriptomic profiles and was retained as an explicitly documented biological shift.

Raw expression data are not distributed with this repository. Users should obtain the public datasets from their original repositories and reconstruct the expected local input structure described below.

Input data requirements

Metadata

The common analytical framework requires:

Column

Description

sample_id

Unique biological-sample identifier

cohort

Source study or cohort identifier

outcome

Binary outcome encoded as case or control

Repository accession identifiers should also be retained as provenance fields where available.

Count data

TRACE-Omics expects genes in rows and samples in columns after import.

Example:

data/raw/TB/
├── counts/
│   ├── GC6_counts.csv
│   ├── ACS_counts.csv
│   └── Catalysis_counts.csv
└── metadata.rds

Each count table should contain a gene-symbol field and one column per sample. Metadata sample identifiers must match count-matrix column names.

Gene harmonisation

Independent cohorts may use Ensembl, Entrez, or other gene identifiers. TRACE-Omics converts cohort-specific identifiers to HGNC symbols, identifies the genes shared across cohorts required for a given analysis, aligns matrices to an identical gene order, and preserves the harmonised raw-count matrix as the canonical molecular input.

Alternative transcriptomic representations are derived from this canonical object.

Transcriptomic representations

The general representation stage evaluates six candidate molecular inputs:

Raw counts
VST only
VST + ComBat
VST + ComBat-Ref
VST + fSVA
ComBat-Seq + log1p

Data-dependent transformations are fitted from development or source observations according to the corresponding experiment.

Representation quality is characterised using:

principal component analysis (PCA);

cohort- and outcome-based silhouette analysis;

residual sequencing-depth dependence.

Feature selection and feature-space characterisation

The general feature-selection stage evaluates eight methods:

ANOVA
AUC
RF importance
Boruta
CMIM
mRMR
Information gain
RF variance

The default comparison retains the 50 highest-ranked genes from each selector.

Feature-set agreement is quantified using pairwise Jaccard similarity. Statistical conditioning is characterised using mean VIF, maximum VIF, and effective rank, with bootstrap uncertainty where specified.

Experiment 1: Sample-support perturbation

Question

Does useful transcriptomic prediction remain recoverable when different diseases are constrained to the same small learning budget?

Design

TB, RA, UC, and CD were each cohort-stratified to:

41 cases + 41 controls = 82 observations

Each equal-size population was divided into:

Development: 66 observations (33 case, 33 control)
Held-out:    16 observations (8 case, 8 control)

T1D was excluded from this experiment because its available sample size was below the common 82-observation budget.

Analytical comparison

Three representations were carried forward:

VST + fSVA
VST + ComBat-Ref
ComBat-Seq + log1p

Five feature-selection methods were evaluated:

Information gain
RF importance
Boruta
mRMR
CMIM

The benchmark included:

Featureless
Naive Bayes
Linear SVM
RBF SVM
C5.0
Random Forest
LightGBM
XGBoost

Repeated cross-validation used 10 folds and 2 repeats. Representation fitting and feature selection were repeated inside each training fold before validation predictions were generated.

The representation-sensitivity analysis fixed RF importance, 50 genes, and LightGBM so that development-to-held-out behaviour could be compared directly across the three retained representations.

Confirmatory analyses

The balanced 16-observation held-out sets permit exact enumeration of:

choose(16, 8) = 12,870

possible case-label allocations.

Exact one-sided permutation tests were therefore used for held-out AUROC testing. Held-out representation comparisons used paired DeLong tests, and 5,000 bootstrap replicates quantified uncertainty in development-to-held-out AUROC differences.

Experiment 2: Disease-composition perturbation

Question

Does transcriptomic information learned across other diseases transfer to a previously unseen target disease, and how does target-disease exposure alter that relationship?

Common feature space

A shared cross-disease gene universe was defined from the five development datasets:

4,450 genes

This feature space was held fixed across the principal pooled-model comparison.

Pooled development populations

Model A
RA + UC + CD + T1D
n = 1,413

Model B
TB + RA + UC + CD + T1D
n = 1,736

Model B therefore added 323 TB development observations while predictor availability remained fixed.

Development pipeline

Representation: VST + fSVA
Feature selection: mRMR / CMIM
Class handling: none / training-only SMOTE
Learners: Random Forest, LightGBM, XGBoost, RBF SVM

Model optimisation used 5 outer folds, 10 inner folds, 20 model-based optimisation evaluations, and BACC as the tuning objective.

Leave-one-disease-out evaluation

LODO evaluation excluded each disease in turn from every data-dependent stage. The four remaining diseases formed the source distribution and the omitted disease formed the target distribution.

The locked LODO configuration was:

VST + fSVA | mRMR | 50 genes | no class balancing | LightGBM

The target disease entered only when probabilities were generated.

Frozen TB comparison

Model A and Model B were refitted using their complete development populations and evaluated on the same frozen TB test population.

Because Model B simultaneously introduced TB-specific information and increased development sample size, this comparison represents the joint effect of target-phenotype exposure and additional sample support.

Repeated permutation importance and Reactome pathway enrichment were used to place the cross-disease predictive structure in molecular context.

Experiment 3: Cohort-composition perturbation

Question

How does strong sample-level learnability compare with transportability to a completely omitted study population?

Development and frozen held-out split

For each disease, an approximately 80:20 development-to-held-out split was created using joint cohort-and-outcome stratification. Held-out sample identifiers were frozen in manifests before disease-specific optimisation.

Initial screening

The disease-specific screen considered:

Representations:
VST only
log2(CPM + 1)

Learners:
RBF SVM
Random Forest

RF permutation importance was used for training-specific gene ranking.

Definitive LOCO optimisation

The final analytical systems converged on Random Forest across all five diseases.

For each disease, development-only LOCO optimisation evaluated:

Feature counts:
25
50
100

Class weighting:
unweighted
class weighted

This produced six candidate configurations for the selected disease-specific representation.

Candidate pipelines were ranked according to:

1. Highest worst-cohort AUROC
2. Highest mean LOCO AUROC
3. Lowest between-cohort AUROC SD
4. Highest mean LOCO BACC

Operating threshold

Each final disease-specific classifier received a threshold derived only from development LOCO predictions.

Candidate thresholds ranged from 0.05 to 0.95 in increments of 0.01 and were ranked by:

1. Highest mean cohort BACC
2. Highest worst-cohort BACC
3. Closest value to 0.50

The selected threshold was frozen before final evaluation.

Final evaluation and molecular reliability

For every omitted target cohort, representation fitting, RF feature ranking, and classifier fitting used only the remaining source cohorts. The target cohort entered at prediction.

Evaluation included:

AUROC with DeLong 95% CI
AUPRC relative to target prevalence
BACC
Sensitivity
Specificity
Brier score
Calibration intercept and slope
Calibration curves

The locked disease-specific system was then reconstructed from the complete development population and applied to the frozen 20% held-out sample population.

Repeated delta-AUROC permutation importance was calculated 100 times per selected gene. First-order ALE profiles were generated for up to six leading genes to describe fitted expression-response functions.

Experiment 4: Model-architecture perturbation

Question

Does changing learner architecture materially alter predictive transportability or the molecular explanation supporting prediction?

Matched molecular input

Representation: log2(CPM + 1)
Feature selection: RF permutation importance
Selected genes: 50
Architectures: LightGBM and TabPFN

The two learners received identical training observations, validation observations, selected genes, and feature order.

Development and cohort-level comparison

LightGBM configurations were selected using development resampling. TabPFN used the pretrained default classifier.

Both architectures were evaluated using matched outer-CV predictions and the corresponding disease-specific LOCO design.

ROC, PR, and calibration analyses were generated from the same prediction-level outputs, with bootstrap uncertainty retained for the curves.

Frozen final evaluation

The final LightGBM and TabPFN models were fitted on the complete development population and evaluated on the same frozen diagnostic test observations.

Reported metrics included:

AUROC
AUPRC
BACC
Sensitivity
Specificity
Brier score

Model restoration

Final Python model objects were serialised after fitting. The experiment state also retained:

Disease identity
Representation
Selected 50 genes
Feature order
Selected LightGBM configuration
Class definitions
Development sample identifiers
Test sample identifiers
TabPFN execution device
Original prediction probabilities

Before molecular interpretation, restored models were required to reproduce the stored sample-level probabilities within the prespecified numerical tolerance.

Cross-architecture molecular reliability

Each of the 50 genes was permuted 100 times per architecture.

Gene-level concordance was evaluated using:

Spearman correlation of delta-AUROC importance ranks
Top-20 Jaccard overlap

ALE analysis then compared fitted expression-response functions across architectures.

Performance metrics

TRACE-Omics stores probability-level predictions and derives the principal metrics from those predictions.

Metric

Role

AUROC

Case-control ranking performance

AUPRC

Precision-recall performance interpreted relative to case prevalence

BACC

Threshold-dependent performance under class imbalance

Sensitivity

Case detection at the frozen operating threshold

Specificity

Control classification at the frozen operating threshold

Brier score

Probability prediction error

Calibration is assessed separately because ranking performance and probability scale can change differently across target populations.

Reproducibility controls

TRACE-Omics preserves the complete fit-and-apply relationship between development and evaluation data.

The workflow retains, where applicable:

Raw-count matrices
Harmonised metadata
Development and evaluation manifests
Fold assignments
Fitted preprocessing states
Feature rankings
Selected gene vectors
Class-handling records
Tuning archives
Prediction-level outputs
Fold-level performance summaries
Statistical-test outputs
Permutation distributions
ALE coordinates
Final fitted models
Session information

The global project seed is:

7832

Deterministic sub-seeds are generated for disease-, fold-, configuration-, gene-, and permutation-specific operations where required.

Evaluation observations enter only after the corresponding preprocessing, feature-selection, class-handling, tuning, and model-estimation states have been fitted from development or source observations.

Repository structure

trace-omics/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── config/
│   ├── example_config.yml
│   ├── hyperparameter_search_spaces.yml
│   └── lightgbm_architecture_grid.yml
├── data/
│   └── README.md
├── docs/
│   └── workflow_overview.md
├── R/
│   └── utils_io.R
├── scripts/
│   ├── 01_representation_specification/
│   ├── 02_feature_set_characterisation/
│   ├── 03_model_development_optimisation/
│   └── 04_generalisation_reliability/
└── results/

Exact executable scripts may be separated by disease or experiment inside these stage directories.

Recommended execution logic

1. Acquire public RNA-seq data
2. Curate eligible biological observations
3. Harmonise gene identifiers and metadata
4. Create and freeze development/evaluation manifests
5. Specify candidate transcriptomic representations
6. Characterise representation structure
7. Generate and characterise feature sets
8. Run the relevant controlled perturbation experiment
9. Optimise candidate pipelines using development-only evidence
10. Freeze the selected predictive system
11. Evaluate on the designated target population
12. Perform molecular reliability analysis on the locked system
13. Save predictions, model state, audit files, and session information

Large analyses should be run from a clean R session.

The LightGBM-TabPFN architecture experiment should initialise its specified Python environment before model fitting or restoration.

Computational environment

The dissertation workflow was coordinated in R.

RNA-seq preprocessing relied primarily on Bioconductor packages including:

DESeq2
sva
AnnotationDbi
org.Hs.eg.db

Model development used the mlr3 ecosystem together with learner-specific implementations. DALEX and ingredients were used for model-agnostic interpretation and ALE analysis.

The architecture experiment used R through reticulate with a pinned Python 3.10 environment containing:

NumPy 1.22.4
scikit-learn 1.2.2
LightGBM 4.6.0
PyTorch 2.4.1
TabPFN 0.1.11

Public repository guidance

Commit:

README.md
LICENSE
CITATION.cff
.gitignore
config/
docs/
R/
scripts/
data/README.md

Exclude generated or potentially large analytical data:

data/raw/
data/processed/
results/
*.rds
*.csv
*.png
*.tiff
*.pdf
*.html
*.pkl

The raw RNA-seq datasets remain available through their original public repositories and should be downloaded by users according to the accessions documented above.

Citation

If you use TRACE-Omics, please cite the repository together with the associated dissertation or subsequent manuscript when available.

If using TRACE-Omics, cite the repository together with the associated dissertation or subsequent manuscript when available.
