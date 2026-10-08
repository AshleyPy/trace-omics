# Reproduction status

TRACE-Omics separates **source-code availability** from **numerical reproduction of the dissertation**. The repository contains the surviving executable code, reconstructed stage drivers, utility functions and legacy analysis branches. Exact numerical reproduction additionally requires the analytical state that defined the original observations and fitted systems.

## Inputs required for numerical reproduction

- harmonised gene-level count matrices after the dissertation eligibility rules;
- participant and sample crosswalks, including repeated-measure resolution;
- the frozen 70:30 and 80:20 partition manifests;
- equal-size repetition sample identifiers where archived;
- original resampling maps and tuning archives for analyses selected through model-based optimisation;
- the fixed 4,450-gene pooled feature universe in its original order;
- fitted preprocessing states, selected-gene vectors and development-derived thresholds where replay rather than refitting is intended;
- the recorded TabPFN 0.1.11 environment and checkpoint for the matched-architecture experiment.

## Final-design boundary

The final dissertation equal-size sensitivity analysis used **100 repeated study-aware samples of n=82 (41 cases and 41 controls)** for TB, RA, UC and CD. Earlier code branches that split an n=82 sample into a 66-observation development set and 16-observation final evaluation set are retained only under `legacy/` or recovered stage code as provenance from framework development. They must not be used to represent the final equal-size sensitivity estimand.

The final LODO-CV analysis carried the prespecified primary pooled system—VST + fSVA, mRMR with 50 genes, observed class distribution and LightGBM—across each complete disease exclusion. This choice holds analytical specification fixed while changing the disease-population boundary.

## Verification before a tagged release

Run the parsing and synthetic tests, execute each retained analysis in the intended R/Python environments, compare participant IDs and predictor order with the dissertation archive, and compare generated prediction probabilities and tables with the saved Chapter 4 sources. Tag a public release only after those checks have been completed against the author-side analytical state.
