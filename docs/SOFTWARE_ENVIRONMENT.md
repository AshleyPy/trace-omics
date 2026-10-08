# Software and computing environment

The dissertation workflow was coordinated in R with the mlr3 ecosystem. The recorded project seed was 7832 and the R random-number generator specification was Mersenne–Twister / Inversion / Rejection.

Principal packages represented in the surviving source include DESeq2, sva, edgeR, neuroCombat, mlr3, mlr3filters, mlr3learners, mlr3extralearners, mlr3mbo, paradox, bbotk, ranger, lightgbm, xgboost, e1071, C50, pROC, PRROC, DALEX and ingredients.

The matched architecture workflow used Python through reticulate. The recorded environment included Python 3.10, NumPy 1.22.4, scikit-learn 1.2.2, LightGBM 4.6.0, PyTorch 2.4.1 and TabPFN 0.1.11.

Interactive development, debugging, molecular-data processing and intermediate validation were performed on the in-house Gaia server with 512 GB RAM, 28 TB RAID 60 storage and dual Intel Xeon Gold processors.

A generated dependency lock should be created from the successfully verified release environment instead of inferring package versions from source code alone.
