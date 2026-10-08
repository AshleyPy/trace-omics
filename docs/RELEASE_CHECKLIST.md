# Release checklist

Before making the repository public:

- [ ] Select a software licence and replace `LICENSE_NOTICE.md` with the final licence text.
- [ ] Remove machine-specific paths, credentials, tokens and private participant crosswalks.
- [ ] Confirm every public study accession against Chapter 3.
- [ ] Restore or archive the frozen partition manifests used for the dissertation where sharing is permitted.
- [ ] Run `Rscript tests/parse_all.R` in the release environment.
- [ ] Run `Rscript tests/test_core.R`.
- [ ] Run the final analytical entry points against the author archive and compare prediction-level outputs.
- [ ] Record R and package versions from the successful release run.
- [ ] Verify the Python 3.10 / TabPFN 0.1.11 environment used for the architecture experiment.
- [ ] Confirm that generated outputs are excluded from Git unless intentionally released.
- [ ] Create a Git tag and archive the release with a persistent identifier if available.
