# The supplied code retains original cohort/outcome names and disease-specific paths.
# Review source_manifest.json and RUN_ORDER.md before selecting chunks.
list(workflow="equal_size",chunks=1L,
 output="results/prepared_equal_size.R",
 replacements=c(
   'file.path(user_home, "Downloads", "RawIBDdata")'='"/absolute/path/to/RawIBDdata"',
   '"CD_EQUAL_SIZE_41v41_CLEAN_BENCHMARK"'='"CD_EQUAL_SIZE_REPLAY_001"'
 ))
