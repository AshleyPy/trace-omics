#!/usr/bin/env Rscript
# TRACE-Omics entry point. The comprehensive implementation is kept in the
# stage subdirectory so later stages can import its function definitions.
source(file.path("scripts", "01_representation_specification", "01_representation_specification.R"), chdir = FALSE)
