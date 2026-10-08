# Run order and source mapping

## Portable scripts

1. Check the environment and run the synthetic tests.
2. Prepare harmonised inputs and a frozen manifest; configure one disease per run for disease-specific analyses.
3. Run stages 01 and 02 for development-only representation and feature-set characterisation.
4. Run stage 03 for the supplied fixed learner benchmark, or select the clearly labelled reconstructed nested-grid mode with an explicit candidate list.
5. Run stage 05 independently for RF study exclusion and its frozen holdout prediction. Select the RF ranking rule explicitly.
6. Run stage 08 independently for the matched architecture experiment. It uses its own inner tuning and source-specific feature selection.
7. Run stage 04 on one set of saved holdout probabilities and its previously selected threshold.
8. Stage 06 constructs a new equal-size manifest. Configure a fresh run on its outputs; recover the original IDs for a historical comparison.
9. Stage 07 audits the pooled population and disease exclusions. It produces manifests, not the missing principal pooled model.

## Supplied source collections

`docs/source_manifest.json` maps every R file to its exact source text and line range.

- `legacy/equal_size/01.R`: supplied CD equal-size workflow; includes representations, FS, fixed learner comparison and later evaluation/plotting blocks.
- `legacy/study_models/01.R`–`29.R`: supplied full-model notebook. Several headings retain TB/T1D wording while paths identify CD inputs. Later chunks contain alternative architecture workflows. These are evidence of the surviving source, not proof that each branch produced a retained dissertation table. Inspect each chunk before replay.
- `legacy/matched_architectures/01.R`–`05.R`: supplied matched-model workflow and interpretation/plotting extensions. Later chunks can depend on objects from earlier chunks and can duplicate or replace earlier outputs.
- `legacy/pooled_combatref/01.R`–`02.R`: supplied alternative pooled VST + ComBat-Ref/Boruta/LightGBM experiment and selected-gene stability summary. The subsequent supplied cross-disease architecture proposal was excluded because it was outside the retained Results.

Use `scripts/90_prepare_legacy.R` to combine selected chunks and apply explicitly configured literal replacements. It writes a script for review and parses it, then stops. It does not execute the script. Preserve original chunks unchanged. Adapted scripts should have their own configuration and output directory.

The raw workflows may clear their interactive session, install packages, create Python environments, download a checkpoint, overwrite their own output files or require Windows desktop paths. Run adapted copies in a separate R process after reviewing their input identities and output paths. A global replacement of disease names cannot establish a valid disease-specific input manifest.
