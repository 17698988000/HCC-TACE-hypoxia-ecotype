# Excluded uncertain files

This staging folder is intentionally conservative.

The script copies only the explicitly whitelisted final locked figures and final locked tables. It deliberately excludes all non-whitelisted files, including:

- Supplementary Figure S5 and any attempted S5 objects.
- Figure 2 v12/v13/v14 or any non-v11 Figure 2 object.
- Draft, candidate, failed, temporary, old, or backup files.
- ROC/AUC/cutoff/nomogram/calibration/decision-curve/KM high-low/TCGA-validation artifacts.
- Raw GEO expression matrices.
- Private or identifiable participant-level data.

Reason: the active manuscript is an endpoint-first translational association study, not a prediction model, classifier, treatment-selection tool, LASSO/KM risk model, or pan-HCC prognostic model.

<!-- A18B_ANALYSIS_CODE_STAGING_BEGIN -->

## A18b analysis-code staging exclusions

Generated: 2026-06-21 13:59:30 CST

Scope: only R/Python scripts from the explicitly allowed source directories were considered. Files with DRAFT/candidate/failed/temp/old/backup-like path terms, blocked Figure 2 versions, Supplementary Figure S5, prediction-model/performance terms, KM high-low, or TCGA-validation wording were not copied.

### Missing allowed source directories

- `a07_final_tables`

### Excluded R/Python scripts

_None._

<!-- A18B_ANALYSIS_CODE_STAGING_END -->

