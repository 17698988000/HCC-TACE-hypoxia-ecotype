# Reproducibility notes

Target repository: https://github.com/17698988000/HCC-TACE-hypoxia-ecotype

## Data redistribution boundary

Raw GEO matrices are not redistributed. Users should download public GEO datasets directly from GEO accessions cited in the manuscript and then apply the documented local workflows.

## Fixed ecotype boundary

The fixed ecotype was not endpoint-optimized. Response labels, outcome data, or disease-context externalization results were not used to add, remove, reweight, replace, or optimize ecotype genes.

## Cohort-role boundary

GSE104580 is the pretreatment mRECIST-defined response-association cohort. GSE14520/Fako is not a response-validation cohort; it is used for treatment-context post-TACE outcome association and a resection-only specificity-boundary comparison.

## Supplementary Table S5/S6 boundary

S5/S6 are audit/feasibility tables, not validation or adjusted inference. Supplementary Table S5 documents GSE104580 sample-level endpoint harmonization and inclusion. Supplementary Table S6 documents clinical covariate availability and adjustment-feasibility boundaries.

## Clinical interpretation boundary

This is not a prediction model or treatment-selection tool. The score is retained as a continuous biological association measure and should not be interpreted as a cutoff, classifier, nomogram, calibration framework, treatment-allocation rule, or clinical decision-support system.

## Files deliberately excluded from this staging folder

Non-whitelisted figures/tables, draft candidates, failed attempts, temporary objects, old backups, Supplementary Figure S5 objects, raw matrices, and performance-model artifacts are excluded.

## a18c staging clarification

The staged scripts are organized by analysis stage. Score transfer is implemented within the cohort-specific association scripts where applicable rather than as a separate standalone `03_score_transfer` directory. The fixed ecotype object is constructed before endpoint testing and is transferred without response- or outcome-guided gene revision.

Two staged copies of `generate_a09_tables_and_inputs.py` may be present because identical copies existed in the final table/visual-asset provenance folders. They are retained as provenance copies for auditability; the files have identical content in the a18b staging manifest. 

 This repository supports transparency and reproducibility of the reported computational workflow. It does not redistribute large raw GEO expression matrices. Users should download public datasets directly from GEO and run the documented scripts with the project-specific inputs described in the manuscript and supplementary files.

The repository is not intended to define a clinical prediction model, response classifier, treatment-selection tool, ROC/AUC workflow, cut-point optimization procedure, nomogram, decision-curve analysis, Kaplan-Meier high/low score grouping, or TCGA validation analysis.

GSE14520/Fako is treated as a treatment-context outcome cohort, not as a radiological response-validation cohort. Supplementary Tables S5 and S6 are endpoint-harmonization and covariate-feasibility audit tables, not validation analyses or covariate-adjusted inference.

## Repository interpretation and script notes

The staged scripts are organized by analysis stage. Score transfer is implemented within the cohort-specific association scripts where applicable rather than as a separate standalone `03_score_transfer` directory. The fixed ecotype object is constructed before endpoint testing and is transferred without response- or outcome-guided gene revision.

Two staged copies of `generate_a09_tables_and_inputs.py` may be present because identical copies existed in the final table/visual-asset provenance folders. They are retained as provenance copies for auditability; the files have identical content in the a18b staging manifest.

This repository supports transparency and reproducibility of the reported computational workflow. It does not redistribute large raw GEO expression matrices. Users should download public datasets directly from GEO and run the documented scripts with the project-specific inputs described in the manuscript and supplementary files.

The repository is not intended to define a clinical prediction model, response classifier, treatment-selection tool, ROC/AUC workflow, cut-point optimization procedure, nomogram, decision-curve analysis, Kaplan-Meier high/low score grouping, or TCGA validation analysis.

GSE14520/Fako is treated as a treatment-context outcome cohort, not as a radiological response-validation cohort. Supplementary Tables S5 and S6 are endpoint-harmonization and covariate-feasibility audit tables, not validation analyses or covariate-adjusted inference.

