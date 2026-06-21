# HCC-TACE-hypoxia-ecotype

Target repository: https://github.com/17698988000/HCC-TACE-hypoxia-ecotype

This local staging folder supports the manuscript:

> A single-cell and spatially anchored hypoxia-adaptive tumor-endothelial-macrophage ecotype associated with mRECIST-defined TACE non-response and adverse post-TACE outcome in hepatocellular carcinoma.

## Repository scope

This staging folder contains only final locked manuscript-supporting figures, final locked tables, reproducibility notes, session information, and a file manifest prepared for repository upload.

Raw GEO matrices are not redistributed. Users should download public GEO datasets directly from GEO and then run verified local workflows against those public inputs.

## Interpretation boundaries

- The fixed ecotype was not endpoint-optimized.
- GSE14520/Fako is not a response-validation cohort.
- S5/S6 are audit/feasibility tables, not validation or adjusted inference.
- This is not a prediction model or treatment-selection tool.
- The repository is intended for transparent reproduction of reported association analyses from public datasets, not for clinical deployment.

## Included final figures

- `outputs/final_figures/Figure_1_endpoint_first_workflow_v14.pdf`
- `outputs/final_figures/Figure_2_ecotype_construction_v11.pdf`
- `outputs/final_figures/Figure_3_gse104580_response_association_v11.pdf`
- `outputs/final_figures/Figure_4_gse14520_treatment_context_forest_v11.pdf`
- `outputs/final_figures/Supplementary_Figure_S1_gse76427_rescue_limited_externalization_v6.pdf`
- `outputs/final_figures/Supplementary_Figure_S2_gse104580_component_level_behavior_v4.pdf`
- `outputs/final_figures/Supplementary_Figure_S3_gene_transferability_duplicate_dll4_v4.pdf`
- `outputs/final_figures/Supplementary_Figure_S4_gse14520_individual_followup_map_final_WPS_compatible.pdf`

## Included final tables

- `outputs/final_tables/Table_1_clinical_bulk_cohorts_and_prespecified_endpoints_LOCKED.xlsx`
- `outputs/final_tables/Table_2_fixed_ecotype_components_and_transfer_rules_LOCKED.xlsx`
- `outputs/final_tables/Table_3_locked_association_estimates_across_cohorts_LOCKED.xlsx`
- `outputs/final_tables/Supplementary_Table_S1_full_component_provenance_and_transfer_rules_LOCKED.xlsx`
- `outputs/final_tables/Supplementary_Table_S2_cohort_manifest_endpoint_completeness_LOCKED.xlsx`
- `outputs/final_tables/Supplementary_Table_S3_score_coverage_descriptive_summaries_LOCKED.xlsx`
- `outputs/final_tables/Supplementary_Table_S4_expanded_association_and_boundary_outputs_LOCKED.xlsx`
- `outputs/final_tables/Supplementary_Table_S5_GSE104580_sample_level_endpoint_harmonization_inclusion_audit_LOCKED.xlsx`
- `outputs/final_tables/Supplementary_Table_S6_clinical_covariate_availability_adjustment_feasibility_audit_LOCKED.xlsx`

## Not redistributed

- Large raw expression matrices.
- Private data.
- Identifiable participant-level information.
- Draft, candidate, failed, temporary, old, or backup figure/table objects.
- Supplementary Figure S5.

## Notes for users

The `scripts/` directory is created for repository organization. This staging script does not auto-copy non-whitelisted historical analysis scripts because older project folders may contain draft, failed, or version-conflicting objects. Add only verified final scripts after a separate code audit.

## Script organization and provenance notes

The staged scripts are organized by analysis stage. Score transfer is implemented within the cohort-specific association scripts where applicable rather than as a separate standalone `03_score_transfer` directory. The fixed ecotype object is constructed before endpoint testing and is transferred without response- or outcome-guided gene revision.

Two staged copies of `generate_a09_tables_and_inputs.py` may be present because identical copies existed in the final table/visual-asset provenance folders. They are retained as provenance copies for auditability; the files have identical content in the a18b staging manifest.

## Repository interpretation and script notes

The staged scripts are organized by analysis stage. Score transfer is implemented within the cohort-specific association scripts where applicable rather than as a separate standalone `03_score_transfer` directory. The fixed ecotype object is constructed before endpoint testing and is transferred without response- or outcome-guided gene revision.

Two staged copies of `generate_a09_tables_and_inputs.py` may be present because identical copies existed in the final table/visual-asset provenance folders. They are retained as provenance copies for auditability; the files have identical content in the a18b staging manifest.

This repository supports transparency and reproducibility of the reported computational workflow. It does not redistribute large raw GEO expression matrices. Users should download public datasets directly from GEO and run the documented scripts with the project-specific inputs described in the manuscript and supplementary files.

The repository is not intended to define a clinical prediction model, response classifier, treatment-selection tool, ROC/AUC workflow, cut-point optimization procedure, nomogram, decision-curve analysis, Kaplan-Meier high/low score grouping, or TCGA validation analysis.

GSE14520/Fako is treated as a treatment-context outcome cohort, not as a radiological response-validation cohort. Supplementary Tables S5 and S6 are endpoint-harmonization and covariate-feasibility audit tables, not validation analyses or covariate-adjusted inference.

