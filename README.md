# HCC-TACE-hypoxia-ecotype

Analysis code and reproducibility materials for a single-cell and spatially anchored hypoxia-adaptive tumor-endothelial-macrophage ecotype associated with mRECIST-defined TACE non-response and adverse post-TACE outcome in hepatocellular carcinoma.

## Repository scope

This repository supports transparent reproduction of the reported computational workflow, including fixed ecotype construction, score transfer, cohort-level association analyses, and final figure/table generation records.

The repository contains:

* analysis scripts;
* reproducibility notes;
* session information;
* final locked manuscript-supporting figures;
* final locked tables;
* file manifests and provenance records.

Raw GEO expression matrices are not redistributed. Users should download public datasets directly from GEO and run the documented scripts with the project-specific inputs described in the manuscript and supplementary files.

## Public datasets

The analyses use publicly available, de-identified datasets from the NCBI Gene Expression Omnibus, including:

* GSE149614
* GSE151530
* GSE277104
* GSE104580
* GSE14520
* GSE76427

## Repository structure

```text
docs/
  reproducibility_notes.md

metadata/
  FINAL_ANALYSIS_SCRIPT_MANIFEST.md
  FINAL_GITHUB_STAGING_FILE_MANIFEST.md

outputs/
  final_figures/
  final_tables/

scripts/
  01_data_preparation/
  02_ecotype_construction/
  04_association_analyses/
  06_final_table_formatting/

session_info/
  R_session_info.txt
```

Score transfer is implemented within the cohort-specific association scripts where applicable rather than as a separate standalone directory. The fixed ecotype object was constructed before endpoint testing and transferred without response- or outcome-guided gene revision.

## Included final figures

* `outputs/final_figures/Figure_1_endpoint_first_workflow_v14.pdf`
* `outputs/final_figures/Figure_2_ecotype_construction_v11.pdf`
* `outputs/final_figures/Figure_3_gse104580_response_association_v11.pdf`
* `outputs/final_figures/Figure_4_gse14520_treatment_context_forest_v11.pdf`
* `outputs/final_figures/Supplementary_Figure_S1_gse76427_rescue_limited_externalization_v6.pdf`
* `outputs/final_figures/Supplementary_Figure_S2_gse104580_component_level_behavior_v4.pdf`
* `outputs/final_figures/Supplementary_Figure_S3_gene_transferability_duplicate_dll4_v4.pdf`
* `outputs/final_figures/Supplementary_Figure_S4_gse14520_individual_followup_map_final_WPS_compatible.pdf`

## Included final tables

* `outputs/final_tables/Table_1_clinical_bulk_cohorts_and_prespecified_endpoints_LOCKED.xlsx`
* `outputs/final_tables/Table_2_fixed_ecotype_components_and_transfer_rules_LOCKED.xlsx`
* `outputs/final_tables/Table_3_locked_association_estimates_across_cohorts_LOCKED.xlsx`
* `outputs/final_tables/Supplementary_Table_S1_full_component_provenance_and_transfer_rules_LOCKED.xlsx`
* `outputs/final_tables/Supplementary_Table_S2_cohort_manifest_endpoint_completeness_LOCKED.xlsx`
* `outputs/final_tables/Supplementary_Table_S3_score_coverage_descriptive_summaries_LOCKED.xlsx`
* `outputs/final_tables/Supplementary_Table_S4_expanded_association_and_boundary_outputs_LOCKED.xlsx`
* `outputs/final_tables/Supplementary_Table_S5_GSE104580_sample_level_endpoint_harmonization_inclusion_audit_LOCKED.xlsx`
* `outputs/final_tables/Supplementary_Table_S6_clinical_covariate_availability_adjustment_feasibility_audit_LOCKED.xlsx`

## Interpretation boundaries

The fixed ecotype was not endpoint-optimized.

GSE14520/Fako was used as a treatment-context outcome cohort, not as a radiological response-validation cohort.

Supplementary Tables S5 and S6 are endpoint-harmonization and covariate-feasibility audit tables, not validation analyses or covariate-adjusted inference.

This repository is not intended to define a clinical prediction model, response classifier, treatment-selection tool, ROC/AUC workflow, cut-point optimization procedure, nomogram, decision-curve analysis, Kaplan-Meier high/low score grouping, or TCGA validation analysis.

## Not redistributed

This repository does not redistribute:

* large raw expression matrices;
* private data;
* identifiable participant-level information;
* draft, candidate, failed, temporary, old, or backup figure/table objects;
* Supplementary Figure S5.

## License

This repository is released under the MIT License.
