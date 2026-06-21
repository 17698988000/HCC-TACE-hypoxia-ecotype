# a18b analysis-code staging QA report

Generated: 2026-06-21 13:59:30 CST

Final verdict: **PASS**

## Scope

This a18b step adds analysis-related R/Python scripts into the existing local GitHub staging folder. It does not rerun analyses, alter locked results, modify final figures/tables, or upload to GitHub.

Staging folder: `D:/Documents/肿瘤—内皮—巨噬细胞生态位的论文书写/results/a18_github_repository_staging`
Target repository URL recorded for documentation only: `https://github.com/17698988000/HCC-TACE-hypoxia-ecotype`

## QA checks

| check | observation | result |
| --- | --- | --- |
| Existing a18 staging folder exists | D:/Documents/肿瘤—内皮—巨噬细胞生态位的论文书写/results/a18_github_repository_staging | PASS |
| Allowed source directories searched only | a00_1_response_definition_rescue; a00_pre_analysis_risk_control; a01_data_preparation; a02_ecotype_construction; a03_gse104580_response_association; a04_gse14520_treatment_context_validation; a05_secondary_externalization; a07_final_tables; a08_publication_figures; a09_submission_ready_assets; a09_visual_assets_final; a12_final_locked_tables; a16_supplementary_table_S5_S6_defense_audit; a17_supplementary_table_S5_S6_publication_ready_fixed_v3_S4_style | PASS |
| R/Python candidate scripts found | 10 candidate R/Python scripts found; 10 copied; 0 excluded/uncertain | PASS |
| Scripts folder is not only create_github_repository_staging_a18_v2.R | 11 R/Python scripts currently under staging/scripts | PASS |
| Data preparation script(s) staged | source_dir match: TRUE | PASS |
| Ecotype construction script(s) staged | source_dir/name match: TRUE | PASS |
| GSE104580 response association script(s) staged | source_dir/name match: TRUE | PASS |
| GSE14520 outcome/context script(s) staged | source_dir/name match: TRUE | PASS |
| GSE76427 externalization script(s) staged | source_dir/name match: TRUE | PASS |
| Final table/figure generation or formatting script(s) staged | source_dir match: TRUE | PASS |
| No forbidden filename/path among staged scripts | 0 forbidden staged filenames | PASS |
| outputs/final_figures unchanged | 8 files before; 8 files after | PASS |
| outputs/final_tables unchanged | 9 files before; 9 files after | PASS |
| GitHub upload not performed | No git command, GitHub API call, or upload routine is used by this script | PASS |

## Staged script count by destination subfolder

| Var1 | Freq |
| --- | --- |
| 01_data_preparation | 2 |
| 02_ecotype_construction | 1 |
| 04_association_analyses | 5 |
| 06_final_table_formatting | 2 |

## Staged scripts

| destination_subdir | destination_file | source_relative_path | source_dir | size_bytes | md5 |
| --- | --- | --- | --- | --- | --- |
| 01_data_preparation | 01_data_preparation/a01_data_preparation__script_A02_0_entry_gate.py | a01_data_preparation/script_A02_0_entry_gate.py | a01_data_preparation | 17525 | 7ef94d66732c82beef83795420d1b2b0 |
| 01_data_preparation | 01_data_preparation/a01_data_preparation__a01_prepare_data.py | a01_data_preparation/scripts/a01_prepare_data.py | a01_data_preparation | 27420 | 99d9b307b3a21679fd85bfedf53b9331 |
| 02_ecotype_construction | 02_ecotype_construction/a02_ecotype_construction__script_A02_construct_ecotype.py | a02_ecotype_construction/script_A02_construct_ecotype.py | a02_ecotype_construction | 22500 | 19a0884003da08d3fcdc71dda12db989 |
| 04_association_analyses | 04_association_analyses/a03_gse104580_response_association__script_A03_gse104580_response_association.py | a03_gse104580_response_association/script_A03_gse104580_response_association.py | a03_gse104580_response_association | 19568 | 628576572a72bf349ff9019cc4ea50a6 |
| 04_association_analyses | 04_association_analyses/a04_gse14520_treatment_context_validation__script_A04_gse14520_treatment_context_validation.py | a04_gse14520_treatment_context_validation/script_A04_gse14520_treatment_context_validation.py | a04_gse14520_treatment_context_validation | 17244 | 667e1eb7fb4376cd443d9f6db90764cb |
| 04_association_analyses | 04_association_analyses/a05_secondary_externalization__script_A05_fast_rescue_gse76427.py | a05_secondary_externalization/script_A05_fast_rescue_gse76427.py | a05_secondary_externalization | 16065 | fe1696ddb2e84153b1d5a4347195eaf2 |
| 04_association_analyses | 04_association_analyses/a05_secondary_externalization__script_A05_rescue_mode.py | a05_secondary_externalization/script_A05_rescue_mode.py | a05_secondary_externalization | 17079 | 50f22eac68df506d121913a0131528c1 |
| 04_association_analyses | 04_association_analyses/a05_secondary_externalization__script_A05_secondary_externalization.py | a05_secondary_externalization/script_A05_secondary_externalization.py | a05_secondary_externalization | 22811 | 650536103da8a381c24591ba7e5d7c17 |
| 06_final_table_formatting | 06_final_table_formatting/a09_submission_ready_assets__generate_a09_tables_and_inputs.py | a09_submission_ready_assets/scripts/generate_a09_tables_and_inputs.py | a09_submission_ready_assets | 24687 | 05e7a9cfddb253d675ae367465f13dc1 |
| 06_final_table_formatting | 06_final_table_formatting/a09_visual_assets_final__generate_a09_tables_and_inputs.py | a09_visual_assets_final/scripts/generate_a09_tables_and_inputs.py | a09_visual_assets_final | 24687 | 05e7a9cfddb253d675ae367465f13dc1 |

## Excluded or uncertain scripts

See `docs/excluded_uncertain_files.md` for detailed exclusion records.

## Output integrity boundary

- `outputs/final_figures`: unchanged
- `outputs/final_tables`: unchanged

## GitHub upload status

No upload was performed. This script only stages local files.

