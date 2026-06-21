# HCC-TACE-hypoxia-ecotype

Analysis code for a single-cell and spatially anchored hypoxia-adaptive tumor-endothelial-macrophage ecotype associated with mRECIST-defined TACE non-response and adverse post-TACE outcome in hepatocellular carcinoma.

## Project overview

This repository contains analysis scripts and reproducibility materials for an endpoint-first translational bioinformatics study of hepatocellular carcinoma treated with transarterial chemoembolization (TACE).

The study evaluates a fixed pretreatment hypoxia-adaptive tumor-endothelial-macrophage ecotype across predefined public cohort roles, including:

- biological anchoring using single-cell and spatial transcriptomic resources;
- pretreatment mRECIST-defined TACE response association in GSE104580;
- treatment-context post-TACE outcome association and resection-only specificity-boundary assessment in GSE14520/Fako;
- rescue-limited secondary externalization in GSE76427.

The repository is intended to support transparent reproduction of the reported score construction, platform mapping, cohort processing, association analyses, and figure/table generation workflows.

## Public datasets

The analyses use publicly available, de-identified datasets from the NCBI Gene Expression Omnibus (GEO), including:

- GSE149614
- GSE151530
- GSE277104
- GSE104580
- GSE14520
- GSE76427

Raw expression matrices are not redistributed in this repository. Users should download public datasets directly from GEO and process them using the documented scripts.

## Repository contents

Planned structure:

```text
scripts/
  01_data_preparation/
  02_ecotype_construction/
  03_score_transfer/
  04_association_analyses/
  05_tables_and_figures/

metadata/
  cohort_manifests/
  platform_mapping_records/
  endpoint_harmonization_audits/

outputs/
  README_outputs.md

session_info/
  R_session_info.txt
