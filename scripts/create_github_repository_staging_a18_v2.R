# create_github_repository_staging_a18_v2.R
# Purpose:
#   Create a clean local GitHub repository staging folder for the HCC TACE ecotype manuscript.
#   This script does NOT upload to GitHub and does NOT rerun any biological/statistical analysis.
#
# Why this v2 exists:
#   v1 was too strict when the same whitelisted basename existed in more than one old project folder.
#   This v2 still copies ONLY the exact whitelisted final filenames, but if multiple exact-basename
#   copies are found it selects one deterministically by preferred final-asset folders and logs all
#   candidates in the manifest/QA report. Missing whitelisted files still cause a clear FAIL.
#
# Run in RStudio:
#   1) Save this file locally.
#   2) Open it in RStudio.
#   3) Click Source.
#
# Important behavior:
#   - Existing a18_github_repository_staging will be removed and recreated.
#   - Only whitelisted final locked figures and tables are copied.
#   - Table 1 is corrected only in the staging copy if "not analysed" is found.
#   - Raw GEO matrices are not copied or redistributed.
#   - No GitHub upload is performed.

options(stringsAsFactors = FALSE)

# -----------------------------
# 0. User-fixed paths and config
# -----------------------------

project_root <- "D:/Documents/肿瘤—内皮—巨噬细胞生态位的论文书写/results"
staging_dir_name <- "a18_github_repository_staging"
repo_url <- "https://github.com/17698988000/HCC-TACE-hypoxia-ecotype"

figure_whitelist <- c(
  "Figure_1_endpoint_first_workflow_v14.pdf",
  "Figure_2_ecotype_construction_v11.pdf",
  "Figure_3_gse104580_response_association_v11.pdf",
  "Figure_4_gse14520_treatment_context_forest_v11.pdf",
  "Supplementary_Figure_S1_gse76427_rescue_limited_externalization_v6.pdf",
  "Supplementary_Figure_S2_gse104580_component_level_behavior_v4.pdf",
  "Supplementary_Figure_S3_gene_transferability_duplicate_dll4_v4.pdf",
  "Supplementary_Figure_S4_gse14520_individual_followup_map_final_WPS_compatible.pdf"
)

table_whitelist <- c(
  "Table_1_clinical_bulk_cohorts_and_prespecified_endpoints_LOCKED.xlsx",
  "Table_2_fixed_ecotype_components_and_transfer_rules_LOCKED.xlsx",
  "Table_3_locked_association_estimates_across_cohorts_LOCKED.xlsx",
  "Supplementary_Table_S1_full_component_provenance_and_transfer_rules_LOCKED.xlsx",
  "Supplementary_Table_S2_cohort_manifest_endpoint_completeness_LOCKED.xlsx",
  "Supplementary_Table_S3_score_coverage_descriptive_summaries_LOCKED.xlsx",
  "Supplementary_Table_S4_expanded_association_and_boundary_outputs_LOCKED.xlsx",
  "Supplementary_Table_S5_GSE104580_sample_level_endpoint_harmonization_inclusion_audit_LOCKED.xlsx",
  "Supplementary_Table_S6_clinical_covariate_availability_adjustment_feasibility_audit_LOCKED.xlsx"
)

# Preferred folders are used ONLY when the exact same whitelisted filename appears in multiple places.
# This does not infer or change figure/table versions; the basename whitelist remains the controlling rule.
preferred_source_patterns <- c(
  "/a09_submission_ready_assets/",
  "/a09_visual_assets_final/",
  "/a12_final_locked_tables/",
  "/a17_supplementary_table_s5_s6_publication_ready_fixed_v2/",
  "/a17_supplementary_table_s5_s6_publication_ready_fixed_v3_s4_style/",
  "/a17_supplementary_table_s5_s6_publication_ready/",
  "/a12_submission_tables/",
  "/a12_submission_tables_rebuilt/",
  "/a08_publication_figures/",
  "/a07_tables_final/"
)

# -----------------------------
# 1. Helper functions
# -----------------------------

normalize_slash <- function(x, mustWork = FALSE) {
  normalizePath(x, winslash = "/", mustWork = mustWork)
}

write_utf8 <- function(lines, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  con <- file(path, open = "w", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(enc2utf8(lines), con = con, useBytes = TRUE)
}

md_escape <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("\\|", "\\\\|", x)
  x <- gsub("\r?\n", "<br>", x)
  x
}

md_table <- function(df) {
  if (nrow(df) == 0) {
    return(c("| Check | Observation | Result |", "|---|---|---|"))
  }
  header <- paste0("| ", paste(md_escape(names(df)), collapse = " | "), " |")
  sep <- paste0("| ", paste(rep("---", ncol(df)), collapse = " | "), " |")
  rows <- apply(df, 1, function(z) paste0("| ", paste(md_escape(z), collapse = " | "), " |"))
  c(header, sep, rows)
}

qa <- data.frame(
  Check = character(),
  Observation = character(),
  Result = character(),
  stringsAsFactors = FALSE
)

add_qa <- function(check, observation, pass) {
  qa[nrow(qa) + 1, ] <<- list(
    Check = check,
    Observation = as.character(observation),
    Result = if (isTRUE(pass)) "PASS" else "FAIL"
  )
}

write_qa_report <- function(staging_dir, final_note = NULL) {
  qa_path <- file.path(staging_dir, "a18_github_repository_staging_QA_report.md")
  final_status <- if (nrow(qa) > 0 && all(qa$Result == "PASS")) "PASS" else "FAIL"
  lines <- c(
    "# a18 GitHub repository staging QA report",
    "",
    paste0("Generated: `", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "`"),
    "",
    paste0("Project root: `", project_root, "`"),
    paste0("Staging folder: `", file.path(project_root, staging_dir_name), "`"),
    paste0("Target GitHub repository: `", repo_url, "`"),
    "",
    paste0("Final status: **", final_status, "**"),
    ""
  )
  if (!is.null(final_note)) {
    lines <- c(lines, "## Note", "", final_note, "")
  }
  lines <- c(lines, "## QA checks", "", md_table(qa), "")
  write_utf8(lines, qa_path)
  invisible(qa_path)
}

fail_and_stop <- function(staging_dir, message) {
  add_qa("Fatal staging stop", message, FALSE)
  write_qa_report(staging_dir, final_note = "Staging stopped before completion. Fix the failed item(s), then rerun the script.")
  stop(message, call. = FALSE)
}

ensure_openxlsx <- function(staging_dir) {
  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    message("R package 'openxlsx' not found. Attempting to install from CRAN...")
    install_ok <- tryCatch({
      utils::install.packages("openxlsx", repos = "https://cloud.r-project.org")
      TRUE
    }, error = function(e) {
      message("Failed to install openxlsx: ", conditionMessage(e))
      FALSE
    })
    if (!isTRUE(install_ok) || !requireNamespace("openxlsx", quietly = TRUE)) {
      fail_and_stop(
        staging_dir,
        "The R package 'openxlsx' is required to patch/check XLSX staging copies and could not be installed. Install openxlsx manually, then rerun."
      )
    }
  }
  invisible(TRUE)
}

get_this_script_path <- function() {
  cmd <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd, value = TRUE)
  if (length(file_arg) > 0) {
    p <- sub("^--file=", "", file_arg[1])
    if (file.exists(p)) return(normalize_slash(p, mustWork = TRUE))
  }

  frames <- sys.frames()
  for (i in rev(seq_along(frames))) {
    ofile <- frames[[i]]$ofile
    if (!is.null(ofile) && is.character(ofile) && length(ofile) == 1 && file.exists(ofile)) {
      return(normalize_slash(ofile, mustWork = TRUE))
    }
  }

  return(NA_character_)
}

source_priority <- function(path) {
  p <- tolower(normalize_slash(path, mustWork = FALSE))
  pats <- tolower(preferred_source_patterns)
  hit <- which(vapply(pats, function(z) grepl(z, p, fixed = TRUE), logical(1)))
  if (length(hit) == 0) return(length(pats) + 100L)
  min(hit)
}

select_best_candidate <- function(matches) {
  matches <- normalize_slash(matches, mustWork = TRUE)
  info <- file.info(matches)
  md5s <- as.character(tools::md5sum(matches))
  priorities <- vapply(matches, source_priority, integer(1))
  # Deterministic rule: preferred folder first, then newest mtime, then shortest path, then alphabetical path.
  ord <- order(priorities, -as.numeric(info$mtime), nchar(matches), matches)
  selected <- matches[ord[1]]

  data.frame(
    candidate_path = matches,
    basename = basename(matches),
    source_priority = priorities,
    mtime = format(info$mtime, "%Y-%m-%d %H:%M:%S"),
    size_bytes = as.integer(info$size),
    md5 = md5s,
    selected = matches == selected,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

find_whitelist_files <- function(root_dir, whitelist, staging_dir) {
  root_norm <- normalize_slash(root_dir, mustWork = TRUE)
  staging_norm <- normalize_slash(staging_dir, mustWork = FALSE)

  all_files <- list.files(
    root_norm,
    recursive = TRUE,
    full.names = TRUE,
    all.files = TRUE,
    no.. = TRUE,
    include.dirs = FALSE
  )
  all_files <- normalize_slash(all_files, mustWork = FALSE)
  all_files <- all_files[file.exists(all_files)]
  all_files <- all_files[file.info(all_files)$isdir == FALSE]
  all_files <- all_files[!startsWith(all_files, paste0(staging_norm, "/"))]

  found <- setNames(vector("list", length(whitelist)), whitelist)
  missing <- character()
  candidate_rows <- list()

  for (nm in whitelist) {
    matches <- all_files[tolower(basename(all_files)) == tolower(nm)]
    matches <- matches[file.info(matches)$isdir == FALSE]

    if (length(matches) == 0) {
      missing <- c(missing, nm)
      next
    }

    cand <- select_best_candidate(matches)
    candidate_rows[[length(candidate_rows) + 1]] <- data.frame(
      whitelist_name = nm,
      cand,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    found[[nm]] <- cand$candidate_path[cand$selected][1]
  }

  candidate_df <- if (length(candidate_rows) == 0) {
    data.frame(
      whitelist_name = character(), candidate_path = character(), basename = character(),
      source_priority = integer(), mtime = character(), size_bytes = integer(), md5 = character(),
      selected = logical(), stringsAsFactors = FALSE, check.names = FALSE
    )
  } else {
    do.call(rbind, candidate_rows)
  }

  list(found = found, missing = missing, candidates = candidate_df)
}

summarize_candidates <- function(candidate_df) {
  if (nrow(candidate_df) == 0) return("No candidates recorded")
  pieces <- character()
  for (nm in unique(candidate_df$whitelist_name)) {
    sub <- candidate_df[candidate_df$whitelist_name == nm, , drop = FALSE]
    n_all <- nrow(sub)
    n_md5 <- length(unique(sub$md5))
    sel <- sub$candidate_path[sub$selected][1]
    pieces <- c(pieces, paste0(nm, ": candidates=", n_all, ", md5_groups=", n_md5, ", selected=", sel))
  }
  paste(pieces, collapse = "\n")
}

safe_copy_whitelist <- function(found_list, whitelist, dest_dir) {
  copied <- setNames(character(length(whitelist)), whitelist)
  dir.create(dest_dir, recursive = TRUE, showWarnings = FALSE)

  for (nm in whitelist) {
    src <- found_list[[nm]]
    if (is.null(src) || length(src) == 0 || is.na(src) || !file.exists(src)) {
      stop(paste0("Internal error: selected source missing for ", nm), call. = FALSE)
    }
    dst <- file.path(dest_dir, nm)
    ok <- file.copy(src, dst, overwrite = TRUE, copy.date = TRUE)
    if (!isTRUE(ok)) stop(paste0("Failed to copy: ", src, " -> ", dst), call. = FALSE)
    copied[[nm]] <- dst
  }

  copied
}

xlsx_all_text <- function(path) {
  sheets <- openxlsx::getSheetNames(path)
  vals <- character()
  for (sh in sheets) {
    dat <- tryCatch(
      openxlsx::read.xlsx(path, sheet = sh, colNames = FALSE, detectDates = FALSE),
      error = function(e) data.frame()
    )
    if (length(dat) > 0) {
      vals <- c(vals, as.character(unlist(dat, use.names = FALSE)))
    }
  }
  vals <- vals[!is.na(vals)]
  vals
}

replace_not_analysed_in_xlsx <- function(path) {
  wb <- openxlsx::loadWorkbook(path)
  sheets <- openxlsx::getSheetNames(path)
  changed_total <- FALSE

  for (sh in sheets) {
    dat <- tryCatch(
      openxlsx::read.xlsx(path, sheet = sh, colNames = FALSE, detectDates = FALSE),
      error = function(e) data.frame()
    )

    if (length(dat) == 0 || nrow(dat) == 0) next

    dat2 <- dat
    changed_sheet <- FALSE

    for (j in seq_along(dat2)) {
      x <- dat2[[j]]
      y <- as.character(x)
      is_na <- is.na(x)
      y2 <- y
      y2[!is_na] <- gsub("not analysed", "not analyzed", y[!is_na], ignore.case = TRUE, perl = TRUE)
      y2[is_na] <- NA_character_
      if (!identical(y, y2)) changed_sheet <- TRUE
      dat2[[j]] <- y2
    }

    if (isTRUE(changed_sheet)) {
      changed_total <- TRUE
      openxlsx::writeData(
        wb,
        sheet = sh,
        x = dat2,
        startCol = 1,
        startRow = 1,
        colNames = FALSE,
        rowNames = FALSE,
        keepNA = FALSE
      )
    }
  }

  if (isTRUE(changed_total)) {
    openxlsx::saveWorkbook(wb, path, overwrite = TRUE)
  }

  changed_total
}

read_xlsx_first_sheet <- function(path) {
  sheets <- openxlsx::getSheetNames(path)
  openxlsx::read.xlsx(path, sheet = sheets[1], colNames = TRUE, detectDates = FALSE)
}

count_response_groups <- function(df) {
  best <- list(column = NA_character_, response = NA_integer_, non_response = NA_integer_)
  best_total <- -1L

  for (nm in names(df)) {
    vals <- tolower(trimws(as.character(df[[nm]])))
    vals[is.na(vals)] <- ""
    n_response <- sum(vals == "response")
    n_non_response <- sum(vals %in% c("non_response", "non-response", "non response"))

    if ((n_response + n_non_response) > best_total) {
      best_total <- n_response + n_non_response
      best <- list(column = nm, response = n_response, non_response = n_non_response)
    }
  }

  best
}

count_duplicate_geo_ids <- function(df) {
  nms <- names(df)
  id_candidates <- grep("(geo.*sample.*id|geo_sample|geo.accession|geo_accession|gsm|sample.*id)", nms, ignore.case = TRUE, value = TRUE)

  for (nm in id_candidates) {
    vals <- trimws(as.character(df[[nm]]))
    vals <- vals[!is.na(vals) & vals != ""]
    if (length(vals) > 0 && any(grepl("^GSM[0-9]+", vals, ignore.case = TRUE))) {
      return(sum(duplicated(vals)))
    }
  }

  flag_candidates <- grep("duplicate", nms, ignore.case = TRUE, value = TRUE)
  for (nm in flag_candidates) {
    vals <- tolower(trimws(as.character(df[[nm]])))
    vals[is.na(vals)] <- ""
    positive <- vals %in% c("true", "yes", "y", "1", "duplicate", "duplicated")
    positive <- positive | (grepl("duplicat", vals) & !grepl("no|none|not|unique|false|0", vals))
    if (length(vals) > 0) return(sum(positive))
  }

  return(NA_integer_)
}

relative_to_staging <- function(path, staging_dir) {
  p <- normalize_slash(path, mustWork = TRUE)
  s <- normalize_slash(staging_dir, mustWork = TRUE)
  prefix <- paste0(s, "/")
  if (identical(p, s)) return("")
  if (startsWith(p, prefix)) return(substr(p, nchar(prefix) + 1, nchar(p)))
  return(p)
}

file_manifest_rows <- function(paths, type, source = "generated", staging_dir) {
  data.frame(
    Type = type,
    Path = vapply(paths, relative_to_staging, character(1), staging_dir = staging_dir),
    Source = source,
    Size_bytes = as.integer(file.info(paths)$size),
    MD5 = as.character(tools::md5sum(paths)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

write_manifest <- function(manifest_path, manifest_df, candidate_df) {
  selected_candidate_df <- candidate_df
  if (nrow(selected_candidate_df) > 0) {
    selected_candidate_df$selected <- ifelse(selected_candidate_df$selected, "YES", "NO")
  }
  lines <- c(
    "# FINAL_GITHUB_STAGING_FILE_MANIFEST",
    "",
    paste0("Generated: `", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z"), "`"),
    "",
    paste0("Target GitHub repository: `", repo_url, "`"),
    "",
    "This manifest records the files copied or generated into the local GitHub repository staging folder.",
    "It does not imply that the repository has been uploaded.",
    "",
    "## Staged files",
    "",
    md_table(manifest_df),
    "",
    "## Source candidates for whitelisted copied files",
    "",
    "When multiple exact-basename candidates were found, the script selected one deterministically using the preferred final-asset folder order, then newest modification time, then path length/alphabetical order. This section is an audit trail only.",
    "",
    md_table(selected_candidate_df),
    ""
  )
  write_utf8(lines, manifest_path)
}

# -----------------------------
# 2. Create clean staging tree
# -----------------------------

if (!dir.exists(project_root)) {
  stop(
    paste0(
      "Project root does not exist: ", project_root,
      "\nCheck the D: drive path and rerun."
    ),
    call. = FALSE
  )
}

project_root <- normalize_slash(project_root, mustWork = TRUE)
staging_dir <- file.path(project_root, staging_dir_name)

if (dir.exists(staging_dir)) {
  unlink(staging_dir, recursive = TRUE, force = TRUE)
}
dir.create(staging_dir, recursive = TRUE, showWarnings = FALSE)

dir_scripts <- file.path(staging_dir, "scripts")
dir_metadata <- file.path(staging_dir, "metadata")
dir_final_figures <- file.path(staging_dir, "outputs", "final_figures")
dir_final_tables <- file.path(staging_dir, "outputs", "final_tables")
dir_session <- file.path(staging_dir, "session_info")
dir_docs <- file.path(staging_dir, "docs")

for (d in c(dir_scripts, dir_metadata, dir_final_figures, dir_final_tables, dir_session, dir_docs)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

add_qa("Staging folder created", staging_dir, dir.exists(staging_dir))
add_qa(
  "Required folder structure",
  paste(c("scripts/", "metadata/", "outputs/final_figures/", "outputs/final_tables/", "session_info/", "docs/"), collapse = "; "),
  all(dir.exists(c(dir_scripts, dir_metadata, dir_final_figures, dir_final_tables, dir_session, dir_docs)))
)

# Copy this staging script into scripts/ when R can detect its source path.
this_script <- get_this_script_path()
self_script_dst <- NA_character_
if (!is.na(this_script) && file.exists(this_script)) {
  self_script_dst <- file.path(dir_scripts, "create_github_repository_staging_a18_v2.R")
  file.copy(this_script, self_script_dst, overwrite = TRUE)
}

# -----------------------------
# 3. Find whitelisted final files
# -----------------------------

fig_search <- find_whitelist_files(project_root, figure_whitelist, staging_dir)
tab_search <- find_whitelist_files(project_root, table_whitelist, staging_dir)
source_candidate_df <- rbind(fig_search$candidates, tab_search$candidates)

if (length(fig_search$missing) > 0 || length(tab_search$missing) > 0) {
  if (length(fig_search$missing) > 0) {
    add_qa("Missing whitelisted figure file(s)", paste(fig_search$missing, collapse = "; "), FALSE)
  } else {
    add_qa("Missing whitelisted figure file(s)", "None", TRUE)
  }

  if (length(tab_search$missing) > 0) {
    add_qa("Missing whitelisted table file(s)", paste(tab_search$missing, collapse = "; "), FALSE)
  } else {
    add_qa("Missing whitelisted table file(s)", "None", TRUE)
  }

  write_qa_report(staging_dir, final_note = "One or more whitelisted files were missing. No GitHub upload was attempted.")
  stop("Staging failed because at least one whitelisted file was missing. See QA report.", call. = FALSE)
}

add_qa("Missing whitelisted figure file(s)", "None", TRUE)
add_qa("Missing whitelisted table file(s)", "None", TRUE)
add_qa("All whitelisted figures found", paste(figure_whitelist, collapse = "; "), TRUE)
add_qa("All whitelisted tables found", paste(table_whitelist, collapse = "; "), TRUE)

nonidentical_figure_dup <- 0L
nonidentical_table_dup <- 0L
if (nrow(fig_search$candidates) > 0) {
  for (nm in unique(fig_search$candidates$whitelist_name)) {
    sub <- fig_search$candidates[fig_search$candidates$whitelist_name == nm, , drop = FALSE]
    if (nrow(sub) > 1 && length(unique(sub$md5)) > 1) nonidentical_figure_dup <- nonidentical_figure_dup + 1L
  }
}
if (nrow(tab_search$candidates) > 0) {
  for (nm in unique(tab_search$candidates$whitelist_name)) {
    sub <- tab_search$candidates[tab_search$candidates$whitelist_name == nm, , drop = FALSE]
    if (nrow(sub) > 1 && length(unique(sub$md5)) > 1) nonidentical_table_dup <- nonidentical_table_dup + 1L
  }
}

add_qa(
  "Exact-basename duplicate figure candidates handled",
  paste0("non-identical duplicate whitelist basenames=", nonidentical_figure_dup, "\n", summarize_candidates(fig_search$candidates)),
  TRUE
)
add_qa(
  "Exact-basename duplicate table candidates handled",
  paste0("non-identical duplicate whitelist basenames=", nonidentical_table_dup, "\n", summarize_candidates(tab_search$candidates)),
  TRUE
)

# -----------------------------
# 4. Copy only whitelisted figures/tables
# -----------------------------

copied_figures <- safe_copy_whitelist(fig_search$found, figure_whitelist, dir_final_figures)
copied_tables <- safe_copy_whitelist(tab_search$found, table_whitelist, dir_final_tables)

add_qa("Whitelisted figures copied", paste(basename(copied_figures), collapse = "; "), TRUE)
add_qa("Whitelisted tables copied", paste(basename(copied_tables), collapse = "; "), TRUE)

# -----------------------------
# 5. XLSX patch/check dependencies
# -----------------------------

ensure_openxlsx(staging_dir)

# Patch Table 1 staging copy only.
table1_path <- file.path(dir_final_tables, "Table_1_clinical_bulk_cohorts_and_prespecified_endpoints_LOCKED.xlsx")
table1_changed <- replace_not_analysed_in_xlsx(table1_path)
add_qa(
  "Table 1 staging spelling patch",
  if (isTRUE(table1_changed)) "Patched staging copy: not analysed -> not analyzed" else "No British spelling found; no patch needed",
  TRUE
)

# -----------------------------
# 6. Generate README/docs/session info
# -----------------------------

readme_path <- file.path(staging_dir, "README.md")
repro_path <- file.path(dir_docs, "reproducibility_notes.md")
excluded_path <- file.path(dir_docs, "excluded_uncertain_files.md")
session_path <- file.path(dir_session, "R_session_info.txt")
manifest_path <- file.path(dir_metadata, "FINAL_GITHUB_STAGING_FILE_MANIFEST.md")

readme_lines <- c(
  "# HCC-TACE-hypoxia-ecotype",
  "",
  paste0("Target repository: ", repo_url),
  "",
  "This local staging folder supports the manuscript:",
  "",
  "> A single-cell and spatially anchored hypoxia-adaptive tumor-endothelial-macrophage ecotype associated with mRECIST-defined TACE non-response and adverse post-TACE outcome in hepatocellular carcinoma.",
  "",
  "## Repository scope",
  "",
  "This staging folder contains only final locked manuscript-supporting figures, final locked tables, reproducibility notes, session information, and a file manifest prepared for repository upload.",
  "",
  "Raw GEO matrices are not redistributed. Users should download public GEO datasets directly from GEO and then run verified local workflows against those public inputs.",
  "",
  "## Interpretation boundaries",
  "",
  "- The fixed ecotype was not endpoint-optimized.",
  "- GSE14520/Fako is not a response-validation cohort.",
  "- S5/S6 are audit/feasibility tables, not validation or adjusted inference.",
  "- This is not a prediction model or treatment-selection tool.",
  "- The repository is intended for transparent reproduction of reported association analyses from public datasets, not for clinical deployment.",
  "",
  "## Included final figures",
  "",
  paste0("- `outputs/final_figures/", figure_whitelist, "`"),
  "",
  "## Included final tables",
  "",
  paste0("- `outputs/final_tables/", table_whitelist, "`"),
  "",
  "## Not redistributed",
  "",
  "- Large raw expression matrices.",
  "- Private data.",
  "- Identifiable participant-level information.",
  "- Draft, candidate, failed, temporary, old, or backup figure/table objects.",
  "- Supplementary Figure S5.",
  "",
  "## Notes for users",
  "",
  "The `scripts/` directory is created for repository organization. This staging script does not auto-copy non-whitelisted historical analysis scripts because older project folders may contain draft, failed, or version-conflicting objects. Add only verified final scripts after a separate code audit."
)

write_utf8(readme_lines, readme_path)

repro_lines <- c(
  "# Reproducibility notes",
  "",
  paste0("Target repository: ", repo_url),
  "",
  "## Data redistribution boundary",
  "",
  "Raw GEO matrices are not redistributed. Users should download public GEO datasets directly from GEO accessions cited in the manuscript and then apply the documented local workflows.",
  "",
  "## Fixed ecotype boundary",
  "",
  "The fixed ecotype was not endpoint-optimized. Response labels, outcome data, or disease-context externalization results were not used to add, remove, reweight, replace, or optimize ecotype genes.",
  "",
  "## Cohort-role boundary",
  "",
  "GSE104580 is the pretreatment mRECIST-defined response-association cohort. GSE14520/Fako is not a response-validation cohort; it is used for treatment-context post-TACE outcome association and a resection-only specificity-boundary comparison.",
  "",
  "## Supplementary Table S5/S6 boundary",
  "",
  "S5/S6 are audit/feasibility tables, not validation or adjusted inference. Supplementary Table S5 documents GSE104580 sample-level endpoint harmonization and inclusion. Supplementary Table S6 documents clinical covariate availability and adjustment-feasibility boundaries.",
  "",
  "## Clinical interpretation boundary",
  "",
  "This is not a prediction model or treatment-selection tool. The score is retained as a continuous biological association measure and should not be interpreted as a cutoff, classifier, nomogram, calibration framework, treatment-allocation rule, or clinical decision-support system.",
  "",
  "## Files deliberately excluded from this staging folder",
  "",
  "Non-whitelisted figures/tables, draft candidates, failed attempts, temporary objects, old backups, Supplementary Figure S5 objects, raw matrices, and performance-model artifacts are excluded."
)

write_utf8(repro_lines, repro_path)

excluded_lines <- c(
  "# Excluded uncertain files",
  "",
  "This staging folder is intentionally conservative.",
  "",
  "The script copies only the explicitly whitelisted final locked figures and final locked tables. It deliberately excludes all non-whitelisted files, including:",
  "",
  "- Supplementary Figure S5 and any attempted S5 objects.",
  "- Figure 2 v12/v13/v14 or any non-v11 Figure 2 object.",
  "- Draft, candidate, failed, temporary, old, or backup files.",
  "- ROC/AUC/cutoff/nomogram/calibration/decision-curve/KM high-low/TCGA-validation artifacts.",
  "- Raw GEO expression matrices.",
  "- Private or identifiable participant-level data.",
  "",
  "Reason: the active manuscript is an endpoint-first translational association study, not a prediction model, classifier, treatment-selection tool, LASSO/KM risk model, or pan-HCC prognostic model."
)

write_utf8(excluded_lines, excluded_path)

session_lines <- c(
  "R session information for a18 GitHub repository staging",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  paste0("Project root: ", project_root),
  paste0("Staging folder: ", staging_dir),
  paste0("Target repository: ", repo_url),
  "",
  "sessionInfo():",
  capture.output(sessionInfo())
)
write_utf8(session_lines, session_path)

# -----------------------------
# 7. Manifest generation
# -----------------------------

manifest_entries <- list()

manifest_entries[[length(manifest_entries) + 1]] <- file_manifest_rows(
  copied_figures,
  type = "final_figure",
  source = unname(unlist(fig_search$found[figure_whitelist])),
  staging_dir = staging_dir
)

manifest_entries[[length(manifest_entries) + 1]] <- file_manifest_rows(
  copied_tables,
  type = "final_table",
  source = unname(unlist(tab_search$found[table_whitelist])),
  staging_dir = staging_dir
)

generated_paths <- c(readme_path, repro_path, excluded_path, session_path)
if (!is.na(self_script_dst) && file.exists(self_script_dst)) {
  generated_paths <- c(generated_paths, self_script_dst)
}

manifest_entries[[length(manifest_entries) + 1]] <- file_manifest_rows(
  generated_paths,
  type = "generated_metadata_or_note",
  source = "generated by create_github_repository_staging_a18_v2.R",
  staging_dir = staging_dir
)

manifest_df <- do.call(rbind, manifest_entries)
write_manifest(manifest_path, manifest_df, source_candidate_df)

# -----------------------------
# 8. Required QA checks
# -----------------------------

final_figure_files <- list.files(dir_final_figures, full.names = FALSE, recursive = FALSE)
final_table_files <- list.files(dir_final_tables, full.names = FALSE, recursive = FALSE)

add_qa(
  "Figure 1 = v14",
  paste(final_figure_files[grepl("^Figure_1_", final_figure_files)], collapse = "; "),
  identical(sort(final_figure_files[grepl("^Figure_1_", final_figure_files)]), "Figure_1_endpoint_first_workflow_v14.pdf")
)

add_qa(
  "Figure 2 = v11",
  paste(final_figure_files[grepl("^Figure_2_", final_figure_files)], collapse = "; "),
  identical(sort(final_figure_files[grepl("^Figure_2_", final_figure_files)]), "Figure_2_ecotype_construction_v11.pdf")
)

sup_fig_s5_count <- sum(grepl("Supplementary[_ -]?Figure[_ -]?S5|Supplementary_Figure_S5|Figure_S5", final_figure_files, ignore.case = TRUE))
add_qa("Supplementary Figure S5 file count", sup_fig_s5_count, sup_fig_s5_count == 0)

add_qa(
  "Only whitelisted final figures present",
  paste(final_figure_files, collapse = "; "),
  setequal(final_figure_files, figure_whitelist)
)

add_qa(
  "Only whitelisted final tables present",
  paste(final_table_files, collapse = "; "),
  setequal(final_table_files, table_whitelist)
)

main_table_present <- all(c(
  "Table_1_clinical_bulk_cohorts_and_prespecified_endpoints_LOCKED.xlsx",
  "Table_2_fixed_ecotype_components_and_transfer_rules_LOCKED.xlsx",
  "Table_3_locked_association_estimates_across_cohorts_LOCKED.xlsx"
) %in% final_table_files)

add_qa("Table 1-3 exist", paste(final_table_files[grepl("^Table_[123]_", final_table_files)], collapse = "; "), main_table_present)

supp_table_required <- paste0("Supplementary_Table_S", 1:6)
supp_table_present <- vapply(
  supp_table_required,
  function(prefix) any(startsWith(final_table_files, prefix)),
  logical(1)
)

add_qa(
  "Supplementary Table S1-S6 exist",
  paste(names(supp_table_present), supp_table_present, sep = "=", collapse = "; "),
  all(supp_table_present)
)

table1_text <- xlsx_all_text(table1_path)
table1_not_analysed_count <- sum(grepl("not analysed", table1_text, ignore.case = TRUE))
add_qa("Table 1 has no 'not analysed'", paste0("count=", table1_not_analysed_count), table1_not_analysed_count == 0)

s5_path <- file.path(dir_final_tables, "Supplementary_Table_S5_GSE104580_sample_level_endpoint_harmonization_inclusion_audit_LOCKED.xlsx")
s5_df <- read_xlsx_first_sheet(s5_path)
s5_counts <- count_response_groups(s5_df)
s5_dup <- count_duplicate_geo_ids(s5_df)

s5_pass <- nrow(s5_df) == 147 &&
  identical(as.integer(s5_counts$response), 81L) &&
  identical(as.integer(s5_counts$non_response), 66L) &&
  !is.na(s5_dup) &&
  identical(as.integer(s5_dup), 0L)

add_qa(
  "S5 audit counts",
  paste0(
    "rows=", nrow(s5_df),
    "; response=", s5_counts$response,
    "; non_response=", s5_counts$non_response,
    "; response_column=", s5_counts$column,
    "; duplicate_geo_sample_ids=", s5_dup
  ),
  s5_pass
)

s6_path <- file.path(dir_final_tables, "Supplementary_Table_S6_clinical_covariate_availability_adjustment_feasibility_audit_LOCKED.xlsx")
s6_df <- read_xlsx_first_sheet(s6_path)
add_qa("S6 row count", paste0("rows=", nrow(s6_df)), nrow(s6_df) == 58)

all_stage_files_rel <- list.files(staging_dir, recursive = TRUE, full.names = FALSE, all.files = TRUE, no.. = TRUE)
all_stage_basenames <- basename(all_stage_files_rel)

# These checks apply to staged filenames only. The notes may discuss excluded terms in text, which is intentional.
bad_draft_files <- all_stage_files_rel[
  grepl("(DRAFT|candidate|failed|temp|old|backup)", all_stage_basenames, ignore.case = TRUE)
]
add_qa(
  "No DRAFT/candidate/failed/temp/old/backup filenames",
  if (length(bad_draft_files) == 0) "None" else paste(bad_draft_files, collapse = "; "),
  length(bad_draft_files) == 0
)

bad_model_files <- all_stage_files_rel[
  grepl("(ROC|AUC|cutoff|cut\\-off|cut[ _-]?point|nomogram|calibration|decision[ _-]?curve|KM[ _-]?high[ _-]?low|TCGA[ _-]?validation)", all_stage_basenames, ignore.case = TRUE)
]
add_qa(
  "No ROC/AUC/cutoff/nomogram/calibration/decision curve/KM high-low/TCGA validation filenames",
  if (length(bad_model_files) == 0) "None" else paste(bad_model_files, collapse = "; "),
  length(bad_model_files) == 0
)

readme_text <- paste(readLines(readme_path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
readme_raw_pass <- grepl("Raw GEO matrices are not redistributed", readme_text, fixed = TRUE) &&
  !grepl("raw GEO matrices are redistributed|raw matrices are included|includes raw expression matrices|contains raw expression matrices", readme_text, ignore.case = TRUE)

add_qa(
  "README does not claim raw matrices are included",
  if (readme_raw_pass) "Contains explicit non-redistribution statement and no affirmative raw-matrix inclusion claim" else "README raw-matrix boundary failed",
  readme_raw_pass
)

readme_model_pass <- grepl("not a prediction model or treatment-selection tool", readme_text, ignore.case = TRUE) &&
  !grepl("\\bis a prediction model\\b|\\bis a treatment-selection tool\\b|\\bclinical prediction model\\b", readme_text, ignore.case = TRUE)

add_qa(
  "README does not claim prediction model or treatment-selection tool",
  if (readme_model_pass) "Contains explicit negative clinical-tool boundary and no affirmative model/tool claim" else "README model/tool boundary failed",
  readme_model_pass
)

required_docs <- c(readme_path, repro_path, excluded_path, session_path, manifest_path)
add_qa(
  "Required generated repository documents exist",
  paste(vapply(required_docs, relative_to_staging, character(1), staging_dir = staging_dir), collapse = "; "),
  all(file.exists(required_docs))
)

# Refresh manifest after QA report will be written separately. The manifest remains the staged-file audit.
write_qa_report(staging_dir, final_note = "No GitHub upload was performed. This script only created the local staging folder.")

if (any(qa$Result != "PASS")) {
  stop("Staging completed with QA failure(s). See a18_github_repository_staging_QA_report.md.", call. = FALSE)
}

message("")
message("PASS: GitHub repository staging folder created successfully.")
message("Staging folder: ", staging_dir)
message("Target repository URL recorded in README: ", repo_url)
message("No GitHub upload was performed.")
