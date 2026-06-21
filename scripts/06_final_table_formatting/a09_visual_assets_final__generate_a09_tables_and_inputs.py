import csv
import html
import math
import statistics
import zipfile
from pathlib import Path

import pandas as pd


ROOT = Path.cwd()
OUT = ROOT / "results" / "a09_visual_assets_final"
DATA = OUT / "data"
TABLES = OUT / "tables"
SCRIPTS = OUT / "scripts"
AUDIT = OUT / "audit"


def ensure_dirs():
    for path in [OUT, DATA, TABLES, SCRIPTS, AUDIT]:
        path.mkdir(parents=True, exist_ok=True)


def read_csv(rel_path):
    return pd.read_csv(ROOT / rel_path)


def write_csv(df, path):
    df.to_csv(path, index=False, encoding="utf-8-sig")


def col_letter(n):
    out = ""
    while n:
        n, rem = divmod(n - 1, 26)
        out = chr(65 + rem) + out
    return out


def write_xlsx(df, path, sheet_name="Sheet1"):
    rows = [list(df.columns)] + df.astype(object).where(pd.notna(df), "").values.tolist()
    sheet_rows = []
    for r_idx, row in enumerate(rows, 1):
        cells = []
        for c_idx, val in enumerate(row, 1):
            ref = f"{col_letter(c_idx)}{r_idx}"
            if isinstance(val, (int, float)) and not isinstance(val, bool) and not pd.isna(val):
                cells.append(f'<c r="{ref}"><v>{val}</v></c>')
            else:
                txt = html.escape(str(val), quote=True)
                cells.append(f'<c r="{ref}" t="inlineStr"><is><t>{txt}</t></is></c>')
        sheet_rows.append(f'<row r="{r_idx}">{"".join(cells)}</row>')
    max_col = len(rows[0]) if rows else 1
    dim = f"A1:{col_letter(max_col)}{len(rows)}"
    sheet_xml = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">
<dimension ref="{dim}"/><sheetData>{"".join(sheet_rows)}</sheetData></worksheet>'''
    content_types = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">
<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>
<Default Extension="xml" ContentType="application/xml"/>
<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>
<Override PartName="/xl/worksheets/sheet1.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>
<Override PartName="/docProps/core.xml" ContentType="application/vnd.openxmlformats-package.core-properties+xml"/>
<Override PartName="/docProps/app.xml" ContentType="application/vnd.openxmlformats-officedocument.extended-properties+xml"/>
</Types>'''
    rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>
</Relationships>'''
    wb = f'''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main" xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">
<sheets><sheet name="{html.escape(sheet_name[:31])}" sheetId="1" r:id="rId1"/></sheets></workbook>'''
    wb_rels = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">
<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet1.xml"/>
</Relationships>'''
    app = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><Properties xmlns="http://schemas.openxmlformats.org/officeDocument/2006/extended-properties"><Application>Codex</Application></Properties>'''
    core = '''<?xml version="1.0" encoding="UTF-8" standalone="yes"?><cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" xmlns:dc="http://purl.org/dc/elements/1.1/"><dc:creator>Codex</dc:creator></cp:coreProperties>'''
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as z:
        z.writestr("[Content_Types].xml", content_types)
        z.writestr("_rels/.rels", rels)
        z.writestr("xl/workbook.xml", wb)
        z.writestr("xl/_rels/workbook.xml.rels", wb_rels)
        z.writestr("xl/worksheets/sheet1.xml", sheet_xml)
        z.writestr("docProps/app.xml", app)
        z.writestr("docProps/core.xml", core)


def write_table(name, df, legend):
    write_csv(df, TABLES / f"{name}.csv")
    write_xlsx(df, TABLES / f"{name}.xlsx", sheet_name=name[:31])
    (TABLES / f"{name}_legend.md").write_text(legend + "\n", encoding="utf-8")


def fmt_p(x):
    if pd.isna(x) or x == "":
        return ""
    x = float(x)
    if x < 0.001:
        return f"{x:.3g}"
    return f"{x:.6g}"


def summarize_scores(df, score_col):
    vals = pd.to_numeric(df[score_col], errors="coerce").dropna().tolist()
    if not vals:
        return {"n": 0, "mean": "", "sd": "", "median": "", "iqr": ""}
    q1 = pd.Series(vals).quantile(0.25)
    q3 = pd.Series(vals).quantile(0.75)
    return {
        "n": len(vals),
        "mean": round(statistics.mean(vals), 6),
        "sd": round(statistics.stdev(vals), 6) if len(vals) > 1 else "",
        "median": round(statistics.median(vals), 6),
        "iqr": f"{q1:.6g}-{q3:.6g}",
    }


def main():
    ensure_dirs()

    a02 = read_csv("results/a02_ecotype_construction/A02_final_ecotype_gene_set.csv")
    fig2_matrix = read_csv("results/a07_R_plotting_inputs/figure2_component_gene_score_matrix.csv")
    fig3_scores = read_csv("results/a07_R_plotting_inputs/figure3_gse104580_sample_scores.csv")
    fig4 = read_csv("results/a07_R_plotting_inputs/figure4_gse14520_outcome_forest.csv")
    gse76427_sample = read_csv("results/a07_R_plotting_inputs/suppfig_s1_gse76427_sample_level_scores.csv")
    gse76427_gene = read_csv("results/a07_R_plotting_inputs/suppfig_s1_gse76427_gene_recoverability.csv")
    gse76427_cox = read_csv("results/a07_R_plotting_inputs/suppfig_s1_gse76427_cox_results.csv")
    gse76427_stage = read_csv("results/a07_R_plotting_inputs/suppfig_s1_gse76427_stage_association_results.csv")
    a03_stats = read_csv("results/a03_gse104580_response_association/A03_response_association_statistics.csv")
    a03_model = read_csv("results/a03_gse104580_response_association/A03_association_model_results.csv")
    a03_sens = read_csv("results/a03_gse104580_response_association/A03_component_consistency_and_sensitivity_results.csv")
    a04_scores = read_csv("results/a04_gse14520_treatment_context_validation/A04_GSE14520_ecotype_scores.csv")
    a04_missing = read_csv("results/a04_gse14520_treatment_context_validation/A04_endpoint_missingness_audit.csv")
    a05_scores = read_csv("results/a05_secondary_externalization/A05_ecotype_scores.csv")
    a05_gene = read_csv("results/a05_secondary_externalization/A05_gene_availability.csv")

    # Supplementary Figure S4 input.
    s4 = fig4.copy()
    s4["plot_group"] = ["A", "A", "B"]
    s4["boundary_note"] = [
        "Primary post-TACE outcome association; not a response cohort.",
        "Secondary/exploratory recurrence/RFS field.",
        "Boundary-setting comparator; not uniquely specific to the TACE-treated setting.",
    ]
    write_csv(s4, DATA / "Supplementary_Figure_S4_expanded_gse14520_outcome_data.csv")
    (DATA / "Supplementary_Figure_S4_panel_manifest.md").write_text(
        "# Supplementary Figure S4 panel/data manifest\n\n"
        "- Panel A: endpoint availability for GSE14520/Fako TACE-treated and resection-only subsets.\n"
        "- Panel B: continuous effect-size display for TACE-treated OS, exploratory TACE-treated recurrence/RFS, and resection-only OS specificity-boundary check.\n"
        "- Source: results/a07_R_plotting_inputs/figure4_gse14520_outcome_forest.csv.\n"
        "- Boundary: no risk groups, no threshold optimization, and no performance display.\n",
        encoding="utf-8",
    )

    # Supplementary Figure S5 input.
    s5_rows = [
        ["Reference anchor layer", "GSE149614/GSE151530/GSE277104", "28 component-gene rows; 27 unique genes", "not applicable", "single-cell/spatial biological anchoring", "Reference-derived anchoring only; not treatment-response evidence"],
        ["Primary response association", "GSE104580", "147/147 scored", "81 responders; 66 non-responders", "mRECIST-defined TACE non-response association", "Association only; no classification framing"],
        ["Treatment-context association", "GSE14520/Fako TACE-treated", "105 scored", "104 OS-complete; 104 recurrence/RFS available for exploratory field", "post-TACE outcome association", "Not a response validation cohort"],
        ["Specificity-boundary check", "GSE14520/Fako resection-only", "86 scored", "85 OS-complete", "specificity-boundary comparator", "Positive association indicates limited treatment-context specificity"],
        ["Rescue-limited externalization", "GSE76427", "27/27 unique genes covered; 167 scored", "115 OS-complete; 108 RFS-complete; 114 TNM; 115 BCLC", "secondary disease-context externalization", "Recoverability and boundary checks only"],
        ["Blocked resource", "TCGA-LIHC", "not scored", "not analyzed", "blocked", "No TCGA-derived claim"],
    ]
    s5 = pd.DataFrame(s5_rows, columns=[
        "evidence_layer", "dataset_or_resource", "score_or_gene_coverage", "endpoint_availability",
        "predefined_role", "boundary_note"
    ])
    write_csv(s5, DATA / "Supplementary_Figure_S5_cohort_coverage_endpoint_role_data.csv")
    (DATA / "Supplementary_Figure_S5_panel_manifest.md").write_text(
        "# Supplementary Figure S5 panel/data manifest\n\n"
        "- Matrix rows: reference layer, GSE104580, GSE14520/Fako TACE-treated, GSE14520/Fako resection-only, GSE76427, TCGA-LIHC.\n"
        "- Columns: coverage/scoring, endpoint availability, predefined analytical role, and interpretation boundary.\n"
        "- Sources: A02-A05 locked outputs and a07 plotting input files.\n"
        "- Boundary: evidence availability map only; no performance or broad validation framing.\n",
        encoding="utf-8",
    )

    # Main tables.
    table1 = pd.DataFrame([
        ["GSE149614/GSE151530/GSE277104", "single-cell pseudobulk / spatial AOI", "HCC biological reference layer", "pre-treatment biological context", "not a clinical scoring cohort", "cell-state/tissue-context anchors", "reference-derived ecotype anchoring", "supports ecotype construction", "not treatment-response evidence"],
        ["GSE104580", "GPL570 bulk expression", "pretreatment TACE biopsy cohort", "pretreatment", "147 scored", "mRECIST-defined response/non-response", "primary response association cohort", "fixed ecotype associated with TACE non-response", "association only; no classification framing"],
        ["GSE14520/Fako TACE-treated", "GPL3921/GPL571 bulk expression plus Fako clinical table", "adjuvant TACE-treated HCC subset", "resection specimen with post-treatment context", "105 scored", "OS primary; recurrence/RFS secondary/exploratory", "treatment-context outcome association", "fixed ecotype linked to adverse post-TACE outcome", "not radiological response validation"],
        ["GSE14520/Fako resection-only", "GPL3921/GPL571 bulk expression plus Fako clinical table", "resection-only HCC subset", "resection-only context", "86 scored", "OS specificity-boundary check", "specificity-boundary comparator", "tests whether association is unique to the TACE-treated setting", "positive result means limited treatment-context specificity"],
        ["GSE76427", "GPL10558 bulk expression", "HCC disease-context cohort", "surgical HCC context", "167 scored; 27/27 genes covered", "TNM/BCLC/OS/RFS rescue-limited checks", "secondary externalization only", "fixed ecotype can be recovered in broader HCC data", "not broad HCC outcome validation"],
        ["TCGA-LIHC", "blocked", "not analyzed", "not analyzed", "0", "not available", "blocked/not analyzed", "none", "no TCGA-derived claim"],
    ], columns=[
        "Dataset", "platform/data type", "disease/treatment context", "pretreatment status", "n scored",
        "available endpoint(s)", "analytical role", "allowed interpretation", "bounded/forbidden interpretation"
    ])
    write_table("Table_1_study_cohorts_endpoint_hierarchy_and_roles", table1,
                "Table 1. Study cohorts, endpoint hierarchy, and analytical role definitions.")

    component_labels = {
        "malignant_hypoxia_adaptive": "malignant hypoxia-adaptive tumor program",
        "endothelial_angiogenesis_interaction": "endothelial/angiogenesis interaction program",
        "spp1_mmp9_tam_matrix_angiogenic_support": "SPP1/MMP9-like TAM matrix/angiogenic support program",
    }
    table2_rows = []
    for comp, sub in a02.groupby("component", sort=False):
        genes = ", ".join(sub["gene_symbol"].tolist())
        unique_note = "VEGFA also appears in another component; total score uses predefined duplicate handling." if "VEGFA" in set(sub["gene_symbol"]) else "No cross-component duplicate in this component except as defined by final set."
        table2_rows.append([
            comp,
            component_labels.get(comp, comp),
            genes,
            len(sub),
            unique_note,
            "DLL4 excluded before downstream transfer; no endpoint-guided gene replacement."
        ])
    table2 = pd.DataFrame(table2_rows, columns=[
        "Component", "biological compartment/state", "final gene membership (collapsed readable display)",
        "number of component-gene rows", "unique-gene accounting", "transfer/preprocessing note"
    ])
    write_table("Table_2_final_ecotype_definition_summary", table2,
                "Table 2. Final ecotype definition summary. The fixed A02 ecotype contains 28 component-gene rows and 27 unique genes.")

    total_stats = a03_stats[a03_stats["score_name"].eq("total_ecotype_equal_component_score")].iloc[0]
    total_model = a03_model[a03_model["score_name"].eq("total_ecotype_equal_component_score")].iloc[0]
    table3_rows = [
        ["GSE104580", "TACE non-response", "Non-response versus response; fixed total score", f"median difference = {total_stats['non_response_minus_response_median']:.7f}; Cliff's delta = {total_stats['cliffs_delta_non_response_vs_response']:.6f}; OR per 1 SD = {total_model['OR_per_1sd_score']:.5f}", f"{total_model['OR_95CI_low']:.5f}-{total_model['OR_95CI_high']:.5f}", f"Wilcoxon P = {total_stats['p_value']:.8f}; logistic P = {total_model['p_value']:.8f}", "associated with TACE non-response"],
        ["GSE14520/Fako TACE-treated", "OS", "Cox association per 1 SD fixed total score", "HR = 1.55533", "1.13482-2.13167", "0.00602555", "linked to adverse post-TACE outcome"],
        ["GSE14520/Fako TACE-treated", "recurrence/RFS", "Secondary/exploratory Cox association per 1 SD", "HR = 1.27168", "1.00029-1.61670", "0.0497186", "secondary/exploratory only"],
        ["GSE14520/Fako resection-only", "OS", "Specificity-boundary Cox association per 1 SD", "HR = 1.57707", "1.13725-2.18698", "0.00631371", "limited treatment-context specificity"],
        ["GSE76427", "TNM stage", "Spearman correlation with fixed score", "r = 0.0223715", "", "0.813637", "rescue-limited boundary output"],
        ["GSE76427", "BCLC stage", "Spearman correlation with fixed score", "r = -0.18947", "", "0.0423888", "rescue-limited boundary output"],
        ["GSE76427", "OS", "Supportive Cox association per 1 SD", "HR = 0.750805", "0.513673-1.09741", "0.138868", "rescue-limited secondary externalization"],
        ["GSE76427", "RFS", "Supportive Cox association per 1 SD", "HR = 0.875147", "0.655357-1.16865", "0.3661", "rescue-limited secondary externalization"],
    ]
    table3 = pd.DataFrame(table3_rows, columns=["Cohort", "endpoint", "analysis scale/contrast", "estimate", "95% CI where applicable", "P value", "interpretation label"])
    write_table("Table_3_locked_association_estimates_across_cohorts", table3,
                "Table 3. Locked association estimates across cohorts. Estimates are reported for association and boundary interpretation only.")

    # Supplementary tables.
    sup1 = fig2_matrix.rename(columns={
        "component": "Component",
        "gene_symbol": "Gene symbol",
        "scRNA_support_score": "Single-cell support score",
        "spatial_support_score": "Spatial-context support score",
        "bulk_transfer_available": "Bulk transfer available",
        "duplicate_gene": "Duplicate gene",
        "excluded_before_transfer": "Excluded before transfer",
    }).copy()
    sup1["Row category"] = sup1["Excluded before transfer"].map(lambda x: "excluded before transfer" if str(x).lower() == "true" else "final component-gene row")
    sup1["Transfer rule"] = sup1.apply(
        lambda r: "DLL4 excluded before downstream transfer" if r["Row category"] == "excluded before transfer"
        else ("VEGFA duplicate handled by predefined rule" if str(r["Duplicate gene"]).lower() == "true" else "included in fixed bulk transfer"),
        axis=1,
    )
    write_table("Supplementary_Table_S1_full_component_provenance_and_transfer_rules", sup1,
                "Supplementary Table S1. Full component-gene provenance and preprocessing rules, including the DLL4 excluded-before-transfer row.")

    sup2 = pd.DataFrame([
        ["GSE104580", "147", "81 responders; 66 non-responders", "pretreatment TACE biopsy", "primary response association", "analyzed"],
        ["GSE14520/Fako TACE-treated", "105", "104 OS-complete; one TACE-treated sample scored but OS missing", "adjuvant TACE-treated", "treatment-context outcome association", "analyzed"],
        ["GSE14520/Fako resection-only", "86", "85 OS-complete; one resection-only sample scored but OS missing", "resection-only", "specificity-boundary check", "analyzed"],
        ["GSE76427", "167", "115 OS-complete; 108 RFS-complete; 114 TNM; 115 BCLC", "broader HCC disease context", "rescue-limited secondary externalization", "analyzed"],
        ["TCGA-LIHC", "0", "not available", "not analyzed", "blocked/not analyzed", "blocked"],
    ], columns=["Cohort", "n scored", "endpoint completeness", "pretreatment/treatment label", "analytical role", "blocked/not analyzed status"])
    write_table("Supplementary_Table_S2_cohort_manifest_endpoint_completeness", sup2,
                "Supplementary Table S2. Cohort manifest and endpoint completeness.")

    a03_sum = summarize_scores(fig3_scores.rename(columns={"fixed_total_ecotype_score": "total_ecotype_equal_component_score"}), "total_ecotype_equal_component_score")
    a04_tace_sum = summarize_scores(a04_scores[a04_scores["analysis_context"].eq("TACE_treated")], "total_ecotype_equal_component_score")
    a04_res_sum = summarize_scores(a04_scores[a04_scores["analysis_context"].eq("resection_only")], "total_ecotype_equal_component_score")
    a05_sum = summarize_scores(a05_scores, "total_ecotype_equal_component_score")
    sup3 = pd.DataFrame([
        ["GSE104580", "27/27 fixed unique genes available", a03_sum["n"], a03_sum["mean"], a03_sum["sd"], a03_sum["median"], a03_sum["iqr"], "GPL570 transfer; fixed score only"],
        ["GSE14520/Fako TACE-treated", "fixed genes transferable across GPL3921/GPL571 per A04 audit", a04_tace_sum["n"], a04_tace_sum["mean"], a04_tace_sum["sd"], a04_tace_sum["median"], a04_tace_sum["iqr"], "platform-specific transfer; no outcome-guided pruning"],
        ["GSE14520/Fako resection-only", "fixed genes transferable across GPL3921/GPL571 per A04 audit", a04_res_sum["n"], a04_res_sum["mean"], a04_res_sum["sd"], a04_res_sum["median"], a04_res_sum["iqr"], "specificity-boundary subset"],
        ["GSE76427", "27/27 fixed unique genes covered", a05_sum["n"], a05_sum["mean"], a05_sum["sd"], a05_sum["median"], a05_sum["iqr"], "rescue-limited local annotation transfer"],
        ["TCGA-LIHC", "blocked/not analyzed", 0, "", "", "", "", "no TCGA-derived score"],
    ], columns=["Cohort", "gene coverage", "score n", "score mean", "score SD", "score median", "score IQR", "platform/transfer note"])
    write_table("Supplementary_Table_S3_per_cohort_score_coverage_descriptive_summaries", sup3,
                "Supplementary Table S3. Per-cohort score coverage and descriptive summaries. Descriptive summaries do not alter locked association results.")

    sup4_rows = [
        ["GSE14520/Fako resection-only", "OS", "HR per 1 SD", 1.57707, 1.13725, 2.18698, 0.00631371, "specificity-boundary comparator", "positive association indicates limited treatment-context specificity"],
        ["GSE76427", "TNM stage", "Spearman r", 0.0223715, "", "", 0.813637, "rescue-limited disease-context check", "not adverse-direction supportive"],
        ["GSE76427", "BCLC stage", "Spearman r", -0.18947, "", "", 0.0423888, "rescue-limited disease-context check", "negative direction; do not interpret as adverse-direction support"],
        ["GSE76427", "OS", "HR per 1 SD", 0.750805, 0.513673, 1.09741, 0.138868, "rescue-limited supportive outcome check", "not broad HCC outcome validation"],
        ["GSE76427", "RFS", "HR per 1 SD", 0.875147, 0.655357, 1.16865, 0.3661, "rescue-limited supportive outcome check", "not broad HCC outcome validation"],
        ["TCGA-LIHC", "all", "not analyzed", "", "", "", "", "blocked", "no TCGA-derived claim"],
    ]
    sup4 = pd.DataFrame(sup4_rows, columns=["Cohort", "Endpoint/context", "Estimate type", "Estimate", "CI lower", "CI upper", "P value", "Allowed interpretation", "Boundary note"])
    write_table("Supplementary_Table_S4_externalization_and_boundary_outputs", sup4,
                "Supplementary Table S4. Expanded resection-only and GSE76427 externalization/boundary outputs.")

    audit_inputs = [
        "results/a07_R_plotting_inputs/figure4_gse14520_outcome_forest.csv",
        "results/a07_R_plotting_inputs/figure2_component_gene_score_matrix.csv",
        "results/a07_R_plotting_inputs/figure3_gse104580_sample_scores.csv",
        "results/a07_R_plotting_inputs/suppfig_s1_gse76427_sample_level_scores.csv",
        "results/a07_R_plotting_inputs/suppfig_s1_gse76427_gene_recoverability.csv",
        "results/a03_gse104580_response_association/A03_response_association_statistics.csv",
        "results/a03_gse104580_response_association/A03_association_model_results.csv",
        "results/a04_gse14520_treatment_context_validation/A04_GSE14520_ecotype_scores.csv",
        "results/a05_secondary_externalization/A05_ecotype_scores.csv",
        "results/a05_secondary_externalization/A05_gene_availability.csv",
    ]
    outputs = sorted([str(p.relative_to(ROOT)) for p in OUT.rglob("*") if p.is_file()])
    audit = [
        "# Codex visual assets execution audit",
        "",
        "## Scope",
        "Only new Supplementary Figures S4-S5, Tables 1-3, Supplementary Tables S1-S4, source scripts, data manifests, and this audit were generated.",
        "",
        "No manuscript DOCX was edited. Frozen v11 main figures and frozen S1-S3 supplementary figures were not redesigned.",
        "",
        "## Input files used",
    ]
    audit += [f"- {x}" for x in audit_inputs]
    audit += [
        "",
        "## Locked-result consistency",
        "- A02: 28 component-gene rows and 27 unique genes were retained; VEGFA duplicate handling and DLL4 exclusion were represented without modifying the gene set.",
        "- A03: GSE104580 locked response association values were copied from A03 outputs.",
        "- A04: GSE14520/Fako TACE-treated and resection-only locked estimates were copied from A04/a07 plotting inputs.",
        "- A05: GSE76427 locked rescue-limited values were copied from A05/a07 repaired plotting inputs; TCGA-LIHC remains blocked/not analyzed.",
        "",
        "## Output files",
    ]
    audit += [f"- {x}" for x in outputs]
    audit += [
        "",
        "## Confirmations",
        "- No new endpoint, cohort, gene-set, or statistical model was introduced.",
        "- No threshold-optimized or model-performance visual was generated.",
        "- S4 uses continuous effect-size-oriented display only.",
        "- S5 is an evidence/availability matrix and does not imply broad validation.",
        "",
        "## Missing fields/manual author resolution",
        "- The v2_18 manuscript DOCX and v2_18 audit markdown named in the execution instruction were not found under D:/Downloads at runtime; asset generation used locked project results instead.",
        "- True scRNA UMAP coordinates remain unavailable in the plotting inputs and were not fabricated or recomputed.",
    ]
    (AUDIT / "Codex_visual_assets_execution_audit.md").write_text("\n".join(audit) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
