import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
A01 = ROOT / "results" / "a01_data_preparation"
A02 = ROOT / "results" / "a02_ecotype_construction"
A04 = ROOT / "results" / "a04_gse14520_treatment_context_validation"

COMPONENTS = [
    "malignant_hypoxia_adaptive",
    "endothelial_angiogenesis_interaction",
    "spp1_mmp9_tam_matrix_angiogenic_support",
]
PRIMARY_SCORE = "total_ecotype_equal_component_score"


def write_csv(path, header, rows):
    with (A04 / path).open("w", encoding="utf-8", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh))


def read_tsv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load_expression(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader)
        samples = header[1:]
        expr = {}
        for row in reader:
            if row:
                expr[row[0]] = np.array([float(x) for x in row[1:]], dtype=float)
    return samples, expr


def rankdata(values):
    values = list(values)
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg
        i = j + 1
    return np.array(ranks, dtype=float)


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def fmt(x):
    if x is None or x == "":
        return ""
    if isinstance(x, str):
        return x
    if isinstance(x, float) and (math.isnan(x) or math.isinf(x)):
        return ""
    return f"{x:.6g}"


def rank_score(samples, expr, genes):
    measured = [g for g in genes if g in expr]
    if not measured:
        return np.full(len(samples), np.nan), []
    all_genes = list(expr)
    matrix = np.vstack([expr[g] for g in all_genes])
    ranks = np.zeros_like(matrix)
    for j in range(matrix.shape[1]):
        ranks[:, j] = rankdata(matrix[:, j])
    gene_to_idx = {g: i for i, g in enumerate(all_genes)}
    idx = [gene_to_idx[g] for g in measured]
    raw = ranks[idx, :].mean(axis=0)
    return (raw - 1.0) / max(len(all_genes) - 1.0, 1.0), measured


def cox_univariate(time, event, score):
    time = np.asarray(time, dtype=float)
    event = np.asarray(event, dtype=int)
    x = np.asarray(score, dtype=float)
    keep = np.isfinite(time) & np.isfinite(x) & (time > 0) & ((event == 0) | (event == 1))
    time = time[keep]
    event = event[keep]
    x = x[keep]
    if len(time) < 10 or event.sum() < 2 or x.std(ddof=1) == 0:
        return {"n": len(time), "events": int(event.sum()), "beta": "", "se": "", "hr": "", "ci_low": "", "ci_high": "", "p": "", "direction": "not_estimable"}
    z = (x - x.mean()) / x.std(ddof=1)
    beta = 0.0
    for _ in range(80):
        score_u = 0.0
        info = 0.0
        for ti, ei, zi in zip(time, event, z):
            if ei != 1:
                continue
            risk = time >= ti
            xr = z[risk]
            w = np.exp(beta * xr)
            sw = w.sum()
            mean_x = (w * xr).sum() / sw
            mean_x2 = (w * xr * xr).sum() / sw
            score_u += zi - mean_x
            info += mean_x2 - mean_x * mean_x
        if info <= 0:
            break
        step = score_u / info
        beta += step
        if abs(step) < 1e-8:
            break
    # Observed information at final beta
    info = 0.0
    for ti, ei in zip(time, event):
        if ei != 1:
            continue
        xr = z[time >= ti]
        w = np.exp(beta * xr)
        sw = w.sum()
        mean_x = (w * xr).sum() / sw
        mean_x2 = (w * xr * xr).sum() / sw
        info += mean_x2 - mean_x * mean_x
    if info <= 0:
        return {"n": len(time), "events": int(event.sum()), "beta": "", "se": "", "hr": "", "ci_low": "", "ci_high": "", "p": "", "direction": "not_estimable"}
    se = math.sqrt(1.0 / info)
    zstat = beta / se
    p = 2 * norm_sf(abs(zstat))
    return {
        "n": len(time),
        "events": int(event.sum()),
        "beta": beta,
        "se": se,
        "hr": math.exp(beta),
        "ci_low": math.exp(beta - 1.96 * se),
        "ci_high": math.exp(beta + 1.96 * se),
        "p": p,
        "direction": "higher_score_higher_hazard" if beta > 0 else "higher_score_lower_hazard",
    }


def parse_float(x):
    try:
        if x is None or x == "":
            return None
        return float(x)
    except ValueError:
        return None


def parse_event(x):
    if x in ("0", "1"):
        return int(x)
    return None


def main():
    A04.mkdir(parents=True, exist_ok=True)
    geneset = read_csv(A02 / "A02_final_ecotype_gene_set.csv")
    component_genes = defaultdict(list)
    gene_components = defaultdict(list)
    for r in geneset:
        component_genes[r["component"]].append(r["gene_symbol"])
        gene_components[r["gene_symbol"]].append(r["component"])
    unique_genes = sorted(gene_components)
    duplicates = {g: c for g, c in gene_components.items() if len(c) > 1}

    clinical = read_tsv(A01 / "GSE14520_Fako_mapped_clinical_metadata.tsv")
    clinical_by_gsm = {r["gsm"]: r for r in clinical if r["gsm"]}

    platform_inputs = {
        "GPL3921": A01 / "GSE14520_GPL3921_gene_expression.tsv",
        "GPL571": A01 / "GSE14520_GPL571_gene_expression.tsv",
    }
    score_rows = []
    availability_rows = []
    platform_score_arrays = {}
    for platform, path in platform_inputs.items():
        samples, expr = load_expression(path)
        scores = {}
        measured = {}
        for component in COMPONENTS:
            scores[component], measured[component] = rank_score(samples, expr, component_genes[component])
            for g in component_genes[component]:
                availability_rows.append([platform, component, g, "yes" if g in expr else "no"])
        scores["total_ecotype_unique_gene_score"], measured["total_unique"] = rank_score(samples, expr, unique_genes)
        scores[PRIMARY_SCORE] = np.vstack([scores[c] for c in COMPONENTS]).mean(axis=0)
        platform_score_arrays[platform] = (samples, scores)
        for i, gsm in enumerate(samples):
            clin = clinical_by_gsm.get(gsm)
            if not clin:
                continue
            score_rows.append([
                gsm,
                platform,
                clin["lcs_id"],
                clin["fako_therapy_group"],
                clin["analysis_context"],
                clin["os_status"],
                clin["os_months"],
                clin["recurrence_status"],
                clin["recurrence_months"],
                scores["malignant_hypoxia_adaptive"][i],
                scores["endothelial_angiogenesis_interaction"][i],
                scores["spp1_mmp9_tam_matrix_angiogenic_support"][i],
                scores[PRIMARY_SCORE][i],
                scores["total_ecotype_unique_gene_score"][i],
            ])
    write_csv(
        "A04_GSE14520_gene_availability_by_platform.csv",
        ["platform", "component", "gene_symbol", "available"],
        availability_rows,
    )
    score_header = [
        "gsm", "platform", "lcs_id", "fako_therapy_group", "analysis_context",
        "os_status", "os_months", "recurrence_status", "recurrence_months",
        "malignant_hypoxia_adaptive_score", "endothelial_angiogenesis_interaction_score",
        "spp1_mmp9_tam_matrix_angiogenic_support_score", PRIMARY_SCORE,
        "total_ecotype_unique_gene_score",
    ]
    write_csv("A04_GSE14520_ecotype_scores.csv", score_header, score_rows)

    score_dicts = [dict(zip(score_header, r)) for r in score_rows]
    flow_counts = Counter(r["analysis_context"] for r in score_dicts)
    therapy_counts = Counter(r["fako_therapy_group"] for r in score_dicts)
    platform_context_counts = Counter((r["platform"], r["analysis_context"]) for r in score_dicts)

    def subset(context):
        return [r for r in score_dicts if r["analysis_context"] == context]

    def run_cox(rows, score_col, time_col, event_col):
        time = [parse_float(r[time_col]) for r in rows]
        event = [parse_event(r[event_col]) for r in rows]
        score = [parse_float(r[score_col]) for r in rows]
        keep = [t is not None and e is not None and s is not None for t, e, s in zip(time, event, score)]
        return cox_univariate([t for t, k in zip(time, keep) if k], [e for e, k in zip(event, keep) if k], [s for s, k in zip(score, keep) if k])

    score_cols = [
        "malignant_hypoxia_adaptive_score",
        "endothelial_angiogenesis_interaction_score",
        "spp1_mmp9_tam_matrix_angiogenic_support_score",
        PRIMARY_SCORE,
        "total_ecotype_unique_gene_score",
    ]
    tace_rows = subset("TACE_treated")
    resection_rows = subset("resection_only")

    os_rows = []
    for score_col in score_cols:
        res = run_cox(tace_rows, score_col, "os_months", "os_status")
        os_rows.append(["TACE_treated", "OS", score_col, res["n"], res["events"], res["hr"], res["ci_low"], res["ci_high"], res["p"], res["direction"], "Cox_association_not_prediction"])
    for platform in ["GPL3921", "GPL571"]:
        rows = [r for r in tace_rows if r["platform"] == platform]
        res = run_cox(rows, PRIMARY_SCORE, "os_months", "os_status")
        os_rows.append([f"TACE_treated_{platform}", "OS", PRIMARY_SCORE, res["n"], res["events"], res["hr"], res["ci_low"], res["ci_high"], res["p"], res["direction"], "platform_consistency"])
    for leave_out in COMPONENTS:
        kept = [c for c in COMPONENTS if c != leave_out]
        col = f"leave_out_{leave_out}_score"
        for r in score_dicts:
            vals = [parse_float(r[f"{c}_score"]) for c in kept]
            r[col] = sum(vals) / len(vals)
        res = run_cox(tace_rows, col, "os_months", "os_status")
        os_rows.append(["TACE_treated", "OS", col, res["n"], res["events"], res["hr"], res["ci_low"], res["ci_high"], res["p"], res["direction"], "leave_one_component_out"])
    write_csv(
        "A04_TACE_treated_OS_association_results.csv",
        ["analysis_set", "endpoint", "score_name", "n", "events", "HR_per_1SD", "HR_95CI_low", "HR_95CI_high", "p_value", "direction", "analysis_scope"],
        [[r[0], r[1], r[2], r[3], r[4], fmt(r[5]), fmt(r[6]), fmt(r[7]), fmt(r[8]), r[9], r[10]] for r in os_rows],
    )

    secondary_rows = []
    recurrence_complete = sum(1 for r in tace_rows if r["recurrence_months"] != "" and r["recurrence_status"] in ("0", "1"))
    for score_col in [PRIMARY_SCORE, "total_ecotype_unique_gene_score"] + [f"{c}_score" for c in COMPONENTS]:
        res = run_cox(tace_rows, score_col, "recurrence_months", "recurrence_status")
        secondary_rows.append(["TACE_treated", "recurrence_RFS_exploratory", score_col, res["n"], res["events"], res["hr"], res["ci_low"], res["ci_high"], res["p"], res["direction"], "secondary_exploratory_field_retained_from_Fako_table"])
    write_csv(
        "A04_TACE_treated_secondary_outcome_results.csv",
        ["analysis_set", "endpoint", "score_name", "n", "events", "HR_per_1SD", "HR_95CI_low", "HR_95CI_high", "p_value", "direction", "analysis_scope"],
        [[r[0], r[1], r[2], r[3], r[4], fmt(r[5]), fmt(r[6]), fmt(r[7]), fmt(r[8]), r[9], r[10]] for r in secondary_rows],
    )

    specificity_rows = []
    for score_col in score_cols:
        res = run_cox(resection_rows, score_col, "os_months", "os_status")
        specificity_rows.append(["resection_only", "OS", score_col, res["n"], res["events"], res["hr"], res["ci_low"], res["ci_high"], res["p"], res["direction"], "specificity_check_not_main_validation"])
    for score_col in [PRIMARY_SCORE, "total_ecotype_unique_gene_score"]:
        res = run_cox(resection_rows, score_col, "recurrence_months", "recurrence_status")
        specificity_rows.append(["resection_only", "recurrence_RFS_exploratory", score_col, res["n"], res["events"], res["hr"], res["ci_low"], res["ci_high"], res["p"], res["direction"], "specificity_check_secondary"])
    write_csv(
        "A04_resection_only_specificity_check_results.csv",
        ["analysis_set", "endpoint", "score_name", "n", "events", "HR_per_1SD", "HR_95CI_low", "HR_95CI_high", "p_value", "direction", "analysis_scope"],
        [[r[0], r[1], r[2], r[3], r[4], fmt(r[5]), fmt(r[6]), fmt(r[7]), fmt(r[8]), r[9], r[10]] for r in specificity_rows],
    )

    primary_tace = next(r for r in os_rows if r[0] == "TACE_treated" and r[2] == PRIMARY_SCORE)
    primary_resection = next(r for r in specificity_rows if r[0] == "resection_only" and r[2] == PRIMARY_SCORE and r[1] == "OS")
    verdict = "READY_FOR_A05" if primary_tace[5] != "" and primary_tace[5] > 1 else "NO_GO"

    audit = f"""# A04 GSE14520 sample flow and endpoint audit

GSE14520/Fako is used only as treatment-context validation. It is not a TACE response cohort and must not be described as response validation.

## Sample flow after score/clinical mapping

- Total mapped scored clinical samples: {len(score_dicts)}
- TACE-treated: {flow_counts['TACE_treated']}
- Resection-only: {flow_counts['resection_only']}
- Other therapy: {flow_counts['other_therapy']}
- Missing survival data group: {flow_counts['missing_survival_data']}

Therapy grouping counts: {dict(therapy_counts)}

Platform by context counts: {dict(platform_context_counts)}

## Endpoint fields

- OS fields: `os_status`, `os_months`; used as primary post-TACE outcome association in TACE-treated samples.
- Recurrence/RFS fields: `recurrence_status`, `recurrence_months`; retained as secondary/exploratory because this is reconstructed from Fako supplement fields and not treated as primary endpoint.
- TACE-treated recurrence-complete samples: {recurrence_complete}

## Boundary

- No response labels exist in GSE14520/Fako for this use.
- No outcome-guided gene pruning was performed.
- No AUC/ROC/cutoff/train-test/nomogram/calibration/decision curve was computed.
"""
    (A04 / "A04_GSE14520_sample_flow_and_endpoint_audit.md").write_text(audit, encoding="utf-8")

    summary = f"""# A04 treatment-context validation summary

A04 used the fixed A02 gene set and A03 scoring rule. GSE14520 was not used as a response validation cohort.

## Primary TACE-treated OS association

- Analysis set: TACE-treated
- N: {primary_tace[3]}
- Events: {primary_tace[4]}
- Score: `{PRIMARY_SCORE}`
- HR per 1 SD: {fmt(primary_tace[5])}
- 95% CI: {fmt(primary_tace[6])} to {fmt(primary_tace[7])}
- P value: {fmt(primary_tace[8])}
- Direction: {primary_tace[9]}

## Resection-only specificity check

- Analysis set: resection-only
- N: {primary_resection[3]}
- Events: {primary_resection[4]}
- Score: `{PRIMARY_SCORE}`
- HR per 1 SD: {fmt(primary_resection[5])}
- 95% CI: {fmt(primary_resection[6])} to {fmt(primary_resection[7])}
- P value: {fmt(primary_resection[8])}
- Direction: {primary_resection[9]}

## Interpretation boundary

This is treatment-context validation of a fixed ecotype score in post-TACE outcome, plus resection-only specificity check. It is not a response predictor, not a clinical prediction model, and not a generic OS/RFS model.

Final verdict: `{verdict}`.
"""
    (A04 / "A04_treatment_context_validation_summary.md").write_text(summary, encoding="utf-8")

    final = f"""# final_A04_verdict

`{verdict}`

## Primary A04 result

In the GSE14520/Fako TACE-treated set, the fixed A02/A03 ecotype score showed:

- Endpoint: OS
- N: {primary_tace[3]}
- Events: {primary_tace[4]}
- HR per 1 SD: {fmt(primary_tace[5])}
- 95% CI: {fmt(primary_tace[6])} to {fmt(primary_tace[7])}
- P value: {fmt(primary_tace[8])}
- Direction: {primary_tace[9]}

## Specificity check

In the resection-only set:

- Endpoint: OS
- N: {primary_resection[3]}
- Events: {primary_resection[4]}
- HR per 1 SD: {fmt(primary_resection[5])}
- 95% CI: {fmt(primary_resection[6])} to {fmt(primary_resection[7])}
- P value: {fmt(primary_resection[8])}
- Direction: {primary_resection[9]}

## Required boundary

GSE14520 is a treatment-context validation and specificity-check resource only. It must not be written as a TACE response validation cohort.

## A05-ready files

- Scores: `results/a04_gse14520_treatment_context_validation/A04_GSE14520_ecotype_scores.csv`
- Primary OS association: `results/a04_gse14520_treatment_context_validation/A04_TACE_treated_OS_association_results.csv`
- Specificity check: `results/a04_gse14520_treatment_context_validation/A04_resection_only_specificity_check_results.csv`

## Still forbidden

- No outcome-guided edits to A02 gene set.
- No LASSO/RF/SVM/AUC/ROC/cutoff/train-test/nomogram/calibration/decision curve.
- No clinical prediction model claim.
- No OS/RFS-led reframing.
- No ICI/TLS/MVI/perivascular/immune-exclusion/early-recurrence pivot.
"""
    (A04 / "final_A04_verdict.md").write_text(final, encoding="utf-8")


if __name__ == "__main__":
    main()
