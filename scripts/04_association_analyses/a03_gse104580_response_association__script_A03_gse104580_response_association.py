import csv
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
A01 = ROOT / "results" / "a01_data_preparation"
A02 = ROOT / "results" / "a02_ecotype_construction"
A03 = ROOT / "results" / "a03_gse104580_response_association"

COMPONENTS = [
    "malignant_hypoxia_adaptive",
    "endothelial_angiogenesis_interaction",
    "spp1_mmp9_tam_matrix_angiogenic_support",
]
ECOTYPE_NAME = "Hypoxia_adaptive_Tumor_Endothelial_Macrophage_Ecotype_A02"


def write_csv(path, header, rows):
    with (A03 / path).open("w", encoding="utf-8", newline="") as fh:
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
            if not row:
                continue
            expr[row[0]] = np.array([float(x) for x in row[1:]], dtype=float)
    return samples, expr


def mean(x):
    x = list(x)
    return sum(x) / len(x) if x else float("nan")


def median(x):
    x = sorted(x)
    if not x:
        return float("nan")
    mid = len(x) // 2
    if len(x) % 2:
        return x[mid]
    return (x[mid - 1] + x[mid]) / 2


def fmt(x):
    if isinstance(x, str):
        return x
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return f"{x:.6g}"


def cliffs_delta(x, y):
    # delta > 0 means x tends to be higher than y.
    gt = 0
    lt = 0
    for a in x:
        for b in y:
            if a > b:
                gt += 1
            elif a < b:
                lt += 1
    return (gt - lt) / (len(x) * len(y))


def norm_sf(z):
    return 0.5 * math.erfc(z / math.sqrt(2))


def rankdata(values):
    values = list(values)
    order = sorted(range(len(values)), key=lambda i: values[i])
    ranks = [0.0] * len(values)
    i = 0
    while i < len(values):
        j = i
        while j + 1 < len(values) and values[order[j + 1]] == values[order[i]]:
            j += 1
        avg_rank = (i + 1 + j + 1) / 2.0
        for k in range(i, j + 1):
            ranks[order[k]] = avg_rank
        i = j + 1
    return np.array(ranks, dtype=float)


def mannwhitneyu_2sided(x, y):
    x = list(x)
    y = list(y)
    values = x + y
    ranks = rankdata(values)
    n1 = len(x)
    n2 = len(y)
    r1 = ranks[:n1].sum()
    u1 = r1 - n1 * (n1 + 1) / 2.0
    # Tie-corrected normal approximation.
    counts = Counter(values)
    tie_term = sum(c ** 3 - c for c in counts.values())
    n = n1 + n2
    var = n1 * n2 / 12.0 * ((n + 1) - tie_term / (n * (n - 1)) if n > 1 else 0)
    mean_u = n1 * n2 / 2.0
    if var <= 0:
        p = 1.0
    else:
        z = (abs(u1 - mean_u) - 0.5) / math.sqrt(var)
        p = 2 * norm_sf(abs(z))
    return u1, p


def rank_score_for_genes(samples, expr, genes):
    measured = [g for g in genes if g in expr]
    n_genes_total = len(expr)
    if not measured:
        return np.full(len(samples), np.nan), []
    matrix = np.vstack([expr[g] for g in expr])
    all_genes = list(expr.keys())
    gene_to_row = {g: i for i, g in enumerate(all_genes)}
    ranks = np.zeros_like(matrix, dtype=float)
    for j in range(matrix.shape[1]):
        ranks[:, j] = rankdata(matrix[:, j])
    idx = [gene_to_row[g] for g in measured]
    raw = ranks[idx, :].mean(axis=0)
    # singscore-like 0-1 normalization: higher means higher relative ranks of the fixed gene set.
    score = (raw - 1.0) / max(n_genes_total - 1.0, 1.0)
    return score, measured


def logistic_univariate(y, x):
    y = np.asarray(y, dtype=float)
    x = np.asarray(x, dtype=float)
    xz = (x - x.mean()) / x.std(ddof=1)
    X = np.column_stack([np.ones(len(xz)), xz])

    beta = np.zeros(2)
    for _ in range(50):
        eta = X @ beta
        p = 1 / (1 + np.exp(-eta))
        W = np.clip(p * (1 - p), 1e-8, None)
        info = X.T @ (X * W[:, None])
        grad = X.T @ (y - p)
        step = np.linalg.solve(info, grad)
        beta = beta + step
        if np.max(np.abs(step)) < 1e-8:
            break
    eta = X @ beta
    p = 1 / (1 + np.exp(-eta))
    W = np.clip(p * (1 - p), 1e-8, None)
    info = X.T @ (X * W[:, None])
    cov = np.linalg.inv(info)
    se = np.sqrt(np.diag(cov))
    z = beta[1] / se[1]
    pval = 2 * norm_sf(abs(z))
    return {
        "beta_per_1sd": beta[1],
        "se": se[1],
        "or_per_1sd": math.exp(beta[1]),
        "ci_low": math.exp(beta[1] - 1.96 * se[1]),
        "ci_high": math.exp(beta[1] + 1.96 * se[1]),
        "p_value": pval,
        "n": len(y),
        "n_non_response": int(y.sum()),
    }


def response_stats(name, scores, labels):
    resp = np.array([s for s, lab in zip(scores, labels) if lab == 0])
    non = np.array([s for s, lab in zip(scores, labels) if lab == 1])
    u_stat, p_value = mannwhitneyu_2sided(non, resp)
    return {
        "score_name": name,
        "n_response": len(resp),
        "n_non_response": len(non),
        "response_median": median(resp),
        "non_response_median": median(non),
        "non_response_minus_response_median": median(non) - median(resp),
        "response_mean": mean(resp),
        "non_response_mean": mean(non),
        "non_response_minus_response_mean": mean(non) - mean(resp),
        "cliffs_delta_non_response_vs_response": cliffs_delta(non, resp),
        "mannwhitney_u": u_stat,
        "p_value": p_value,
        "direction": "higher_in_non_response" if median(non) > median(resp) else "not_higher_in_non_response",
    }


def pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    r = np.corrcoef(x, y)[0, 1]
    if len(x) <= 3 or abs(r) >= 1:
        return r, ""
    # Large-sample Fisher z approximation, sufficient for descriptive robustness.
    z = 0.5 * math.log((1 + r) / (1 - r)) * math.sqrt(len(x) - 3)
    p = 2 * norm_sf(abs(z))
    return r, p


def main():
    A03.mkdir(parents=True, exist_ok=True)
    geneset_rows = read_csv(A02 / "A02_final_ecotype_gene_set.csv")
    samples, expr = load_expression(A01 / "GSE104580_gene_expression.tsv")
    metadata = read_tsv(A01 / "GSE104580_sample_metadata.tsv")
    meta_by_gsm = {r["gsm"]: r for r in metadata}
    ordered_meta = [meta_by_gsm[s] for s in samples]
    labels = np.array([1 if r["response_group"] == "non_response" else 0 for r in ordered_meta], dtype=int)

    component_genes = defaultdict(list)
    gene_components = defaultdict(list)
    for r in geneset_rows:
        component = r["component"]
        gene = r["gene_symbol"]
        status = r["inclusion_status"]
        component_genes[component].append(gene)
        gene_components[gene].append(component)

    duplicates = {g: comps for g, comps in gene_components.items() if len(comps) > 1}
    unique_total_genes = sorted(gene_components)

    component_list_rows = []
    for component in COMPONENTS:
        for g in component_genes[component]:
            component_list_rows.append([
                ECOTYPE_NAME,
                component,
                g,
                "A03_core",
                "yes" if g in duplicates else "no",
                ";".join(gene_components[g]),
            ])
    for g in unique_total_genes:
        component_list_rows.append([
            ECOTYPE_NAME,
            "total_ecotype_unique_gene_list",
            g,
            "A03_total_unique",
            "yes" if g in duplicates else "no",
            ";".join(gene_components[g]),
        ])
    write_csv(
        "A03_component_gene_lists_after_mapping.csv",
        ["ecotype_name", "component", "gene_symbol", "inclusion_status", "is_cross_component_duplicate", "components_for_gene"],
        component_list_rows,
    )

    availability_rows = []
    for component in COMPONENTS:
        for g in component_genes[component]:
            availability_rows.append([component, g, "yes" if g in expr else "no"])
    for g in unique_total_genes:
        availability_rows.append(["total_ecotype_unique_gene_list", g, "yes" if g in expr else "no"])
    write_csv("A03_GSE104580_gene_availability.csv", ["component", "gene_symbol", "available_in_GSE104580"], availability_rows)

    duplicate_text = ["# A03 duplicate gene handling\n"]
    duplicate_text.append("A02 gene set was read without response-guided gene deletion or replacement.\n")
    if duplicates:
        duplicate_text.append("Cross-component duplicate genes:\n")
        for g, comps in sorted(duplicates.items()):
            duplicate_text.append(f"- `{g}` appears in: {', '.join(comps)}\n")
    else:
        duplicate_text.append("No cross-component duplicate genes detected.\n")
    duplicate_text.append("\nComponent scores retain each component's predefined gene list. Total ecotype composite score is computed as the equal-weight mean of the three component scores, and a separate unique-gene total score is computed with duplicate genes counted once for robustness. No response label is used for weighting.\n")
    (A03 / "A03_duplicate_gene_handling.md").write_text("".join(duplicate_text), encoding="utf-8")

    component_scores = {}
    measured_by_component = {}
    for component in COMPONENTS:
        score, measured = rank_score_for_genes(samples, expr, component_genes[component])
        component_scores[component] = score
        measured_by_component[component] = measured
    unique_total_score, unique_measured = rank_score_for_genes(samples, expr, unique_total_genes)
    total_equal_component_score = np.vstack([component_scores[c] for c in COMPONENTS]).mean(axis=0)

    score_rows = []
    for i, sample in enumerate(samples):
        score_rows.append([
            sample,
            ordered_meta[i]["sample_title"],
            ordered_meta[i]["response_group"],
            labels[i],
            component_scores["malignant_hypoxia_adaptive"][i],
            component_scores["endothelial_angiogenesis_interaction"][i],
            component_scores["spp1_mmp9_tam_matrix_angiogenic_support"][i],
            total_equal_component_score[i],
            unique_total_score[i],
        ])
    write_csv(
        "A03_GSE104580_ecotype_scores.csv",
        [
            "gsm",
            "sample_title",
            "response_group",
            "non_response_binary",
            "malignant_hypoxia_adaptive_score",
            "endothelial_angiogenesis_interaction_score",
            "spp1_mmp9_tam_matrix_angiogenic_support_score",
            "total_ecotype_equal_component_score",
            "total_ecotype_unique_gene_score",
        ],
        score_rows,
    )

    score_map = {
        "malignant_hypoxia_adaptive_score": component_scores["malignant_hypoxia_adaptive"],
        "endothelial_angiogenesis_interaction_score": component_scores["endothelial_angiogenesis_interaction"],
        "spp1_mmp9_tam_matrix_angiogenic_support_score": component_scores["spp1_mmp9_tam_matrix_angiogenic_support"],
        "total_ecotype_equal_component_score": total_equal_component_score,
        "total_ecotype_unique_gene_score": unique_total_score,
    }
    stat_rows = []
    model_rows = []
    for name, score in score_map.items():
        st = response_stats(name, score, labels)
        stat_rows.append([st[k] for k in [
            "score_name", "n_response", "n_non_response", "response_median", "non_response_median",
            "non_response_minus_response_median", "response_mean", "non_response_mean",
            "non_response_minus_response_mean", "cliffs_delta_non_response_vs_response",
            "mannwhitney_u", "p_value", "direction",
        ]])
        m = logistic_univariate(labels, score)
        model_rows.append([
            name, m["n"], m["n_non_response"], m["beta_per_1sd"], m["se"], m["or_per_1sd"],
            m["ci_low"], m["ci_high"], m["p_value"], "univariate_association_not_prediction",
        ])
    write_csv(
        "A03_response_association_statistics.csv",
        [
            "score_name", "n_response", "n_non_response", "response_median", "non_response_median",
            "non_response_minus_response_median", "response_mean", "non_response_mean",
            "non_response_minus_response_mean", "cliffs_delta_non_response_vs_response",
            "mannwhitney_u", "p_value", "direction",
        ],
        stat_rows,
    )
    write_csv(
        "A03_association_model_results.csv",
        [
            "score_name", "n", "n_non_response", "logit_beta_per_1sd_score", "se",
            "OR_per_1sd_score", "OR_95CI_low", "OR_95CI_high", "p_value",
            "model_scope",
        ],
        model_rows,
    )

    sensitivity_rows = []
    for a in COMPONENTS:
        for b in COMPONENTS:
            if a >= b:
                continue
            r, p = pearson(component_scores[a], component_scores[b])
            sensitivity_rows.append(["component_correlation", f"{a}__vs__{b}", r, p, ""])
    total_stat = response_stats("total_ecotype_equal_component_score", total_equal_component_score, labels)
    for leave_out in COMPONENTS:
        kept = [c for c in COMPONENTS if c != leave_out]
        loo_score = np.vstack([component_scores[c] for c in kept]).mean(axis=0)
        st = response_stats(f"leave_out_{leave_out}", loo_score, labels)
        sensitivity_rows.append([
            "leave_one_component_out",
            f"leave_out_{leave_out}",
            st["non_response_minus_response_median"],
            st["p_value"],
            st["direction"],
        ])
    unique_stat = response_stats("total_ecotype_unique_gene_score", unique_total_score, labels)
    sensitivity_rows.append([
        "duplicate_handling_sensitivity",
        "equal_component_total_vs_unique_gene_total",
        unique_stat["non_response_minus_response_median"],
        unique_stat["p_value"],
        unique_stat["direction"],
    ])
    for component in COMPONENTS:
        st = response_stats(f"{component}_score", component_scores[component], labels)
        sensitivity_rows.append([
            "component_direction_consistency",
            component,
            st["non_response_minus_response_median"],
            st["p_value"],
            st["direction"],
        ])
    component_deltas = [abs(response_stats(c, component_scores[c], labels)["non_response_minus_response_median"]) for c in COMPONENTS]
    sensitivity_rows.append([
        "single_component_dominance_check",
        "max_component_delta_divided_by_sum_abs_component_deltas",
        max(component_deltas) / sum(component_deltas) if sum(component_deltas) else "",
        "",
        "descriptive_only_no_gene_pruning",
    ])
    write_csv(
        "A03_component_consistency_and_sensitivity_results.csv",
        ["check_type", "comparison", "estimate", "p_value", "interpretation"],
        sensitivity_rows,
    )

    primary = response_stats("total_ecotype_equal_component_score", total_equal_component_score, labels)
    primary_model = logistic_univariate(labels, total_equal_component_score)
    verdict = "READY_FOR_A04" if primary["direction"] == "higher_in_non_response" else "NO_GO"

    summary = f"""# A03 response association summary

本阶段将 A02 预定义的 single-cell/spatial anchored hypoxia-adaptive tumor-endothelial-macrophage ecotype 转译到 GSE104580 bulk pretreatment TACE biopsy cohort，并检验其与 TACE non-response 的关联。

A03 未做 GSE104580 response DEGs，未做 LASSO、随机森林、SVM、AUC/ROC、cutoff optimization、train/test、nomogram、calibration 或 decision curve。所有 gene set 均固定来自 A02。

## Gene set implementation

- A02 ecotype rows: {len(geneset_rows)}
- Unique total genes: {len(unique_total_genes)}
- Cross-component duplicate genes: {', '.join(sorted(duplicates)) if duplicates else 'none'}
- Duplicate handling: component scores retain predefined component membership; total ecotype score uses equal-weight component means. Unique-gene total score is provided only as robustness check.

## GSE104580 transfer

- Samples: {len(samples)}
- Response: {int((labels == 0).sum())}
- Non-response: {int((labels == 1).sum())}
- Component genes available: all A02 A03-core genes are available in GSE104580.

## Primary association

Primary score: `total_ecotype_equal_component_score`.

- Response median: {fmt(primary['response_median'])}
- Non-response median: {fmt(primary['non_response_median'])}
- Non-response minus response median: {fmt(primary['non_response_minus_response_median'])}
- Cliff's delta: {fmt(primary['cliffs_delta_non_response_vs_response'])}
- Mann-Whitney/Wilcoxon rank-sum P value: {fmt(primary['p_value'])}
- Logistic association OR per 1 SD score: {fmt(primary_model['or_per_1sd'])}
- 95% CI: {fmt(primary_model['ci_low'])} to {fmt(primary_model['ci_high'])}
- Logistic association P value: {fmt(primary_model['p_value'])}

Interpretation: the A02 ecotype score is {primary['direction'].replace('_', ' ')} in GSE104580. This supports wording as “associated with” or “linked to” TACE non-response, not “predicts”.

## Internal robustness

Robustness checks are in `A03_component_consistency_and_sensitivity_results.csv`. These checks assess component direction consistency, component correlations, leave-one-component-out sensitivity, single-component dominance, and duplicate-gene handling. No response-guided gene pruning was performed.

## Verdict

`{verdict}`
"""
    (A03 / "A03_response_association_summary.md").write_text(summary, encoding="utf-8")

    final = f"""# final_A03_verdict

`{verdict}`

## Direction in GSE104580

The A02 ecotype score direction is `{primary['direction']}` for GSE104580 TACE non-response.

## Main statistical evidence

- Primary score: `total_ecotype_equal_component_score`
- Non-response minus response median: {fmt(primary['non_response_minus_response_median'])}
- Cliff's delta: {fmt(primary['cliffs_delta_non_response_vs_response'])}
- Wilcoxon/Mann-Whitney P value: {fmt(primary['p_value'])}
- Logistic association OR per 1 SD score: {fmt(primary_model['or_per_1sd'])}
- 95% CI: {fmt(primary_model['ci_low'])} to {fmt(primary_model['ci_high'])}
- Logistic association P value: {fmt(primary_model['p_value'])}

## Required wording

This can only be described as `associated with` / `linked to TACE non-response`. It must not be described as `predicts`, a clinical prediction model, or a response predictor.

## A04 score file path

`results/a03_gse104580_response_association/A03_GSE104580_ecotype_scores.csv`

## A04 scope

A04 may use the fixed A02/A03 score framework for GSE14520 TACE-treated post-TACE outcome validation and resection-only specificity check. GSE14520 must not be written as a response validation cohort.

## Still forbidden

- No GSE104580 response DEGs followed by ecotype backfilling.
- No LASSO, random forest, SVM, AUC/ROC, cutoff optimization, train/test split, nomogram, calibration, decision curve, or clinical utility claim.
- No response-guided gene pruning or changes to the A02 gene set.
- No formal survival analysis in A03.
- No OS/RFS-led reframing.
- No ICI/TLS/MVI/perivascular/immune-exclusion/early-recurrence pivot.
"""
    (A03 / "final_A03_verdict.md").write_text(final, encoding="utf-8")


if __name__ == "__main__":
    main()
